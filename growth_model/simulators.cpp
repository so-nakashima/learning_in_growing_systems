#include "headers/simulators.h"
#include <cassert>
#include <algorithm>
#include <iterator>

MBPRE::MBPRE(){
}

MBPRE::~MBPRE(){

}

void MBPRE::excecute(){
    //check offspring_distributions is valid
    assert(offspring_distributions.size() == env->cardinality());
    for(int y = 0; y != env->cardinality(); y++){
        assert(offspring_distributions[y].size() == pop->cardinality());
        for(int x = 0; x != pop->cardinality(); x++){
            assert(!offspring_distributions[y][x].empty());
            for(auto z : offspring_distributions[y][x]){
                assert(z >= 0);
            }
        }
    }


    //run simulation
    for(int t = 0; t != end_time; t++){
        //record (this is first in order to record t = 0)
        record();
        if(t == 0)
            pop->record(out_pop_full);

        time_evolution();

    }
    record();
}

void MBPRE::record(){
    env->record(out_env);
    pop->record(out_pop);
}

void MBPRE::time_evolution(){
    //population
    pop->time_evolution(offspring_distributions[env->current_state()]);
    pop->record(out_pop_full); //record to-be-discored cells in advance of selection.
    pop->selection();

    //environment
    env->next_state(time);
    time++;
}



void MBPRE::set_population(Population* population){
    pop = population;
}

void MBPRE::set_offspring_distributions(const std::vector<std::vector<std::vector<double>>>& offspring_dist){
    offspring_distributions = offspring_dist;
}


Markov_Environments::~Markov_Environments(){
}

void Markov_Environments::set_initial_distribution(const std::vector<double>& init){
    assert(init.size() == m_cardinality);
    for(int i = 0; i != m_cardinality; i++){
        assert(init[i] >= 0);
    }

    std::discrete_distribution<int> dist(init.begin(), init.end());
    m_current_state = dist(mt);
}

void Markov_Environments::set_cardinality(int n){
    assert(n > 0);
    m_cardinality = n;
}

void Markov_Environments::set_transition(const std::vector<std::vector<double>>& tran_mat){
    assert(tran_mat.size() == m_cardinality);

    std::vector<std::vector<std::function<double(int)>>> transit_func_mat;

    for(int i = 0; i != m_cardinality; i++){
        assert(tran_mat[i].size() == m_cardinality);
        std::vector<std::function<double(int)>> temp;
        for(int j = 0; j != m_cardinality; j++){
            assert(tran_mat[i][j] >= 0);
            temp.push_back([=](int _){return tran_mat[i][j];});
        }
        transit_func_mat.push_back(temp);
    }

    transition = transit_func_mat;
}

void Markov_Environments::set_transition(const std::vector<std::vector<std::function<double(int)>>>& tran_func_mat){
    assert(tran_func_mat.size() == m_cardinality);
    for(int i = 0; i != m_cardinality; i++){
        assert(tran_func_mat[i].size() == m_cardinality);
    }

    transition = tran_func_mat;
}



int Markov_Environments::next_state(int t){

    std::vector<double> distribution_vec;
    for(int i = 0; i != m_cardinality; i++){
        distribution_vec.push_back(std::max(transition[current_state()][i](t), 0.0));
    }

    std::discrete_distribution<int> dist(distribution_vec.begin(), distribution_vec.end());

    m_current_state = dist(mt);
    return m_current_state;
}



Markov_Environments::Markov_Environments(const std::vector<std::vector<double>>& transit_mat ,const int cardinality , const std::vector<double>& initial){
    //set random generator
    std::random_device rnd;
    mt.seed(rnd());

    set_cardinality(cardinality);

    if(!initial.empty()){
        set_initial_distribution(initial);
    }
    else{
        m_current_state = 0;
    }

    //set transition matrix if valid; otherwise initialize by uniform transition
    if(transit_mat.size() != 0){
        set_transition(transit_mat);
    }
    else
    {
        set_transition(std::vector<std::vector<double>>(cardinality, std::vector<double>(cardinality, 1.0)));
    }
    
}

 Markov_Environments::Markov_Environments(const std::vector<std::vector<std::function<double(int)>>>& transit_func_mat ,const int cardinality , const std::vector<double>& initial){
    //set random generator
    std::random_device rnd;
    mt.seed(rnd());

    set_cardinality(cardinality);

    if(!initial.empty()){
        set_initial_distribution(initial);
    }
    else{
        m_current_state = 0;
    }

    //set transition matrix if valid; otherwise initialize by uniform transition
    set_transition(transit_func_mat);
} 

void Markov_Environments::record(std::ofstream* out){
    *out << std::to_string(current_state()) << std::endl;
}


Cell::Cell(int type , std::string ID)
{
    //validity of argument type is assured in Cells class (i.e. daughters method)
    m_type = type;
    m_ID = ID;
}

Cell::~Cell()
{
}

std::vector<Cell> Cell::daughters(const std::vector<std::vector<double>>& type_transition, const std::vector<std::vector<double>>& offspring_distribution, std::mt19937_64& mt) const {

    int no_offsprings = std::discrete_distribution<int>(offspring_distribution[m_type].begin(), offspring_distribution[m_type].end())(mt);

    std::vector<Cell> res;
    for(int i = 0; i != no_offsprings; i++){
        int next_type = std::discrete_distribution<int>(type_transition[m_type].begin(), type_transition[m_type].end())(mt);
        res.push_back(Cell(next_type, m_ID + std::to_string(i)));
    }

    return res;
}

void Cell::record(std::ofstream* out_population){
    *out_population << m_ID << " " << m_type << std::endl;
}

Cells::Cells(int type_no, 
    const std::vector<std::vector<double>>& type_tran, 
    const std::vector<Cell>& initial_population,int max_pop_size)
{
    set_maximum_population_size(max_pop_size);
    set_type_cardinality(type_no);
    set_initial_population(initial_population);

    if(type_tran.empty()){//default argument (uniform transition)
        set_type_transition(std::vector<std::vector<double>>(type_cardinality, std::vector<double>(type_cardinality, 1.0)));
    }
    else
    {
        set_type_transition(type_tran);
    }

    //random generator
    std::random_device rnd;
    mt.seed(rnd());
}

Cells::~Cells()
{
}

void Cells::set_type_cardinality(int n){
    assert(n > 0);
    type_cardinality = n;
}

void Cells::set_type_transition(const std::vector<std::vector<double>>& transit){
    assert(transit.size() == type_cardinality);
    for(int i = 0; i != type_cardinality; i++){
        assert(transit[i].size() == type_cardinality);
    }

    type_transition = transit;
}


void Cells::set_maximum_population_size(int n){
    assert(n > 0);
    maximum_population_size = n;
}

void Cells::set_initial_population(const std::vector<Cell>& initial_population){
    for(auto cell : initial_population){
        assert(cell.type() < type_cardinality);
    }
    current_population = initial_population;
}

void Cells::time_evolution(const std::vector<std::vector<double>>& offspring_distribution){
    std::vector<Cell> new_population;
    for(auto cell : current_population){
        std::vector<Cell> daughters = cell.daughters(type_transition, offspring_distribution, mt);

        new_population.insert(new_population.end(), daughters.begin(), daughters.end());
    }

    current_population.clear();
    current_population = new_population;
}

void Cells::selection(){
    if(current_population.size() > maximum_population_size){
        std::vector<Cell> temp;

        std::sample(current_population.begin(), current_population.end(), std::back_inserter(temp), maximum_population_size, mt);
        current_population = temp;
    }
}

void Cells::record(std::ofstream* out_population){
    for (auto cell : current_population){
        cell.record(out_population);
    }
}


std::vector<Cell_Learn> const Cell_Learn::daughters(const std::vector<std::vector<double>>& offspring_distribution,
    const std::function<void(int, int, int,  std::vector<std::vector<double>>&, std::vector<std::vector<double>>&, std::vector<double>&,std::vector<double>&, std::mt19937_64&)>& learning_rule,
    std::mt19937_64& mt){

    int no_offsprings = std::discrete_distribution<int>(offspring_distribution[m_type].begin(), offspring_distribution[m_type].end())(mt);

    std::vector<Cell_Learn> res;
    for(int i = 0; i != no_offsprings; i++){
        int next_type = std::discrete_distribution<int>(transition[m_type].begin(), transition[m_type].end())(mt);


        //learning
        std::vector<std::vector<double>> temp_jump = ancestral_jump;
        std::vector<std::vector<double>> temp_tran_mat = transition;
        std::vector<double> temp_rep_hist = replication_history;
        std::vector<double> temp_mem = memory;

        learning_rule(type(), next_type, no_offsprings, temp_tran_mat, temp_jump, temp_rep_hist, temp_mem, mt);

        res.push_back(
            Cell_Learn(next_type, m_ID + std::to_string(i),
            temp_jump,
            temp_tran_mat,
            temp_rep_hist,
            temp_mem
        ));
    }

    return res;
}

Cell_Learn::Cell_Learn(int type, std::string ID, 
    const std::vector<std::vector<double>>& jump, 
    const std::vector<std::vector<double>>& tran_mat,
    const std::vector<double>& rep_hist,
    const std::vector<double>& mem){
    m_type = type;
    m_ID = ID;
    set_ancestral_jump(jump);
    set_transition(tran_mat);
    set_memory(mem);

}

void Cell_Learn::record(std::ofstream* out){
    *out << id() << " " << type() << std::endl;
    out_mat(ancestral_jump, out);
    out_mat(transition, out);
    //out_vec(replication_history, out);
    out_vec(memory, out);
    *out << std::endl;
}

Cells_Learn::Cells_Learn(int type_no, int max_pop_size, 
const std::vector<Cell_Learn>& initial_pop,
const std::function<void(int, int, int, std::vector<std::vector<double>>&, std::vector<std::vector<double>>&, std::vector<double>&,std::vector<double>&, std::mt19937_64&)>& rule){

    set_type_cardinality(type_no);
    set_maximum_population_size(max_pop_size);
    set_learning_rule(rule);
    set_initial_population(initial_pop);

    //random generator
    std::random_device rnd;
    mt.seed(rnd());
}


void Cells_Learn::selection(){
    if(current_population.size() > maximum_population_size){
        std::vector<Cell_Learn> temp;

        std::sample(current_population.begin(), current_population.end(), std::back_inserter(temp), maximum_population_size, mt);
        current_population = temp;
    }
}

void Cells_Learn::time_evolution(const std::vector<std::vector<double>>& offspring_distribution){
    std::vector<Cell_Learn> new_population;
    for(auto cell : current_population){
        std::vector<Cell_Learn> daughters = cell.daughters(offspring_distribution, learning_rule, mt);

        new_population.insert(new_population.end(), daughters.begin(), daughters.end());
    }

    current_population.clear();
    current_population = new_population;
}

void Cells_Learn::record(std::ofstream* out_population){
    for (auto cell : current_population){
        cell.record(out_population);
    }
}