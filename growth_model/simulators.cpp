#include "headers/simulators.h"
#include <cassert>

World_MBPRE::World_MBPRE(){
}

World_MBPRE::~World_MBPRE(){

}

void World_MBPRE::excecute(){
    for(int t = 0; t != end_time; t++){
        //record (this is first in order to record t = 0)
        record();

        time_evolution();

    }
    record();
}
void World_MBPRE::record(){
    env.record(out_env);
    pop->record(out_pop);
}

void World_MBPRE::time_evolution(){
    //environments
    int nowstate = env.next_state();

    //population
    pop->time_evolution(nowstate);
}



void World_MBPRE::set_population(Population* population){
    pop = population;
}

Environments::~Environments(){
}

void Environments::set_initial_distribution(const std::vector<double>& init){
    assert(init.size() == cardinality);
    for(int i = 0; i != cardinality; i++){
        assert(init[i] >= 0);
    }

    std::discrete_distribution<int> dist(init.begin(), init.end());
    m_current_state = dist(mt);
}

void Environments::set_cardinality(int n){
    assert(n > 0);
    cardinality = n;
}

void Environments::set_transition(const std::vector<std::vector<double>>& tran_mat){
    assert(tran_mat.size () == cardinality);
    for(int i = 0; i != cardinality; i++){
        assert(tran_mat[i].size() == cardinality);
        for(int j = 0; j != cardinality; j++){
            assert(tran_mat[i][j] >= 0);
        }
    }
        transition = tran_mat;
}



int Environments::next_state(){
    std::discrete_distribution<int> dist(transition[current_state()].begin(), transition[current_state()].end());

    m_current_state = dist(mt);
    return m_current_state;
}



Environments::Environments(const int cardinality , const std::vector<double>& initial  , const std::vector<std::vector<double>>& transit_mat ){
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
    if(transit_mat.size() != 0)
        transition = transit_mat;
    else
    {
        transition = std::vector<std::vector<double>>(cardinality, std::vector<double>(cardinality, 1.0));
    }
    
}

void Environments::record(std::ofstream* out){
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

std::vector<Cell> Cell::daughters(int current_env, const std::vector<std::vector<double>>& type_transition, const std::vector<std::vector<double>>& offspring_distribution, std::mt19937_64& mt){

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
    const std::vector<std::vector<double>>& offsprings_prob,
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

    if(offsprings_prob.empty()){
        set_offspring_distribution(std::vector<std::vector<double>>(type_cardinality, std::vector<double>(1,0.0)));
    }
    else
    {
        set_offspring_distribution(offsprings_prob);
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

void Cells::set_offspring_distribution(const std::vector<std::vector<double>>& offspring_dist){
    assert(offspring_dist.size() == type_cardinality);
    for(int i = 0; i != type_cardinality; i++){
        assert(0 < offspring_dist[i].size() && offspring_dist[i].size() <= type_cardinality);
    }

    offspring_distribution = offspring_dist;
}

void Cells::set_maximum_population_size(int n){
    assert(n > 0);
    maximum_population_size = n;
}

void Cells::set_initial_population(const std::vector<Cell>& initial_population){
    current_population = initial_population;
}

void Cells::time_evolution(int current_env){
    std::vector<Cell> new_population;
    for(auto cell : current_population){
        std::vector<Cell> daughters = cell.daughters(current_env, type_transition, offspring_distribution, mt);

        new_population.insert(new_population.end(), daughters.begin(), daughters.end());
    }

    current_population = new_population;

    if(current_population.size() > maximum_population_size){
        //to be implemented
    }
}

void Cells::record(std::ofstream* out_population){
    for (auto cell : current_population){
        cell.record(out_population);
    }
}
