#include "headers/simulators.h"
#include <cassert>
#include <algorithm>
#include <iterator>

MBPRE::MBPRE()
{
}

MBPRE::~MBPRE()
{
}

void MBPRE::excecute()
{
    //check offspring_distributions is valid
    assert(offspring_distributions.size() == env->cardinality());
    for (int y = 0; y != env->cardinality(); y++)
    {
        assert(offspring_distributions[y].size() == pop->cardinality());
        for (int x = 0; x != pop->cardinality(); x++)
        {
            assert(!offspring_distributions[y][x].empty());
            for (auto z : offspring_distributions[y][x])
            {
                assert(z >= 0);
            }
        }
    }

    //run simulation
    for (int t = 0; t != end_time; t++)
    {
        //record (this is first in order to record t = 0)
        record();
        //if(t == 0)
        // pop->record_before_selection();

        time_evolution();
    }
    record();
}

void MBPRE::record()
{
    env->record();
    pop->record();
}

void MBPRE::time_evolution()
{
    //environment
    int current_state = env->current_state();
    int next_state = env->next_state(time);

    //population
    pop->time_evolution(offspring_distributions[current_state], offspring_distributions[next_state]);
    //pop->record_before_selection(); //record to-be-discored cells in advance of selection.
    //pop->selection();

    time++;
}

void MBPRE::set_population(Population *population)
{
    pop = population;
}

void MBPRE::set_offspring_distributions(const std::vector<std::vector<std::vector<double>>> &offspring_dist)
{
    offspring_distributions = offspring_dist;
}

Markov_Environments::~Markov_Environments()
{
}

void Markov_Environments::set_initial_distribution(const std::vector<double> &init)
{
    assert(init.size() == m_cardinality);
    for (int i = 0; i != m_cardinality; i++)
    {
        assert(init[i] >= 0);
    }

    std::discrete_distribution<int> dist(init.begin(), init.end());
    m_current_state = dist(mt);
}

void Markov_Environments::set_cardinality(int n)
{
    assert(n > 0);
    m_cardinality = n;
}

void Markov_Environments::set_transition(const std::vector<std::vector<double>> &tran_mat)
{
    assert(tran_mat.size() == m_cardinality);

    std::vector<std::vector<std::function<double(int)>>> transit_func_mat;

    for (int i = 0; i != m_cardinality; i++)
    {
        assert(tran_mat[i].size() == m_cardinality);
        std::vector<std::function<double(int)>> temp;
        for (int j = 0; j != m_cardinality; j++)
        {
            assert(tran_mat[i][j] >= 0);
            temp.push_back([=](int _) { return tran_mat[i][j]; });
        }
        transit_func_mat.push_back(temp);
    }

    transition = transit_func_mat;
}

void Markov_Environments::set_transition(const std::vector<std::vector<std::function<double(int)>>> &tran_func_mat)
{
    assert(tran_func_mat.size() == m_cardinality);
    for (int i = 0; i != m_cardinality; i++)
    {
        assert(tran_func_mat[i].size() == m_cardinality);
    }

    transition = tran_func_mat;
}

int Markov_Environments::next_state(int t)
{

    std::vector<double> distribution_vec;
    for (int i = 0; i != m_cardinality; i++)
    {
        distribution_vec.push_back(std::max(transition[current_state()][i](t), 0.0));
    }

    std::discrete_distribution<int> dist(distribution_vec.begin(), distribution_vec.end());

    m_current_state = dist(mt);
    return m_current_state;
}

Markov_Environments::Markov_Environments(const std::vector<std::vector<double>> &transit_mat, const int cardinality, const std::vector<double> &initial)
{
    //set random generator
    std::random_device rnd;
    mt.seed(rnd());

    set_cardinality(cardinality);

    if (!initial.empty())
    {
        set_initial_distribution(initial);
    }
    else
    {
        m_current_state = 0;
    }

    //set transition matrix if valid; otherwise initialize by uniform transition
    if (transit_mat.size() != 0)
    {
        set_transition(transit_mat);
    }
    else
    {
        set_transition(std::vector<std::vector<double>>(cardinality, std::vector<double>(cardinality, 1.0)));
    }
}

Markov_Environments::Markov_Environments(const std::vector<std::vector<std::function<double(int)>>> &transit_func_mat, const int cardinality, const std::vector<double> &initial)
{
    //set random generator
    std::random_device rnd;
    mt.seed(rnd());

    set_cardinality(cardinality);

    if (!initial.empty())
    {
        set_initial_distribution(initial);
    }
    else
    {
        m_current_state = 0;
    }

    //set transition matrix if valid; otherwise initialize by uniform transition
    set_transition(transit_func_mat);
}

std::vector<int> Markov_Environments::generate(const int n)
{
    std::vector<int> res;
    for (int i = 0; i != n; i++)
    {
        res.push_back(current_state());
        next_state(i);
    }

    return res;
}

void Markov_Environments::record()
{
    if (out_env == nullptr)
        return;
    *out_env << std::to_string(current_state()) << std::endl;
}

Cell::Cell(int type, std::string ID)
{
    //validity of argument type is assured in Cells class (i.e. daughters method)
    m_type = type;
    m_ID = ID;
}

Cell::~Cell()
{
}

std::vector<Cell> Cell::daughters(const std::vector<std::vector<double>> &type_transition, const std::vector<std::vector<double>> &offspring_distribution, std::mt19937_64 &mt) const
{
    return daughters(m_type, type_transition, offspring_distribution, mt);
}

std::vector<Cell> Cell::daughters(const int common_parent_type, const std::vector<std::vector<double>> &type_transition, const std::vector<std::vector<double>> &offspring_distribution, std::mt19937_64 &mt) const
{

    const int p_type = common_parent_type; //for simplification

    //growth first
    int no_offsprings = std::discrete_distribution<int>(offspring_distribution[type()].begin(), offspring_distribution[type()].end())(mt);

    //then type-swithching
    std::vector<Cell> res;
    for (int i = 0; i != no_offsprings; i++)
    {
        int next_type = std::discrete_distribution<int>(type_transition[p_type].begin(), type_transition[p_type].end())(mt);
        res.push_back(Cell(next_type, m_ID + "S" + std::to_string(i)));
    }

    return res;
}

void Cell::record(std::ofstream *out_population)
{
    *out_population << m_ID << " " << m_type << std::endl;
}

Cells::Cells(int type_no,
             const std::vector<std::vector<double>> &type_tran,
             const std::vector<Cell> &initial_population, int max_pop_size)
{
    set_maximum_population_size(max_pop_size);
    set_type_cardinality(type_no);
    set_initial_population(initial_population);

    if (type_tran.empty())
    { //default argument (uniform transition)
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

void Cells::set_type_cardinality(int n)
{
    assert(n > 0);
    type_cardinality = n;
}

void Cells::set_type_transition(const std::vector<std::vector<double>> &transit)
{
    assert(transit.size() == type_cardinality);
    for (int i = 0; i != type_cardinality; i++)
    {
        assert(transit[i].size() == type_cardinality);
    }

    type_transition = transit;
}

void Cells::set_maximum_population_size(int n)
{
    assert(n > 0);
    maximum_population_size = n;
}

void Cells::set_initial_population(const std::vector<Cell> &initial_population)
{
    for (auto cell : initial_population)
    {
        assert(cell.type() < type_cardinality);
    }
    current_population = initial_population;
    before_selection_population = initial_population;
}

void Cells::time_evolution(const std::vector<std::vector<double>> &offspring_distribution, const std::vector<std::vector<double>> &next_offspring_distribution)
{
    std::vector<Cell> new_population;
    for (auto cell : current_population)
    {
        std::vector<Cell> daughters = cell.daughters(type_transition, offspring_distribution, mt);

        new_population.insert(new_population.end(), daughters.begin(), daughters.end());
    }

    current_population.clear();
    current_population = new_population;
    before_selection_population = current_population;

    selection();
}

void Cells::selection()
{

    if (current_population.size() > maximum_population_size)
    {
        std::vector<Cell> temp;

        std::sample(current_population.begin(), current_population.end(), std::back_inserter(temp), maximum_population_size, mt);
        current_population = temp;
    }
}

void Cells::record()
{
    for (auto cell : current_population)
    {
        cell.record(out_pop);
    }
    for (auto cell : before_selection_population)
    {
        cell.record(out_pop_full);
    }
}

Cells_Common::Cells_Common(int type_no, int init_common_p_type,
                           const std::vector<std::vector<double>> &type_tran,
                           const std::vector<Cell> &initial_population,
                           int max_pop_size)
    : Cells(type_no, type_tran, initial_population, max_pop_size)
{
    set_common_p_type(init_common_p_type);
}

void Cells_Common::time_evolution(const std::vector<std::vector<double>> &offspring_distribution, const std::vector<std::vector<double>> &next_offspring_distribution)
{

    //new population
    std::vector<Cell> new_population;
    for (auto cell : current_population)
    {

        std::vector<Cell> daughters = cell.daughters(m_common_p_type, type_transition, offspring_distribution, mt);

        new_population.insert(new_population.end(), daughters.begin(), daughters.end());
    }

    current_population.clear();
    current_population = new_population;

    //new shared parent type

    std::vector<double> offspring_no;
    for (auto cell : current_population)
    {
        int no_offsprings = std::discrete_distribution<int>(next_offspring_distribution[cell.type()].begin(), next_offspring_distribution[cell.type()].end())(mt);

        offspring_no.push_back(no_offsprings);
    }
    if (offspring_no.empty())
    {
        m_common_p_type = 0;
    }
    else
    {
        int itr = std::discrete_distribution<int>(offspring_no.begin(), offspring_no.end())(mt);
        m_common_p_type = current_population[itr].type();
    }

    m_common_p_type = current_population.empty() ? 0 : current_population[0].type();

    //update before_selection_population
    before_selection_population = current_population;

    selection();
}

void Cells_Common::set_initial_population(const std::vector<Cell> &initial_population)
{
    Cells::set_initial_population(initial_population);

    m_common_p_type = initial_population.empty() ? 0 : current_population[0].type();
}

Cells_Common::Cells_Common(const Cells &cells)
    : Cells(cells)
{
    m_common_p_type = current_population.empty() ? 0 : current_population[0].type();
}

void Cells_Common::record()
{
    Cells::record();
    *out_shared_p_type << common_p_type() << std::endl;
}

std::vector<Cell_Learn> const Cell_Learn::daughters(const std::vector<std::vector<double>> &offspring_distribution,
                                                    const std::function<void(int, int, int, std::vector<std::vector<double>> &, std::vector<std::vector<double>> &, std::vector<double> &, std::vector<double> &, std::mt19937_64 &)> &learning_rule,
                                                    std::mt19937_64 &mt)
{

    int no_offsprings = std::discrete_distribution<int>(offspring_distribution[m_type].begin(), offspring_distribution[m_type].end())(mt);

    std::vector<Cell_Learn> res;
    for (int i = 0; i != no_offsprings; i++)
    {
        int next_type = std::discrete_distribution<int>(transition[m_type].begin(), transition[m_type].end())(mt);

        //learning
        std::vector<std::vector<double>> temp_jump = ancestral_jump;
        std::vector<std::vector<double>> temp_tran_mat = transition;
        std::vector<double> temp_rep_hist = replication_history;
        std::vector<double> temp_mem = memory;

        learning_rule(type(), next_type, no_offsprings, temp_tran_mat, temp_jump, temp_rep_hist, temp_mem, mt);

        res.push_back(
            Cell_Learn(next_type, m_ID + "S" + std::to_string(i),
                       temp_jump,
                       temp_tran_mat,
                       temp_rep_hist,
                       temp_mem));
    }

    return res;
}

Cell_Learn::Cell_Learn(int type, std::string ID,
                       const std::vector<std::vector<double>> &jump,
                       const std::vector<std::vector<double>> &tran_mat,
                       const std::vector<double> &rep_hist,
                       const std::vector<double> &mem)
{
    m_type = type;
    m_ID = ID;
    set_ancestral_jump(jump);
    set_transition(tran_mat);
    set_memory(mem);
}

void Cell_Learn::record(std::ofstream *out)
{
    *out << id() << " " << type() << std::endl;
    out_mat(ancestral_jump, out);
    out_mat(transition, out);
    //out_vec(replication_history, out);
    out_vec(memory, out);
    *out << std::endl;
}

Cells_Learn::Cells_Learn(int type_no, int max_pop_size,
                         const std::vector<Cell_Learn> &initial_pop,
                         const std::function<void(int, int, int, std::vector<std::vector<double>> &, std::vector<std::vector<double>> &, std::vector<double> &, std::vector<double> &, std::mt19937_64 &)> &rule)
{

    set_type_cardinality(type_no);
    set_maximum_population_size(max_pop_size);
    set_learning_rule(rule);
    set_initial_population(initial_pop);

    //random generator
    std::random_device rnd;
    mt.seed(rnd());
}

void Cells_Learn::selection()
{
    if (current_population.size() > maximum_population_size)
    {
        std::vector<Cell_Learn> temp;

        std::sample(current_population.begin(), current_population.end(), std::back_inserter(temp), maximum_population_size, mt);
        current_population = temp;
    }
}

void Cells_Learn::time_evolution(const std::vector<std::vector<double>> &offspring_distribution, const std::vector<std::vector<double>> &next_offspring_distribution)
{
    std::vector<Cell_Learn> new_population;
    for (auto cell : current_population)
    {
        std::vector<Cell_Learn> daughters = cell.daughters(offspring_distribution, learning_rule, mt);

        new_population.insert(new_population.end(), daughters.begin(), daughters.end());
    }

    before_selection = new_population;
    current_population = new_population;

    selection();
}

void Cells_Learn::record()
{
    for (auto cell : current_population)
    {
        cell.record(out_pop);
    }
    for (auto cell : before_selection)
    {
        cell.record(out_pop_full);
    }
}

Cells_Learn_Common::Cells_Learn_Common(
    int type_no,
    int max_pop_size,
    const std::vector<Cell_Learn> &initial_pop,
    const std::function<void(int, int, int, std::vector<std::vector<double>> &, std::vector<std::vector<double>> &, std::vector<double> &, std::vector<double> &, std::mt19937_64 &)> &rule)
{
    Cells_Learn_Common(Cells_Learn(type_no, max_pop_size, initial_pop, rule));
}

Cells_Learn_Common::Cells_Learn_Common(const Cells_Learn &cell) : Cells_Learn(cell)
{
    //select the spine cell
    //exceptional case
    if (current_population.empty())
    {
        spine_cell = Cell_Learn();
        return;
    }

    set_spine_cell(Cell_Learn(current_population[0]));
    //set lerning historty to the same one (spine)
    for (auto cell : current_population)
    {
        cell.set_transition(spine_cell.transition);
        cell.set_replication_history(spine_cell.replication_history);
        cell.set_ancestral_jump(spine_cell.ancestral_jump);
        cell.set_memory(spine_cell.memory);
    }
}

void Cells_Learn_Common::set_spine_cell(const Cell_Learn &c)
{
    spine_cell = c;
    spine_cell.set_id("0");
}

void Cells_Learn_Common::time_evolution(const std::vector<std::vector<double>> &offspring_distribution, const std::vector<std::vector<double>> &next_offspring_distribution)
{
    std::vector<Cell_Learn> daughters_population;
    std::vector<double> no_grand_dauthers;

    //generate next population
    for (auto cell : current_population)
    {
        std::vector<Cell_Learn> daughters = cell.daughters(offspring_distribution, learning_rule, mt);

        daughters_population.insert(daughters_population.end(), daughters.begin(), daughters.end());
    }

    //count # of daughters of each cell in daughters_population
    for (auto cell : daughters_population)
    {
        std::vector<Cell_Learn> grand_daughters = cell.daughters(offspring_distribution, learning_rule, mt);

        no_grand_dauthers.push_back(grand_daughters.size());
    }

    //select spine cell
    const std::string parent_spine_id = spine_cell.id();
    if (daughters_population.empty())
    {
        spine_cell = Cell_Learn();
    }
    else
    {
        int itr = std::discrete_distribution<int>(no_grand_dauthers.begin(), no_grand_dauthers.end())(mt);
        spine_cell = daughters_population[itr];
    }
    spine_cell.set_id(parent_spine_id + "S0");

    //set transition matrix (and other elements) to common one learned by the spine cell
    for (auto &cell : daughters_population)
    {
        cell.set_transition(spine_cell.transition);
        cell.set_replication_history(spine_cell.replication_history);
        cell.set_ancestral_jump(spine_cell.ancestral_jump);
        cell.set_memory(spine_cell.memory);
    }

    //update population
    before_selection = daughters_population;
    current_population = daughters_population;

    selection();
}

void Cells_Learn_Common::record()
{
    for (auto cell : current_population)
    {
        cell.record(out_pop);
    }
    if (out_pop_full != nullptr)
    {
        for (auto cell : before_selection)
            cell.record(out_pop_full);
    }
    if (out_spine != nullptr)
    {
        spine_cell.record(out_spine);
    }
}

Cells_Infinite::Cells_Infinite(int type_no, const std::vector<std::vector<double>> &transit, const std::vector<double> &initial_pop)
{
    set_type_cardinality(type_no);
    set_type_transition(transit);
    set_initial_population(initial_pop);
    lambdas.push_back(0.0);
}

void Cells_Infinite::set_type_transition(const std::vector<std::vector<double>> &transit)
{
    assert(transit.size() == type_cardinality);
    for (int i = 0; i != type_cardinality; i++)
    {
        assert(transit[i].size() == type_cardinality);
    }

    type_transition = std::vector<std::vector<double>>(type_cardinality, std::vector<double>(type_cardinality));

    //normalizatoin
    for (int i = 0; i != type_cardinality; i++)
    {
        double sum = 0.0;
        for (int j = 0; j != type_cardinality; j++)
        {
            sum += transit[i][j];
        }
        for (int j = 0; j != type_cardinality; j++)
        {
            type_transition[i][j] = transit[i][j] / sum;
        }
    }
}

void Cells_Infinite::set_type_cardinality(int n)
{
    assert(n > 0);
    type_cardinality = n;
}

void Cells_Infinite::set_initial_population(const std::vector<double> &initial_pop)
{
    assert(initial_pop.size() == cardinality());
    current_pop = initial_pop;
}

void Cells_Infinite::record()
{
    if (out_pop_and_lambda == nullptr)
        return;

    *out_pop_and_lambda << lambdas[lambdas.size() - 1] << " ";
    for (auto i : current_pop)
    {
        *out_pop_and_lambda << i << " ";
    }
    *out_pop_and_lambda << std::endl;
}

void Cells_Infinite::time_evolution(const std::vector<std::vector<double>> &offspring_distribution, const std::vector<std::vector<double>> &next_offspring_distribution)
{
    //calculate expectation of offspring_distribution
    std::vector<double> mean_daughters_no(cardinality(), 0.0);
    for (int i = 0; i != cardinality(); i++)
    {
        double normalizer = 0.0;
        double mean = 0.0;
        for (int j = 0; j != offspring_distribution[i].size(); j++)
        {
            normalizer += offspring_distribution[i][j];
            mean += j * offspring_distribution[i][j];
        }
        mean_daughters_no[i] = mean / normalizer;
    }

    std::vector<double> res(cardinality(), 0.0);
    for (int i = 0; i != cardinality(); i++)
    {
        for (int j = 0; j != cardinality(); j++)
        {
            res[j] += current_pop[i] * mean_daughters_no[i] * type_transition[i][j]; //type_transition is assumed to be normalized
        }
    }

    double growth_ratio = 0.0;
    for (auto f : res)
    {
        growth_ratio += f;
    }
    for (int i = 0; i != cardinality(); i++)
    {
        res[i] /= growth_ratio;
    }

    current_pop = res;
    lambdas.push_back(log(growth_ratio));
}

double Cells_Infinite::lambda() const
{
    double res = 0.0;
    int end_time = lambdas.size() - 1;
    for (auto d : lambdas)
    {
        res += d / end_time;
    }
    return res;
}

Cells_Infinite_Common::Cells_Infinite_Common(Cells_Infinite cells)
    : Cells_Infinite(cells)
{
    //random generator
    std::random_device rnd;
    mt.seed(rnd());

    //initial p-type
    std::discrete_distribution<int> disc(cells.current_pop.begin(), cells.current_pop.end());
    common_p_type = disc(mt);
}

Cells_Infinite_Common::Cells_Infinite_Common(int type_no, const std::vector<std::vector<double>> &transit, const std::vector<double> &initial_pop) : Cells_Infinite_Common(Cells_Infinite(type_no, transit, initial_pop))
{
}

void Cells_Infinite_Common::time_evolution(const std::vector<std::vector<double>> &offspring_distribution, const std::vector<std::vector<double>> &next_offspring_distribution)
{
    //calculate expectation of offspring_distribution
    std::vector<double> mean_daughters_no(cardinality(), 0.0);
    for (int i = 0; i != cardinality(); i++)
    {
        double normalizer = 0.0;
        double mean = 0.0;
        for (int j = 0; j != offspring_distribution[i].size(); j++)
        {
            normalizer += offspring_distribution[i][j];
            mean += j * offspring_distribution[i][j];
        }
        mean_daughters_no[i] = mean / normalizer;
    }

    /*
    //next
    std::vector<double> next_mean_daughters_no(cardinality(), 0.0);
    for(int i = 0; i != cardinality(); i++){
        double normalizer = 0.0;
        double mean = 0.0;
        for(int j = 0; j != next_offspring_distribution[i].size(); j++){
            normalizer += next_offspring_distribution[i][j];
            mean += j * next_offspring_distribution[i][j];
        }
        next_mean_daughters_no[i] = mean / normalizer;
    }*/

    //determine common_p_type
    std::vector<double> growth(cardinality(), 0.0);
    double growth_ratio = 0.0;
    for (int i = 0; i != cardinality(); i++)
    {
        double temp = current_pop[i] * mean_daughters_no[i];
        growth_ratio += temp;
        growth[i] = temp;
    }
    std::discrete_distribution<int> dist(growth.begin(), growth.end());
    common_p_type = dist(mt);

    //evolve current_pop
    std::vector<double> res(cardinality(), 0.0);
    for (int i = 0; i != cardinality(); i++)
    {
        for (int j = 0; j != cardinality(); j++)
        {
            res[j] = type_transition[common_p_type][j]; //type_transition is assumed to be normalized
        }
    }

    current_pop = res;
    lambdas.push_back(log(growth_ratio));
    hist_common_p_type.push_back(common_p_type);
}

void Cells_Infinite_Common::record()
{
    if (out_pop_and_lambda == nullptr)
        return;

    *out_pop_and_lambda << lambdas[lambdas.size() - 1] << " " << common_p_type << " ";
    for (auto i : current_pop)
    {
        *out_pop_and_lambda << i << " ";
    }
    *out_pop_and_lambda << std::endl;
}