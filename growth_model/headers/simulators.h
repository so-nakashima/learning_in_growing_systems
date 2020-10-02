#include<vector>
#include <random>
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <limits>

class Environments
{
private:
    //initial data
    int cardinality = 1;
    std::vector<std::vector<double>> transition; //transition[i][j] = prob. of the transition form i to j.
    //std::vector<int> initial_distribution = std::vector<int>(1, 1);

    //states
    int m_current_state = 0;

    //random number geeerators
    std::mt19937_64 mt;  

public:
    Environments(const int cardinality = 1, const std::vector<double>& initial  = std::vector<double>(), const std::vector<std::vector<double>>& transit_mat = std::vector<std::vector<double>>());
    ~Environments();

    const int current_state(){return m_current_state;};

    //void set_initial(int init);
    void set_initial_distribution(const std::vector<double>& init);
    void set_cardinality(int n);
    void set_transition(const std::vector<std::vector<double>>& tran_mat);
    int next_state();
    void record(std::ofstream* out);
};

class Cell
{
private:
    int m_type;
    std::string m_ID; //<hoge>i means this cell is the i-th (0-origin) daughter of the cell with ID <hoge>
public:
    Cell(int type = 0, std::string ID = "");
    ~Cell();

    std::vector<Cell> daughters(int current_env, const std::vector<std::vector<double>>& type_transition, const std::vector<std::vector<double>>& offspring_distribution,
    std::mt19937_64& mt);

    void record(std::ofstream* out_population);
};



class Population
{
public:

    virtual void time_evolution(int current_env){};
    virtual void record(std::ofstream* out_population){};
};

class Cells : public Population
{
private:
    int type_cardinality;
    std::vector<std::vector<double>> type_transition; //type_transiton[i][j] = prob. of type transition from i to j
    std::vector<std::vector<double>> offspring_distribution; //offspring_distribution[i][j] = prob. of type-i-cell having j daughters.
    int maximum_population_size; 
    std::vector<Cell> current_population;
    std::mt19937_64 mt;

public:
    Cells(int type_no = 1, 
    const std::vector<std::vector<double>>& type_tran = std::vector<std::vector<double>>(), 
    const std::vector<std::vector<double>>& offsprings_prob = std::vector<std::vector<double>>(),
    const std::vector<Cell>& initial_population = std::vector<Cell>(),int max_pop_size = std::numeric_limits<int>::max()
    );
    ~Cells();

    void set_type_cardinality(int n);
    void set_type_transition(const std::vector<std::vector<double>>& transit);
    void set_offspring_distribution(const std::vector<std::vector<double>>& offspring_dist);
    void set_maximum_population_size(int n);
    void set_initial_population(const std::vector<Cell>& initial_population);

    void time_evolution(int current_env);
    void record(std::ofstream* out_population);
};


class World_base{
public:
    virtual void execute() {};   
};

class World_MBPRE : public World_base
{
public:
    World_MBPRE(/* args */);
    ~World_MBPRE();

    //run simulation
    void excecute();

    //set intial state
    void set_environments(const Environments& enviornments){env = enviornments;};
    void set_env_record(std::ofstream* out){out_env = out;};
    void set_end_time(int t){assert(t > 0); end_time = t;};
    void set_population(Population* population);
    void set_pop_record(std::ofstream* out){out_pop = out;};

private:
    Environments env;
    Population* pop;
    void time_evolution();
    int end_time = 10;
    std::ofstream* out_env;
    std::ofstream* out_pop;
    void record();
};



