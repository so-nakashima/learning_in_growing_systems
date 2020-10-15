#pragma once


#include <functional>
#include <vector>
#include <random>
#include <iostream>
#include <fstream>
#include <string>
#include <cassert>
#include <limits>

class Environments_base

{
private:
    virtual int next_state(int t){return 0;};
    virtual void record(std::ofstream* out){};
public:
    Environments_base(){};
    ~Environments_base(){};
};



class Markov_Environments : Environments_base
{
private:
    //initial data
    int m_cardinality = 1;
    //std::vector<std::vector<std::function<double(int)>>> transition; //transition(t)[i][j] = prob. of the transition form i to j at time t.
    std::vector<std::vector< std::function<double(int)> >> transition;
    //std::vector<std::vector<double>> transition;

    //states
    int m_current_state = 0;

    //random number geeerators
    std::mt19937_64 mt;  

public:
    Markov_Environments(const std::vector<std::vector<double>>& transit_mat = std::vector<std::vector<double>>(),
     const int cardinality = 1, const std::vector<double>& initial  = std::vector<double>()); //for time-homogenous Markov chain (construction from transition matrix)
    Markov_Environments(const std::vector<std::vector<std::function<double(int)>>>& tran_func_mat, const int cardinality = 1, const std::vector<double>& initial  = std::vector<double>());  //for time-inhomogenous Markov chain
    ~Markov_Environments();

    const int current_state(){return m_current_state;};
    const int cardinality(){return m_cardinality;}

    //void set_initial(int init);
    void set_initial_distribution(const std::vector<double>& init);
    void set_cardinality(int n);
    void set_transition(const std::vector<std::vector<double>>& tran_mat);
    void set_transition(const std::vector<std::vector< std::function<double(int)> >>& tran_func_mat);
    int next_state(int t);
    void record(std::ofstream* out);
};

class Cell
{
protected:
    int m_type;
    std::string m_ID; //<hoge>i means this cell is the i-th (0-origin) daughter of the cell with ID <hoge>
public:
    Cell(int type = 0, std::string ID = "");
    ~Cell();

    std::vector<Cell> daughters(const std::vector<std::vector<double>>& type_transition, const std::vector<std::vector<double>>& offspring_distribution,
    std::mt19937_64& mt) const ;

    int type() const {return m_type;}
    std::string id() const {return m_ID;}

    void record(std::ofstream* out_population);
};




class Population
{
public:
    virtual void time_evolution(const std::vector<std::vector<double>>& offspring_distribution){};
    virtual void record(std::ofstream* out_population){};
    const virtual int cardinality(){return 0;}; //no of type
    const virtual int size(){return 0;}; //no of cells 
    virtual void selection(){};
};

class Cells : public Population
{
private:
    int type_cardinality;
    std::vector<std::vector<double>> type_transition; //type_transiton[i][j] = prob. of type transition from i to j
    int maximum_population_size; 
    std::mt19937_64 mt;
    std::vector<Cell> current_population;

public:
    Cells(int type_no = 1, 
    const std::vector<std::vector<double>>& type_tran = std::vector<std::vector<double>>(),
    const std::vector<Cell>& initial_population = std::vector<Cell>(),int max_pop_size = std::numeric_limits<int>::max()
    );
    ~Cells();

    void set_type_cardinality(int n);
    void set_type_transition(const std::vector<std::vector<double>>& transit);
    void set_maximum_population_size(int n);
    void set_initial_population(const std::vector<Cell>& initial_population);

    void time_evolution(const std::vector<std::vector<double>>& offspring_distribution);
    void selection();
    void record(std::ofstream* out_population);

    const int cardinality(){return type_cardinality;};
    const int size(){return current_population.size();}
};

class Cell_Learn : Cell
{
protected:
    std::vector<std::vector<double>> ancestral_jump;
    std::vector<std::vector<double>> transition;
    std::vector<double> replication_history;
    std::vector<double> memory;

public:
    Cell_Learn(int type = 0, std::string ID = "", 
    const std::vector<std::vector<double>>& jump = std::vector<std::vector<double>>(), 
    const std::vector<std::vector<double>>& tran_mat = std::vector<std::vector<double>>(),
    const std::vector<double>& rep_hist = std::vector<double>(),
    const std::vector<double>& mem = std::vector<double>());
    ~Cell_Learn(){};

    std::vector<Cell_Learn> const daughters(const std::vector<std::vector<double>>& offspring_distribution,
    const std::function<void(int, int, std::vector<std::vector<double>>&, std::vector<std::vector<double>>&, std::vector<double>&,std::vector<double>&, std::mt19937_64&)>& learning_rule,
    std::mt19937_64& mt);

    void set_ancestral_jump(const std::vector<std::vector<double>>& jump){ancestral_jump = jump;};
    void set_transition(const std::vector<std::vector<double>>& tran){transition = tran;};
    void set_replication_history(const std::vector<double>& hist){replication_history = hist;};
    void set_memory(const std::vector<double>& mem){memory = mem;};

    void record(std::ofstream* out_population);
};

class Cells_Learn : Population
{
private:
    std::function<
    void(int my_type, 
    int no_of_daughters, 
    std::vector<std::vector<double>>&, std::vector<std::vector<double>>&, 
    std::vector<double>&,
    std::vector<double>&, 
    std::mt19937_64&)> learning_rule;

    int m_type_cardinality = 1;
    std::vector<Cell_Learn> current_population;
    int maximum_population_size; 
    std::mt19937_64 mt;

public:
    Cells_Learn(int type_no, int max_pop_size, const std::function<void(int, int, std::vector<std::vector<double>>&, std::vector<std::vector<double>>&, std::vector<double>&,std::vector<double>&, std::mt19937_64&)>& rule 
    = [](int _, int __, std::vector<std::vector<double>>& _a, std::vector<std::vector<double>>& _b, std::vector<double>& _c,std::vector<double>& _d, std::mt19937_64& _e){});
    ~Cells_Learn(){};

    void set_learning_rule(const std::function<void(int, int , std::vector<std::vector<double>>&, std::vector<std::vector<double>>&, std::vector<double>&,std::vector<double>&, std::mt19937_64&)>& rule){learning_rule = rule;};

    //called before replication
    //update of transition matrix (arg3), memory (arbitrary scalars, arg6)
    // retrospective mean (arg4), replication_history(type-#daughters pair, arg5) can be used to update, and they are also updated by using my_type(arg1) and no_of_daughters(arg2)
    //arg 3--arg6 are kept in each cells.

    void set_type_cardinality(const int type_no){m_type_cardinality = type_no;};
    void set_max_pop_size(const int max_size){maximum_population_size = max_size;};
    void set_initial_population(const std::vector<Cell>& initial_population);

    void time_evolution(const std::vector<std::vector<double>>& offspring_distribution);
    void selection();
    void record(std::ofstream* out_population);

    const int cardinality(){return m_type_cardinality;};
    const int size(){return current_population.size();}
};


class World_base{
public:
    virtual void execute(){};   
};

class MBPRE : public World_base
{
public:
    MBPRE(/* args */);
    ~MBPRE();

    //run simulation
    void excecute();

    //set intial state
    void set_environments(Markov_Environments* enviornments){env = enviornments;};
    void set_env_record(std::ofstream* out){out_env = out;};
    //void set_env_record(std::string path_out);
    void set_end_time(int t){assert(t > 0); end_time = t;};
    void set_population(Population* population);
    void set_pop_record(std::ofstream* out){out_pop = out;};
    void set_pop_full_record(std::ofstream* out){out_pop_full = out;};
    //void set_pop_record(std::string paht_out);
    void set_offspring_distributions(const std::vector<std::vector<std::vector<double>>>& offspring_dist); //offspring_dist[y][x][i] = prob of type-x cell having i daughters under env. y.
    const int size_population(){return pop->size();}

    //history of the simulation
    //std::vector<int> enviornments_history; //hoge[i] is the i-the env. state
    //std::vector<std::vector<Cell>> population_history; //hoge[i] is the i-th current population.

private:
    Markov_Environments* env;
    Population* pop;
    void time_evolution();
    int end_time = 10;
    std::ofstream* out_env;
    std::ofstream* out_pop;
    std::ofstream* out_pop_full; //when max_pop_size is introduced, record the discarded cell in addition to the selected cells.
    void record();
    std::vector<std::vector<std::vector<double>>> offspring_distributions;
    int time = 0;
};



//utilities
inline void out_mat(const std::vector<std::vector<double>>& mat, std::ofstream* out){
    for(int i = 0; i != mat.size(); i++){
        for(int j = 0; j != mat[i].size(); j++){
            *out << mat[i][j] << " ";
        }
        *out << std::endl;
    }
}

inline void out_vec(const std::vector<double>& vec, std::ofstream* out){
    for(int i = 0; i != vec.size(); i++){
        *out << vec[i] << " ";
    }
    *out << std::endl;
}
