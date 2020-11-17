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
public:
    Environments_base(){};
    ~Environments_base(){};

    virtual int next_state(int t){return 0;};
    virtual int current_state() const {return 0;};
    virtual int cardinality() const {return 1;};
    virtual void record(){};
};

class Environments_Sequence : public Environments_base
{
private:
    std::vector<int> m_sequence;
    int time = 0;
    int m_cardinality = 1;

public:
    Environments_Sequence(int card, const std::vector<int>& seq){m_sequence = seq; m_cardinality = card;};
    ~Environments_Sequence() = default;

    virtual int next_state(int t){time = t + 1; return m_sequence[t];};
    virtual int current_state() const {return m_sequence[time];};
    virtual void record() {};
    virtual int cardinality() const  {return m_cardinality;};
};


class Markov_Environments : public Environments_base
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

    //recording
    std::ofstream* out_env;

public:
    Markov_Environments(const std::vector<std::vector<double>>& transit_mat = std::vector<std::vector<double>>(),
     const int cardinality = 1, const std::vector<double>& initial  = std::vector<double>()); //for time-homogenous Markov chain (construction from transition matrix)
    Markov_Environments(const std::vector<std::vector<std::function<double(int)>>>& tran_func_mat, const int cardinality = 1, const std::vector<double>& initial  = std::vector<double>());  //for time-inhomogenous Markov chain
    ~Markov_Environments();

    virtual int current_state() const {return m_current_state;};
    virtual int cardinality() const {return m_cardinality;}

    //void set_initial(int init);
    void set_initial_distribution(const std::vector<double>& init);
    void set_cardinality(int n);
    void set_transition(const std::vector<std::vector<double>>& tran_mat);
    void set_transition(const std::vector<std::vector< std::function<double(int)> >>& tran_func_mat);
    int next_state(int t);

    //generate sequence of environments
    std::vector<int> generate(const int n);


    void set_env_record(std::ofstream* out){out_env = out;};
    virtual void record();
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

    std::vector<Cell> daughters(const int common_parent_type, const std::vector<std::vector<double>>& type_transition, const std::vector<std::vector<double>>& offspring_distribution,
    std::mt19937_64& mt) const ;

    int type() const {return m_type;}
    std::string id() const {return m_ID;}

    void record(std::ofstream* out_population);
};




class Population
{
public:
    virtual void time_evolution(const std::vector<std::vector<double>>& offspring_distribution, const std::vector<std::vector<double>>& next_offspring_distribution = std::vector<std::vector<double>>()){};
    virtual void record(){};
    //virtual void record_before_selection(){};
    virtual int cardinality() const {return 0;}; //no of type
    virtual int size() const {return 0;}; //no of cells 
};

class Cells : public Population
{
protected:
    int type_cardinality;
    std::vector<std::vector<double>> type_transition; //type_transiton[i][j] = prob. of type transition from i to j
    int maximum_population_size; 
    std::mt19937_64 mt;
    std::vector<Cell> current_population;
    std::vector<Cell> before_selection_population;

    void selection();

    std::ofstream* out_pop;
    std::ofstream* out_pop_full; //when max_pop_size is introduced, record the discarded cell in addition to the selected cells.
    

public:
    Cells(int type_no = 1, 
    const std::vector<std::vector<double>>& type_tran = std::vector<std::vector<double>>(),
    const std::vector<Cell>& initial_population = std::vector<Cell>(),int max_pop_size = std::numeric_limits<int>::max()
    );
    ~Cells();

    void set_type_cardinality(int n);
    void set_type_transition(const std::vector<std::vector<double>>& transit);
    void set_maximum_population_size(int n);
    virtual void set_initial_population(const std::vector<Cell>& initial_population);

    void set_pop_record(std::ofstream* out){out_pop = out;};
    void set_pop_full_record(std::ofstream* out){out_pop_full = out;};


    void time_evolution(const std::vector<std::vector<double>>& offspring_distribution, const std::vector<std::vector<double>>& next_offspring_distribution = std::vector<std::vector<double>>());
    virtual void record();
    //virtual void record_before_selection();

    const int cardinality(){return type_cardinality;};
    int size() const {return current_population.size();}
};

class Cells_Common : public Cells
{
protected:
    int m_common_p_type = 0;
    std::ofstream* out_shared_p_type = NULL;

public:
    Cells_Common(int type_no = 1, int init_common_p_type = 0, 
    const std::vector<std::vector<double>>& type_tran = std::vector<std::vector<double>>(),
    const std::vector<Cell>& initial_population = std::vector<Cell>(),int max_pop_size = std::numeric_limits<int>::max());
    explicit Cells_Common(const Cells& cells);
    ~Cells_Common(){};

    void set_out_shared_p_type(std::ofstream* out){out_shared_p_type = out;};
    virtual void set_initial_population(const std::vector<Cell>& initial_population);

    int common_p_type() const {return m_common_p_type;};
    void set_common_p_type(int p_type){ m_common_p_type = p_type; };
    void time_evolution(const std::vector<std::vector<double>>& offspring_distribution, const std::vector<std::vector<double>>& next_offspring_distribution );
    virtual void record(); 
};


class Cell_Learn : public Cell
{
protected:

public:
    Cell_Learn(int type = 0, std::string ID = "", 
    const std::vector<std::vector<double>>& jump = std::vector<std::vector<double>>(), 
    const std::vector<std::vector<double>>& tran_mat = std::vector<std::vector<double>>(),
    const std::vector<double>& rep_hist = std::vector<double>(),
    const std::vector<double>& mem = std::vector<double>());
    ~Cell_Learn() = default;


    std::vector<std::vector<double>> transition;
    std::vector<std::vector<double>> ancestral_jump;
    std::vector<double> replication_history;
    std::vector<double> memory;

    std::vector<Cell_Learn> const daughters(const std::vector<std::vector<double>>& offspring_distribution,
    const std::function<void(int, int, int, std::vector<std::vector<double>>&, std::vector<std::vector<double>>&, std::vector<double>&,std::vector<double>&, std::mt19937_64&)>& learning_rule,
    std::mt19937_64& mt);

    void set_ancestral_jump(const std::vector<std::vector<double>>& jump){ancestral_jump = jump;};
    void set_transition(const std::vector<std::vector<double>>& tran){transition = tran;};
    void set_replication_history(const std::vector<double>& hist){replication_history = hist;};
    void set_memory(const std::vector<double>& mem){memory = mem;};

    void record(std::ofstream* out_population);
};





class Cells_Learn : public Population
{
private:
    std::function<
    void(int, 
    int ,
    int, 
    std::vector<std::vector<double>>&, std::vector<std::vector<double>>&, 
    std::vector<double>&,
    std::vector<double>&, 
    std::mt19937_64&)> learning_rule;

    int m_type_cardinality = 1;
    std::vector<Cell_Learn> current_population;
    int maximum_population_size; 
    std::mt19937_64 mt;


    std::ofstream* out_pop;
    std::ofstream* out_pop_full; //when max_pop_size is introduced, record the discarded cell in addition to the selected cells.

public:
    Cells_Learn(int type_no = 0, 
    int max_pop_size = 1, 
    const std::vector<Cell_Learn>& initial_pop = {},
    const std::function<void(int, int, int, std::vector<std::vector<double>>&, std::vector<std::vector<double>>&, std::vector<double>&,std::vector<double>&, std::mt19937_64&)>& rule 
    = [](int _, int __, int _p, std::vector<std::vector<double>>& _a, std::vector<std::vector<double>>& _b, std::vector<double>& _c,std::vector<double>& _d, std::mt19937_64& _e){});
    ~Cells_Learn(){};

    void set_learning_rule(const std::function<void(int, int, int,  std::vector<std::vector<double>>&, std::vector<std::vector<double>>&, std::vector<double>&,std::vector<double>&, std::mt19937_64&)>& rule){learning_rule = rule;};

    //called at the moment of the birth of a daughter cell
    //update of transition matrix (arg4), memory (arbitrary scalars, arg7)
    // retrospective mean (arg5), replication_history(type-#daughters pair, arg6) can be used to update, and they are also updated by using my_type(arg1), daughters_type (arg2) and no_of_daughters(arg3)
    //arg 4--arg7 are kept in each cells.

    void set_type_cardinality(const int type_no){m_type_cardinality = type_no;};
    void set_maximum_population_size(const int max_size){maximum_population_size = max_size;};
    void set_initial_population(const std::vector<Cell_Learn>& initial_population){current_population = initial_population;};

    void time_evolution(const std::vector<std::vector<double>>& offspring_distribution, const std::vector<std::vector<double>>& next_offspring_distribution = std::vector<std::vector<double>>());
    void selection();

    void set_pop_record(std::ofstream* out){out_pop = out;};
    void set_pop_full_record(std::ofstream* out){out_pop_full = out;};
    void record(std::ofstream* out_population);

    const int cardinality(){return m_type_cardinality;};
    const int size(){return current_population.size();}
};

class Cells_Infinite : public Population
{
protected:
    int type_cardinality;
    std::vector<std::vector<double>> type_transition;

    std::ofstream* out_pop_and_lambda = nullptr;

public:
    Cells_Infinite() = default;
    ~Cells_Infinite() = default;
    Cells_Infinite(int type_no, const std::vector<std::vector<double>>& transit, const std::vector<double>& initial_pop);


    virtual int size() const {return std::numeric_limits<int>::max();};
    virtual void record();
    virtual int cardinality() const {return type_cardinality;};
    virtual void time_evolution(const std::vector<std::vector<double>>& offspring_distribution, const std::vector<std::vector<double>>& next_offspring_distribution = std::vector<std::vector<double>>());


    std::vector<double> current_pop;

    void set_type_transition(const std::vector<std::vector<double>>& transit);
    void set_type_cardinality(int n);
    void set_initial_population(const std::vector<double>& initial_pop);
    void set_output(std::ofstream* out){out_pop_and_lambda = out;};

    std::vector<double> lambdas; //lambdas[t] = lambda at [0,T]
};

class Cells_Infinite_Common : public Cells_Infinite
{
private:
    int common_p_type = 0;
    std::mt19937_64 mt;

public:
    Cells_Infinite_Common() = default;
    ~Cells_Infinite_Common() = default;
    Cells_Infinite_Common(int type_no, const std::vector<std::vector<double>>& transit, const std::vector<double>& initial_pop);
    explicit Cells_Infinite_Common(Cells_Infinite cells);

    std::vector<int> hist_common_p_type;

    virtual void time_evolution(const std::vector<std::vector<double>>& offspring_distribution, const std::vector<std::vector<double>>& next_offspring_distribution = std::vector<std::vector<double>>());

    virtual void record();

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
    void set_environments(Environments_base* enviornments){env = enviornments;};
    //void set_env_record(std::string path_out);
    void set_end_time(int t){assert(t > 0); end_time = t;};
    void set_population(Population* population);

    //void set_pop_record(std::string paht_out);
    void set_offspring_distributions(const std::vector<std::vector<std::vector<double>>>& offspring_dist); //offspring_dist[y][x][i] = prob of type-x cell having i daughters under env. y.
    const int size_population(){return pop->size();}

    //history of the simulation
    //std::vector<int> enviornments_history; //hoge[i] is the i-the env. state
    //std::vector<std::vector<Cell>> population_history; //hoge[i] is the i-th current population.

private:
    Environments_base* env;
    Population* pop;
    void time_evolution();
    int end_time = 10;

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
