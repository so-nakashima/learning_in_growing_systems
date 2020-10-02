#include<vector>
#include <random>
#include <iostream>
#include <fstream>

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
    void set_end_time(int t){end_time = t;};


private:
    Environments env;
    void time_evolution();
    void record_environments();
    int end_time = 10;
    std::ofstream* out_env;
};



