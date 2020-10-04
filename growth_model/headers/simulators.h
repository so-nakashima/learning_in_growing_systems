#pragma once

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
    int m_cardinality = 1;
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
    const int cardinality(){return m_cardinality;}

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

    std::vector<Cell> daughters(const std::vector<std::vector<double>>& type_transition, const std::vector<std::vector<double>>& offspring_distribution,
    std::mt19937_64& mt);

    const int type(){return m_type;}

    void record(std::ofstream* out_population);
};



class Population
{
public:
    virtual void time_evolution(const std::vector<std::vector<double>>& offspring_distribution){};
    virtual void record(std::ofstream* out_population){};
    const virtual  int cardinality(){return 0;};
};

class Cells : public Population
{
private:
    int type_cardinality;
    std::vector<std::vector<double>> type_transition; //type_transiton[i][j] = prob. of type transition from i to j
    int maximum_population_size; 
    std::vector<Cell> current_population;
    std::mt19937_64 mt;

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
    void record(std::ofstream* out_population);

    const int cardinality(){return type_cardinality;};
};


class World_base{
public:
    virtual void execute() {};   
};

class MBPRE : public World_base
{
public:
    MBPRE(/* args */);
    ~MBPRE();

    //run simulation
    void excecute();

    //set intial state
    void set_environments(const Environments& enviornments){env = enviornments;};
    void set_env_record(std::ofstream* out){out_env = out;};
    //void set_env_record(std::string path_out);
    void set_end_time(int t){assert(t > 0); end_time = t;};
    void set_population(Population* population);
    void set_pop_record(std::ofstream* out){out_pop = out;};
    //void set_pop_record(std::string paht_out);
    void set_offspring_distributions(const std::vector<std::vector<std::vector<double>>>& offspring_dist); //offspring_dist[y][x][i] = prob of type-x cell having i daughters under env. y.

private:
    Environments env;
    Population* pop;
    void time_evolution();
    int end_time = 10;
    std::ofstream* out_env;
    std::ofstream* out_pop;
    void record();
    std::vector<std::vector<std::vector<double>>> offspring_distributions;
};



