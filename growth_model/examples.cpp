#include "headers/simulators.h"
#include <vector>
#include "headers/simulator_utility.h"

void test_env(){
    MBPRE w;

    w.set_end_time(30);
    
    std::vector<double> init = {1.0, 1.0};
    std::vector<std::vector<double>> tran = {{1.0, 0}, {0.5, 0.5}};
    Environments env(2, init, tran);

    std::ofstream out_env(".//res//environments.dat");
    w.set_env_record(&out_env);
    w.set_environments(env);

    w.excecute();
}

void test_env_cells(){
        MBPRE w;

    w.set_end_time(20);
    
    //env
    std::vector<double> init = {1.0, 1.0};
    std::vector<std::vector<double>> env_tran = {{0.5, 0.5}, {0.5, 0.5}};
    Environments env(2, init, env_tran);

    std::ofstream out_env(".//res//environments.dat");
    w.set_env_record(&out_env);
    w.set_environments(env);


    //cells
    std::vector<std::vector<double>> type_tran = {{0.8, 0.2}, {0.2,0.8}};
    std::vector<Cell> init_pop = {Cell(0, "0")};//, Cell(0, "1")};

    Cells cells(2, type_tran, init_pop);
    std::ofstream out_pop(".//res//population.dat");
    cells.set_maximum_population_size(5);

    w.set_population(&cells);
    w.set_pop_record(&out_pop);


    //replication
    std::vector<std::vector<std::vector<double>>> replication = {{{0.0, 0.0, 1.0}, {0.0, 0.0, 1.0}}, {{0.0, 0.0, 1.0},{0.0, 0.0, 1.0}}};
    w.set_offspring_distributions(replication);

    w.excecute();
}

void file_read_test(){
    MBPRE w;
    std::ifstream in_env(".//experiments//sim_1//env.dat");
    Environments env = read_env(in_env);
    std::ifstream in_cells(".//experiments//sim_1//initial_cells.dat");
    std::ifstream in_cell_tran(".//experiments//sim_1//cell_type_tran.dat");
    Cells cells = read_cells(in_cells, in_cell_tran);

    w.set_end_time(10);
    std::ofstream out_env(".//experiments//sim_1//res//env.dat");
    w.set_env_record(&out_env);
    w.set_environments(env);
    w.set_population(&cells);
    std::ofstream out_pop(".//experiments//sim_1//res//pop.dat");
    w.set_pop_record(&out_pop);

    //replication
    std::vector<std::vector<std::vector<double>>> replication = {{{0.0, 0.0, 1.0}, {0.0, 0.0, 1.0}}, {{0.0, 0.0, 1.0},{0.0, 0.0, 1.0}}};
    w.set_offspring_distributions(replication);

    w.excecute();
}