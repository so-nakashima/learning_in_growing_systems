#include "headers/simulators.h"
#include <vector>

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