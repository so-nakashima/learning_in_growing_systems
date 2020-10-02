#include "headers/simulators.h"
#include <vector>

void test_env(){
    World_MBPRE w;

    w.set_end_time(20);
    
    std::vector<double> init = {1.0, 1.0};
    std::vector<std::vector<double>> tran = {{0.5, 0.5}, {0.5, 0.5}};
    Environments env(2, init, tran);

    std::ofstream out_env(".//res//environments.dat");
    w.set_env_record(&out_env);
    w.set_environments(env);

    w.excecute();
}