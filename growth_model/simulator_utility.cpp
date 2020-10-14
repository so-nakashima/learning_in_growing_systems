#include "headers/simulator_utility.h"

Cells read_cells(std::ifstream& init_cells, std::ifstream& transition){
    //std::ifstream init_cells(path_init_cells);
    //std::ifstream transition(path_transition);

    std::vector<Cell> init_pop;
    int type_no, cell_no, max_cell_no;
    init_cells >> type_no >> cell_no >> max_cell_no;
    for(int i = 0; i != cell_no; i++){
        int temp;
        init_cells >> temp;
        init_pop.emplace_back(temp, std::to_string(i));
    }

    std::vector<std::vector<double>> transit;
    readMat<double>(type_no, type_no, transit, transition);

    Cells res(type_no, transit, init_pop); 
    if(max_cell_no > 0){
        res.set_maximum_population_size(max_cell_no);
    }
    
    return res;
}

Markov_Environments read_env(std::ifstream& in_env){
    //std::ifstream in_env(path_in_env);

    int env_state_no;
    in_env >> env_state_no;

    std::vector<double> initial_env;
    readVec(env_state_no, initial_env, in_env);

    std::vector<std::vector<double>> tran;
    readMat(env_state_no, env_state_no, tran, in_env);



    Markov_Environments env;
    env.set_cardinality(env_state_no);
    env.set_initial_distribution(initial_env);
    env.set_transition(tran);


    return env;
}


MBPRE read_mbpre(std::ifstream& in){
    int endtime;
    std::string dummy;
    in >> dummy >> endtime;
    MBPRE w;
    w.set_end_time(endtime);
    return w;
}
