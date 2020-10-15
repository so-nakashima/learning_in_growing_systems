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

std::map<std::string, std::string> read_parameters(std::ifstream& in){
    std::string temp, temp2;
    std::map<std::string, std::string> res;
    while(in >> temp >> temp2){
        res[temp] = temp2;
    }
    return res;
}

Cells_Learn read_cells_learn(std::ifstream& in){
    //std::ifstream init_cells(path_init_cells);
    //std::ifstream transition(path_transition);

    std::vector<Cell_Learn> init_pop;
    int type_no, cell_no, max_cell_no;
    in >> type_no >> cell_no >> max_cell_no;
    assert(type_no > 0 && cell_no >= 0);

    for(int i = 0; i != cell_no; i++){
        int type; in >> type;
        std::vector<std::vector<double>> ancestral_jump;
        std::vector<std::vector<double>> transit;
        std::vector<double> replication_history;
        std::vector<double> mem;
        readMat<double>(type_no, type_no, ancestral_jump, in);
        readMat<double>(type_no, type_no, transit, in);
        int vec_size; in >> vec_size;
        readVec<double>(vec_size, replication_history, in);
        in >> vec_size;
        readVec<double>(vec_size, mem, in);

        init_pop.emplace_back(type, std::to_string(i), ancestral_jump, transit, replication_history, mem);
    }

    Cells_Learn res(type_no, max_cell_no, init_pop); 
    if(max_cell_no < 0){
        res.set_maximum_population_size(std::numeric_limits<int>::max());
    }
    
    return res;
}