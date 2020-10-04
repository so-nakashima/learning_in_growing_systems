#include "headers/simulator_utility.h"

template<typename T>
void readMat(std::vector<std::vector<T>>& mat, std::ifstream& in){ 
    int n, m;
    in >> n >> m;

    readMat(n,m,mat, in);
}

template<typename T>
void readMat(int n, int m, std::vector<std::vector<T>>& mat, std::ifstream& in){ 
    mat.clear();
    for(int i = 0; i != n; i++){
        std::vector<T> temp;
        for(int j = 0; j != m; j++){
            T elem;
            in >> elem;
            temp.push_back(elem);
        }
        mat.push_back(temp);
    }
}



template<typename T>
void readVec(std::vector<T>& vec, std::ifstream& in){ 

    int n;
    in >> n;
    readVec(n, vec, in);
}

template<typename T>
void readVec(int n, std::vector<T>& vec, std::ifstream& in){ 
    vec.clear();
    for(int i = 0; i != n; i++){
       T temp;
       in >> temp;
       vec.push_back(temp);
    }   
}


Cells read_cells(std::ifstream& init_cells, std::ifstream& transition){
    //std::ifstream init_cells(path_init_cells);
    //std::ifstream transition(path_transition);

    std::vector<Cell> init_pop;
    int type_no, cell_no;
    init_cells >> type_no >> cell_no;
    for(int i = 0; i != cell_no; i++){
        int temp;
        init_cells >> temp;
        init_pop.emplace_back(temp, std::to_string(i));
    }

    std::vector<std::vector<double>> transit;
    readMat<double>(type_no, type_no, transit, transition);


    return Cells(type_no, transit, init_pop);
}

Environments read_env(std::ifstream& in_env){
    //std::ifstream in_env(path_in_env);

    int env_state_no;
    in_env >> env_state_no;

    std::vector<double> initial_env;
    readVec(env_state_no, initial_env, in_env);

    std::vector<std::vector<double>> tran;
    readMat(env_state_no, env_state_no, tran, in_env);



    Environments env;
    env.set_cardinality(env_state_no);
    env.set_initial_distribution(initial_env);
    env.set_transition(tran);


    return env;
}


