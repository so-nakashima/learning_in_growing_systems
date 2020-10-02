#include "headers/simulators.h"
#include <cassert>

World_MBPRE::World_MBPRE(){
}

World_MBPRE::~World_MBPRE(){

}

void World_MBPRE::excecute(){
    for(int t = 0; t != end_time; t++){
        time_evolution();

        //record
        record_environments();
    }
}

void World_MBPRE::time_evolution(){
    //environments
    env.next_state();
}

void World_MBPRE::record_environments(){
    *out_env << std::to_string(env.current_state()) << " ";
}

/* void Environments::set_initial(int init){
    initial_state = init;
} */

Environments::~Environments(){
}

void Environments::set_initial_distribution(const std::vector<double>& init){
    assert(init.size() == cardinality);
    for(int i = 0; i != cardinality; i++){
        assert(init[i] >= 0);
    }

    std::discrete_distribution<int> dist(init.begin(), init.end());
    m_current_state = dist(mt);
}

void Environments::set_cardinality(int n){
    cardinality = n;
}

void Environments::set_transition(const std::vector<std::vector<double>>& tran_mat){
    assert(tran_mat.size () == cardinality);
    for(int i = 0; i != cardinality; i++){
        assert(tran_mat[i].size() == cardinality);
        for(int j = 0; j != cardinality; j++){
            assert(tran_mat[i][j] >= 0);
        }
    }
        transition = tran_mat;
}



int Environments::next_state(){
    std::discrete_distribution<int> dist(transition[current_state()].begin(), transition[current_state()].end());

    m_current_state = dist(mt);
    return m_current_state;
}



Environments::Environments(const int cardinality , const std::vector<double>& initial  , const std::vector<std::vector<double>>& transit_mat ){
    //set random generator
    std::random_device rnd;
    mt.seed(rnd());

    set_cardinality(cardinality);

    if(!std::vector<int>().empty()){
        set_initial_distribution(initial);
    }
    else{
        m_current_state = 0;
    }



    //set transition matrix if valid; otherwise initialize by uniform transition
    if(transit_mat.size() != 0)
        transition = transit_mat;
    else
    {
        transition = std::vector<std::vector<double>>(cardinality, std::vector<double>(cardinality, 1.0));
    }
    
}