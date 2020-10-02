#include "headers/simulators.h"
#include <cassert>

void Environments::set_initial(int init){
    initial_state = init;
}

void Environments::set_cardinality(int n){
    cardinality = n;
    uniform = std::uniform_int_distribution(0, cardinality - 1);
}

void Environments::set_transition(const std::vector<std::vector<int>>& tran_mat){
    assert(transition.size () == cardinality);
    for(int i = 0; i != cardinality; i++){
        assert(transition[i].size() == cardinality);
    }

    transition = tran_mat;
}



int Environments::next_state(){
    std::discrete_distribution dist(transition[current_state].begin(), transition[current_state].end());

    current_state = dist(mt);
    return current_state;
}


Environments::Environments(/* args */)
{
    std::random_device rnd;
    mt.seed(rnd());

    uniform = std::uniform_int_distribution(0, cardinality);
}
