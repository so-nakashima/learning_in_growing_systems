#include "headers/simulators.h"
#include <map>
#include "headers/analyze.h"
#include "headers/simulator_utility.h"

Lineage<Cell> read_lineage(std::ifstream& in){
    std::string id;
    int type;
    std::vector<Cell> current_generation_pop;
    int gen = 1;
    std::vector<std::vector<Cell>> pop;
    std::map<std::string, Cell> pop_map;

    while(in >> id){
        in >> type;
        Cell c(type, id);
        pop_map[id] = c;
        
        if(id.size() != gen){
            pop.push_back(current_generation_pop);
            current_generation_pop.clear();
            gen++;
        }
        current_generation_pop.push_back(c);
    }

    pop.push_back(current_generation_pop);
    return Lineage<Cell>(pop, pop_map);
}

Lineage<Cell_Learn> read_learning_lineage(int type_no, int memory_no, std::ifstream& in){
    std::string id;
    int type;
    std::vector<Cell_Learn> current_generation_pop;
    int gen = 1;
    std::vector<std::vector<Cell_Learn>> pop;
    std::map<std::string, Cell_Learn> pop_map;


    std::vector<std::vector<double>> ancestral_jump;
    std::vector<std::vector<double>> transition;
    std::vector<double> replication_history;
    std::vector<double> mem;

    while(in >> id){
        in >> type;
        readMat(type_no, type_no, ancestral_jump, in);
        readMat(type_no, type_no, transition, in);
        //readVec(type_no, replication_history, in);
        readVec(memory_no, mem, in);

        Cell_Learn c(type, id, ancestral_jump, transition, replication_history, mem);
        pop_map[id] = c;
        
        if(id.size() != gen){
            pop.push_back(current_generation_pop);
            current_generation_pop.clear();
            gen++;
        }
        current_generation_pop.push_back(c);
    }

    pop.push_back(current_generation_pop);
    return Lineage<Cell_Learn>(pop, pop_map);
}
