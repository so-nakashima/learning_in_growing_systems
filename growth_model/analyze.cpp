#include "headers/simulators.h"
#include <map>
#include "headers/analyze.h"

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