#include "headers/simulators.h"
#include <map>

void readPopulation_cell(std::vector<std::vector<Cell>>& pop, std::map<std::string, Cell>& pop_map, std::ifstream& in){
    pop.clear();
    pop_map.clear();

    std::string id;
    int type;
    std::vector<Cell> temp;
    int gen = 1;

    while(in >> id){
        in >> type;
        Cell c(type, id);
        pop_map[id] = c;

    }
}