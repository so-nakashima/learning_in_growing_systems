#include "headers/simulators.h"
#include <map>
#include "headers/analyze.h"
#include "headers/simulator_utility.h"

Lineage<Cell> read_lineage(std::ifstream &in)
{
    std::string id;
    int type;
    std::vector<Cell> current_generation_pop;
    int gen = 1;
    std::vector<std::vector<Cell>> pop;
    std::map<std::string, Cell> pop_map;

    while (in >> id)
    {
        in >> type;
        Cell c(type, id);
        pop_map[id] = c;

        //current generation
        const int no_dots = 0;
        std::count(id.begin(), id.end(), '.');

        if (no_dots + 1 != gen)
        {
            pop.push_back(current_generation_pop);
            current_generation_pop.clear();
            gen++;
        }
        current_generation_pop.push_back(c);
    }

    pop.push_back(current_generation_pop);
    return Lineage<Cell>(pop, pop_map);
}

Lineage<Cell_Learn> read_learning_lineage(int type_no, int memory_no, std::ifstream &in)
{
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

    while (in >> id)
    {
        in >> type;
        readMat(type_no, type_no, ancestral_jump, in);
        readMat(type_no, type_no, transition, in);
        //readVec(type_no, replication_history, in);
        readVec(memory_no, mem, in);

        Cell_Learn c(type, id, ancestral_jump, transition, replication_history, mem);
        pop_map[id] = c;

        const int no_dots = std::count(id.begin(), id.end(), '.');
        if (no_dots + 1 != gen)
        {
            pop.push_back(current_generation_pop);
            current_generation_pop.clear();
            gen++;
        }
        current_generation_pop.push_back(c);
    }

    pop.push_back(current_generation_pop);
    return Lineage<Cell_Learn>(pop, pop_map);
}

namespace Lineage_utility
{
    std::string change_id(const std::string &str, int no_roots)
    {
        std::string res;
        std::string head = str;

        //split str into head and tail if '.' exists (gen >= 1)
        if (std::count(str.begin(), str.end(), '.') != 0)
        {
            //seek the first "."
            const auto fst_dot_index = str.find(".");

            //the head string is constructed later
            res.insert(0, str, fst_dot_index, str.size() - fst_dot_index + 1);

            std::string head;
            head.insert(0, str, 0, fst_dot_index);
        }

        const int new_root_id = std::stoi(head) + no_roots;
        const std::string new_root_id_str = std::to_string(new_root_id);
        res = new_root_id_str + res;

        return res;
    } // namespace Lineage_utility
} // namespace Lineage_utility
