#pragma onece

#include <vector>
#include <map>
#include <string>
#include <functional>
#include <limits>
#include <cmath>
#include "simulators.h"

template<typename T> class Lineage
{
private:
    //the following two data should be private (if possible)
    std::vector<std::vector<T>> m_population; //m_population[i] = cells at time i
    std::map<std::string, T> m_pop_map; //ID to cell
public:
    Lineage(const std::vector<std::vector<T>> m_population, std::map<std::string, T> pop_map);
    ~Lineage();
    int const endtime(){return m_population.size();};
    size_t const pop_size(const int i){return m_population[i].size();};
    T const access(int t, int i){return m_population[t][i];};
    bool const is_id_exists(std::string id){return m_pop_map.find(id) != m_pop_map.end();}
    T const lookup(std::string id){return m_pop_map[id];}

    //for retrospective process
    T const parent(const T& cell);
    std::vector<double> const backward_mean(std::vector<std::function<double(T, T, int)>> funcs, int max_no = std::numeric_limits<int>::max()); //calculate backward mean for all f \in funcs. [0...max_no-1]cell is used at the last time// arguments of func are (parent, current , generation) generation is that of cell (0-origin)
    double const backward_mean(std::function<double(T, T, int)> func, int max_no = std::numeric_limits<int>::max());
    void const graphic(std::function<double(T)> func, std::ofstream& out);
    //coloar plot for  func(cell) over lineage (graphviz)

    //for growth rate
    double const lambda(int max_pop_no);
};

Lineage<Cell> read_lineage(std::ifstream& in);

Lineage<Cell_Learn> read_learning_lineage(int type_no, int memory_no, std::ifstream& in);


template<typename T> Lineage<T>::Lineage(const std::vector<std::vector<T>> population, std::map<std::string, T> pop_map){
    m_population = population;
    m_pop_map = pop_map;
}

template<typename T> Lineage<T>::~Lineage()
{
}

template<typename T> T const Lineage<T>::parent(const T& cell)
{
    std::string id = cell.id();
    assert(id.size() > 1);
    id.pop_back();
    return m_pop_map[id];
}

template<typename T> double const Lineage<T>::backward_mean(std::function<double(T, T, int)> func, int max_no){
    std::vector<std::function<double(T, T, int)>> funcs;
    funcs.push_back(func);
    return backward_mean(funcs, max_no)[0];
}

template<typename T> std::vector<double> const Lineage<T>::backward_mean(std::vector<std::function<double(T, T,int)>> funcs, int max_no){
    size_t d = funcs.size();
    std::vector<double> res(d, 0);

    int endt = endtime() - 1;

    size_t search_no = std::min((int)pop_size(endt) , max_no);
    for(int i = 0; i != search_no; i++){
        T cell = access(endt, i);
        for(int t = endt; t >= 1; t--){
            T parent_cell = parent(cell);
            for(int j = 0; j != d; j++){
                res[j] += funcs[j](parent_cell, cell, t) / endt / search_no;
            }
            if(t != 0){
                cell = parent_cell;
            }
        }
    }
    return res;
}

template<typename T> double const Lineage<T>::lambda(int max_pop_no){
    int endT = endtime();
    double res = 0.0;
    for(int i = 0; i < endT - 1; i++){
         res += std::log( (double) pop_size(i+1) / std::min(max_pop_no, (int)pop_size(i))  ) / (endT - 1);
    }
    return res;
}


template<typename T> void const Lineage<T>::graphic(std::function<double(T)> func, std::ofstream& out){
    out << "digraph lineage { \n";

    //graph
    std::string graph_property = R"(graph [
layout = dot
    ]
    )";
    out << graph_property;

    //node
    std::string node_property = R"(node [
shape = circle,
label = "";
    ]
    )";
    out << node_property;

    //edge
    std::string edge_property = R"(edge [
dir = none
    ]
    )";
    out << edge_property;



    //generate node and edages from lienage data
    //nodes
    for(auto pop_t : m_population){
        for(auto c : pop_t){
            out << c.id() << std::endl;
        }
    }
    //edges
    for(auto pop_t : m_population){
        for(auto c : pop_t){
            if(c.id().size() > 1){
                T parent_cell = parent(c);
                out << parent_cell.id() << " -> " << c.id() << std::endl;
            }
        }
    }

    //rank



    out << "}";
}