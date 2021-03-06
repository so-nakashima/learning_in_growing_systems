#pragma onece

#include <iomanip>
#include <vector>
#include <map>
#include <string>
#include <functional>
#include <limits>
#include <cmath>
#include "simulators.h"
#include <colormap-shaders/include/colormap/colormap.h>
#include <boost/format.hpp>

//just utility, do not use outside of this class
namespace Lineage_utility
{
    std::string change_id(const std::string &str, int no_roots);
    std::string add_quote(const std::string &str);
} // namespace Lineage_utility

template <typename T>
class Lineage
{
private:
    //the following two data should be private (if possible)
    std::vector<std::vector<T>> m_population; //m_population[i] = cells at time i
    std::map<std::string, T> m_pop_map;       //ID to cell
public:
    Lineage() = default;
    Lineage(const std::vector<std::vector<T>> m_population, std::map<std::string, T> pop_map);
    ~Lineage();
    int endtime() const { return m_population.size(); };
    size_t pop_size(const int i) { return m_population[i].size(); };
    T access(int t, int i) const { return m_population[t][i]; };
    bool const is_id_exists(std::string id) { return m_pop_map.find(id) != m_pop_map.end(); }
    T const lookup(std::string id) { return m_pop_map[id]; }
    void push(const Lineage<T> &lineage);

    //for retrospective process
    T parent(const T &cell);
    std::vector<double> const backward_mean(std::vector<std::function<double(T, T, int)>> funcs, int max_no = std::numeric_limits<int>::max()); //calculate backward mean for all f \in funcs. [0...max_no-1]cell is used at the last time// arguments of func are (parent, current , generation) generation is that of cell (0-origin)
    double const backward_mean(std::function<double(T, T, int)> func, int max_no = std::numeric_limits<int>::max());

    //for graphics
    void const graphic(std::function<double(T)> func, std::ofstream &out, std::ofstream &out_other_data);
    void const graphic(std::function<double(T)> func, const double min, const double max, std::ofstream &out, std::ofstream &out_other_data);
    double const max(std::function<double(T)> func);
    double const min(std::function<double(T)> func);
    //coloar plot for  func(cell) over lineage (graphviz)

    //for growth rate
    double const lambda(int max_pop_no);
};

Lineage<Cell> read_lineage(std::ifstream &in);

Lineage<Cell_Learn> read_learning_lineage(int type_no, std::ifstream &in);

template <typename T>
Lineage<T>::Lineage(const std::vector<std::vector<T>> population, std::map<std::string, T> pop_map)
{
    m_population = population;
    m_pop_map = pop_map;
}

template <typename T>
Lineage<T>::~Lineage()
{
}

template <typename T>
T Lineage<T>::parent(const T &cell)
{
    std::string id = cell.id();
    assert(id.size() > 1);
    // compute parent id by eliminating the tail
    // ex 0S2S12S31 -> 0S2S12S
    while (id[id.size() - 1] != 'S')
    {
        id.pop_back();
    }
    //ex 0S2S12S -> 0S2S12
    id.pop_back();

    return m_pop_map[id];
}

template <typename T>
double const Lineage<T>::backward_mean(std::function<double(T, T, int)> func, int max_no)
{
    std::vector<std::function<double(T, T, int)>> funcs;
    funcs.push_back(func);
    return backward_mean(funcs, max_no)[0];
}

template <typename T>
std::vector<double> const Lineage<T>::backward_mean(std::vector<std::function<double(T, T, int)>> funcs, int max_no)
{
    size_t d = funcs.size();
    std::vector<double> res(d, 0);

    int endt = endtime() - 1;

    size_t search_no = std::min((int)pop_size(endt), max_no);
    for (int i = 0; i != search_no; i++)
    {
        T cell = access(endt, i);
        for (int t = endt; t >= 1; t--)
        {
            T parent_cell = parent(cell);
            for (int j = 0; j != d; j++)
            {
                res[j] += funcs[j](parent_cell, cell, t) / endt / search_no;
            }
            if (t != 0)
            {
                cell = parent_cell;
            }
        }
    }
    return res;
}

template <typename T>
double const Lineage<T>::lambda(int max_pop_no)
{
    int endT = endtime();
    double res = 0.0;
    for (int i = 0; i < endT - 1; i++)
    {
        res += std::log((double)pop_size(i + 1) / std::min(max_pop_no, (int)pop_size(i))) / (endT - 1);
    }
    return res;
}

template <typename T>
void const Lineage<T>::graphic(std::function<double(T)> func, const double min_val, const double max_val, std::ofstream &out, std::ofstream &out_other_data)
{
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
label = "",
style = filled,
color = white,
width = 1.3,
height = 1.3
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
    using namespace colormap;
    MATLAB::Jet jet;

    out_other_data << min_val << " " << max_val << std::endl;

    using namespace Lineage_utility; //to use add_quote

    for (auto pop_t : m_population)
    {
        for (auto c : pop_t)
        {
            //id
            out << add_quote(c.id()) << " "; //to avoid automated split of node by graphvis
            //color
            float normalized_func_val = (func(c) - min_val) / (max_val - min_val);
            Color col = jet.getColor(normalized_func_val);

            int r = std::max(0, std::min(255, (int)(col.r * 256)));
            int g = std::max(0, std::min(255, (int)(col.g * 256)));
            int b = std::max(0, std::min(255, (int)(col.b * 256)));
            std::string color_s = (boost::format("%02x%02x%02x") % r % g % b).str();

            out << "[ fillcolor = \"#" << color_s
                << "\" ]"
                << std::endl;
        }
    }

    //edges
    for (int gen = 1; gen < m_population.size(); gen++)
    {
        for (auto c : m_population[gen])
        {
            T parent_cell = parent(c);
            out << add_quote(parent_cell.id()) << " -> " << add_quote(c.id()) << std::endl;
        }
    }

    //rank

    out << "}" << std::flush;
}

template <typename T>
void const Lineage<T>::graphic(std::function<double(T)> func, std::ofstream &out, std::ofstream &out_other_data)
{
    graphic(func, min(func), max(func), out, out_other_data);
}

template <typename T>
double const Lineage<T>::max(std::function<double(T)> func)
{
    double res = std::numeric_limits<double>::min();
    for (auto pop_t : m_population)
    {
        for (auto t : pop_t)
            res = std::max(res, func(t));
    }
    return res;
}

template <typename T>
double const Lineage<T>::min(std::function<double(T)> func)
{
    double res = std::numeric_limits<double>::max();
    for (auto pop_t : m_population)
    {
        for (auto t : pop_t)
            res = std::min(res, func(t));
    }
    return res;
}

template <typename T>
void Lineage<T>::push(const Lineage<T> &lineage)
{

    const int cureent_lineage_height = m_population.size();
    const int no_roots = (cureent_lineage_height == 0) ? 0 : m_population[0].size();

    for (int gen = 0; gen != lineage.m_population.size(); gen++)
    {
        //extend m_population if lineage is higher than current lineage
        if (gen >= cureent_lineage_height)
        {
            m_population.emplace_back();
        }

        //add cells in lineage to current lineage while modifying id
        for (const auto &cell : lineage.m_population[gen])
        {
            //modify id
            T modifying_cell(cell);
            std::string new_id = Lineage_utility::change_id(cell.id(), no_roots);
            modifying_cell.set_id(new_id);

            //update data
            m_population[gen].push_back(modifying_cell);
            m_pop_map[new_id] = modifying_cell;
        }
    }
}