#pragma onece

#include <vector>
#include <map>
#include <string>
#include <functional>
#include <limits>
#include <cmath>

template<typename T> class Lineage
{
private:
    //the following two data should be private (if possible)
    std::vector<std::vector<T>> m_population; //m_population[i] = cells at time i
    std::map<std::string, T> m_pop_map; //ID to cell
public:
    Lineage(const std::vector<std::vector<T>> m_population, std::map<std::string, T> pop_map);
    ~Lineage();
    int const endtime(){return m_population.size()};
    size_t const pop_size(const int){return m_population[i].size()};
    T const access(int t, int i){return m_population[t][i];};
    bool const is_id_exists(std::string id){return m_pop_map.find(id) != m_pop_map.end();}
    T const lookup(std::string id){return m_pop_map[id];}

    //for retrospective process
    T const parent(const T& cell);
    std::vector<double> const backward_mean(std::vector<std::function<double(T)>> funcs, int max_no = std::numeric_limits<int>::max()); //calculate backward mean for all f \in funcs. [0...max_no-1]cell is used at the last time
    double const backward_mean(std::function<double(T)> func, int max_no = std::numeric_limits<int>::max());
    void graphic(){}; //to be implemented via graphviz

    //for growth rate
    double const lambda();
};

template<typename T> Lineage<T>::Lineage(const std::vector<std::vector<T>> population, std::map<std::string, T> pop_map){
    m_population = population;
    m_pop_map = pop_map;
}

template<typename T> Lineage<T>::~Lineage()
{
}

template<typename T> T const Lineage<T>::parent(const T& cell)
{
    assert(cell.id().size() > 1);
    return m_pop_map[cell.id().pop_back()];
}

template<typename T> double const Lineage<T>::backward_mean(std::function<double(T)> func, int max_no){
    std::vector<std::function<double(T)>> funcs;
    funcs.push_back(func);
    return backward_mean(funcs, max_no)[0];
}

template<typename T> std::vector<double> const Lineage<T>::backward_mean(std::vector<std::function<double(T)>> funcs, int max_no){
    size_t d = funct.size();
    std::vector<double> res(d, 0);

    int endt = endtime() - 1;
    size_t search_no = std::min(pop_size(endt) , max_no);
    for(int i = 0; i != search_no; i++){
        T cell = access(endt, i);
        for(int t = endt; t >= 0; t--){
            for(int j = 0; j != d; j++){
                res[j] += funcs[j](cell) / endt / search_no;
            }
            if(t != 0){
                cell = parent(cell);
            }
        }
    }
    return res;
}

template<typename T> double const Lineage<T>::lambda(){
    int endT = endtime();
    double res = 0.0;
    for(int i = 0; i < endT - 1; i++){
         res += std::log( (double) pop_size(i+1) / pop_size(i)  ) / (endT - 1);
    }
    return res;
}
