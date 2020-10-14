#include "matplotlibcpp.h"
#include "../growth_model/headers/simulator_utility.h"
#include <cmath>
#include <functional>
#include <vector>

namespace plt = matplotlibcpp;

class hoge
{
private:
    std::vector<std::vector< std::function<double(int)> >> nyo;
public:
    hoge
(const std::vector<std::vector< std::function<double(int)> >>&);
    ~hoge
();
    void set(const std::vector<std::vector< std::function<double(int)> >>& n);
};

void hoge::set(const std::vector<std::vector< std::function<double(int)> >>& n){
    nyo = n;
}

hoge::hoge(const std::vector<std::vector< std::function<double(int)> >>& n){
    set(n);
}



hoge::~hoge()
{
}


int main() {



    std::vector<std::vector< std::function<double(int)> >> funcs = {{[=](int n){return 0.0;}}};

    hoge fuga(funcs);

    return 0;
}