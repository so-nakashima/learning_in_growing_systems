#include "matplotlibcpp.h"
#include "../growth_model/headers/simulator_utility.h"
#include <cmath>

namespace plt = matplotlibcpp;

int main() {
    std::vector<double> x;
    for(int t = -30; t != 31; t++){
        //x.push_back(log((0.5 + t * 0.01) / 0.5));
        x.push_back(t * 0.01);
    }
    
    std::vector<double> y;
    std::ifstream in("..//growth_model//experiments//sim_1//res//analyze//lambda_curve.dat");
    readVec(61, y, in);

    plt::plot(x,y);

    std::vector<double> z;
    for(int t = -30; t != 31; t++){
        /*
        z.push_back(0.2922 * log((0.5 + t * 0.01) / 0.5)
                    + 0.20733 * log((0.5 - t * 0.01) / 0.5)
                    + 0.20278 * log((0.5 - t * 0.01) / 0.5)
                    + 0.29319 * log((0.5 + t * 0.01) / 0.5)
        + y[10]);*/
        z.push_back(0.2922 *  t * 0.01 / 0.5
                    - 0.20733 * t * 0.01 / 0.5
                    - 0.20278 * t * 0.01 / 0.5
                    + 0.29319 * t * 0.01 / 0.5
        + y[30]);
    }

    plt::plot(x, z);

    plt::save("./basic.png");

    plt::show();
    plt::clf();

    std::vector<double> X(21),Y(21),Z(21);
    copy(x.begin() + 20, x.begin() + 41, X.begin());
    copy(y.begin() + 20, y.begin() + 41, Y.begin());
    copy(z.begin() + 20, z.begin() + 41, Z.begin());

    plt::plot(X,Y);
    plt::plot(X,Z);
    plt::save("./basic_expand.png");
    plt::show();
}