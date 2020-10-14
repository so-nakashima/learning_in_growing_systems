#pragma once

#include <iostream>
#include <vector>
#include "simulators.h"

template<typename T>
void readMat(std::vector<std::vector<T>>& mat, std::ifstream& in);
/*
format:
n  m
a_11 a_12 ... a_1m
a_21 a_22 ... a_2m
.
.
a_n1 a_n2 ... a_nm
*/
template<typename T>
void readMat(int n, int m, std::vector<std::vector<T>>& mat, std::ifstream& in);

template<typename T>
void readVec(std::vector<T>& vec, std::ifstream& in);
/*
format:
n
a_1 a_2 ... a_n
*/

template<typename T>
void read3DTensor(int n, int m,  int l, std::vector< std::vector< std::vector<T> > >& tensor, std::ifstream& in);

template<typename T>
void read3DTensor(std::vector< std::vector< std::vector<T> > >& tensor, std::ifstream& in);

Cells read_cells(std::ifstream& init_cells, std::ifstream& transition);
/*
format
init_cells:
type_no cell_no(= n) max_cell_no(-1 is INF)
x_1 x_2 ... x_n (types of the cells)

transition:
a_11 a_12 ... a_1n
a_21 a_22 ... a_2n
.
.
a_n1 a_n2 ... a_nn
*/


Markov_Environments read_env(std::ifstream& in_env);
/*
format:
env_state_no(=n)
init_dist_0 init_dist_1 ... init_dist_n
a_11(transition) a_12 ... a_1n
a_21 a_22 ... a_2n
.
.
a_n1 a_n2 ... a_nn
*/


MBPRE read_mbpre(std::ifstream& in);
/*
format:
endtime <value-of-endtime>


the first collum is dummy for readability
*/



template<typename T>
void readMat(std::vector<std::vector<T>>& mat, std::ifstream& in){ 
    int n, m;
    in >> n >> m;

    readMat(n,m,mat, in);
}

template<typename T>
void readMat(int n, int m, std::vector<std::vector<T>>& mat, std::ifstream& in){ 
    mat.clear();
    for(int i = 0; i != n; i++){
        std::vector<T> temp;
        for(int j = 0; j != m; j++){
            T elem;
            in >> elem;
            temp.push_back(elem);
        }
        mat.push_back(temp);
    }
}





template<typename T>
void readVec(int n, std::vector<T>& vec, std::ifstream& in){ 
    vec.clear();
    for(int i = 0; i != n; i++){
       T temp;
       in >> temp;
       vec.push_back(temp);
    }   
}

template<typename T>
void readVec(std::vector<T>& vec, std::ifstream& in){ 

    int n;
    in >> n;
    readVec(n, vec, in);
}

template<typename T>
void read3DTensor(int n, int m,  int l, std::vector<std::vector<std::vector<T>>>& tensor, std::ifstream& in){ 
    tensor.clear();
    for(int i = 0; i != n; i++){
        std::vector<std::vector<T>> mat;
        for(int j = 0; j != m; j++){
            std::vector<T> vec;
            for(int k = 0; k != l; k++){
                T temp;
                in >> temp;
                vec.push_back(temp);
            }
            mat.push_back(vec);
        }
        tensor.push_back(mat);
    }
}

template<typename T>
void read3DTensor(std::vector<std::vector<std::vector<T>>>& tensor, std::ifstream& in){ 
    int n,m,l;
    in >> n >> m >> l;
    read3DTensor(n,m,l, tensor, in);
}