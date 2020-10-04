#pragma once

#include <iostream>
#include <vector>
#include "simulators.h"
#include <string>
//#include "headers/simulators.h"

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
void readMat(int n, int, std::vector<std::vector<T>>& mat, std::ifstream& in);

template<typename T>
void readVec(std::vector<T>& vec, std::ifstream& in);
/*
format:
n
a_1 a_2 ... a_n
*/

Cells read_cells(std::ifstream& init_cells, std::ifstream& transition);
/*
format
init_cells:
type_no cell_no(= n)
x_1 x_2 ... x_n (types of the cells)

transition:
a_11 a_12 ... a_1n
a_21 a_22 ... a_2n
.
.
a_n1 a_n2 ... a_nn
*/


Environments read_env(std::ifstream& in_env);
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