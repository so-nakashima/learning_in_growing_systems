#pragma once

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include "simulators.h"

template <typename T>
void readMat(std::vector<std::vector<T>> &mat, std::ifstream &in);
/*
format:
n  m
a_11 a_12 ... a_1m
a_21 a_22 ... a_2m
.
.
a_n1 a_n2 ... a_nm
*/
template <typename T>
void readMat(const int n, const int m, std::vector<std::vector<T>> &mat, std::ifstream &in);

template <typename T>
void readVec(std::vector<T> &vec, std::ifstream &in);
/*
format:
n
a_1 a_2 ... a_n
*/

template <typename T>
void read3DTensor(const int n, const int m, const int l, std::vector<std::vector<std::vector<T>>> &tensor, std::ifstream &in);

template <typename T>
void read3DTensor(std::vector<std::vector<std::vector<T>>> &tensor, std::ifstream &in);

Cells read_cells(std::ifstream &init_cells, std::ifstream &transition);
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

Markov_Environments read_env(std::ifstream &in_env);
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

MBPRE read_mbpre(std::ifstream &in);
/*
format:
endtime <value-of-endtime>


the first collum is dummy for readability
*/

template <typename T>
void readMat(std::vector<std::vector<T>> &mat, std::ifstream &in)
{
    int n, m;
    in >> n >> m;

    readMat(n, m, mat, in);
}

template <typename T>
void readMat(const int n, const int m, std::vector<std::vector<T>> &mat, std::ifstream &in)
{
    mat.clear();
    for (int i = 0; i != n; i++)
    {
        std::vector<T> temp;
        for (int j = 0; j != m; j++)
        {
            T elem;
            in >> elem;
            temp.push_back(elem);
        }
        mat.push_back(temp);
    }
}

template <typename T>
void readVec(const int n, std::vector<T> &vec, std::ifstream &in)
{
    vec.clear();
    for (int i = 0; i != n; i++)
    {
        T temp;
        in >> temp;
        vec.push_back(temp);
    }
}

template <typename T>
void readVec(std::vector<T> &vec, std::ifstream &in)
{

    int n;
    in >> n;
    readVec(n, vec, in);
}

template <typename T>
void read3DTensor(const int n, const int m, const int l, std::vector<std::vector<std::vector<T>>> &tensor, std::ifstream &in)
{
    tensor.clear();
    for (int i = 0; i != n; i++)
    {
        std::vector<std::vector<T>> mat;
        for (int j = 0; j != m; j++)
        {
            std::vector<T> vec;
            for (int k = 0; k != l; k++)
            {
                T temp;
                in >> temp;
                vec.push_back(temp);
            }
            mat.push_back(vec);
        }
        tensor.push_back(mat);
    }
}

template <typename T>
void read3DTensor(std::vector<std::vector<std::vector<T>>> &tensor, std::ifstream &in)
{
    int n, m, l;
    in >> n >> m >> l;
    read3DTensor(n, m, l, tensor, in);
}

std::map<std::string, std::string> read_parameters(std::ifstream &);
//to be converted from string to desired type

Cells_Learn read_cells_learn(std::ifstream &init_cells);
Cells_Learn *new_cells_learn_from_read(std::ifstream &in);
Cells_Learn_Common *new_cells_learn_common_from_read(std::ifstream &in);
/*
format: 
<type_no> <cell_no> <max_cell_no>
//max_cell_no < 0 means no limits

//repeat the following format for one cell for <cell_no> times:
<type>
<ancestral_jump: entries of type_no * type_no matrix>
<transit: entries of type_no * type_no matrix> //see also readMat
<length of replication histroty> <entries of replication history: vector> //see also readVec
<length of memory> <entries of memory: vector>
*/

Cells_Infinite read_cells_inifinite(std::ifstream &init_cells, std::ifstream &transition);
/*
init_cells
<type_no>
<type_1_fraction> <type_2_fraction> ... <type_n_fraction>

not necessary to be normalized

transition: same format as read_cell
*/
