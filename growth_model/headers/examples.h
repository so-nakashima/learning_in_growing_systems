#pragma once

#include <functional>
#include "simulators.h"
#include <vector>
#include "simulator_utility.h"
#include "analyze.h"
#include "linterp.h"
#include <iterator>
#include <algorithm>

void test_env();
void test_env_cells();
void file_read_test();
void test_analyze();
void lambda_curve();
void test_learning();
void graphic_test();     //obsolte; too slow and inaccurate
void graphic_infinite(); //use this instead
void lambda_sample();
void lambda_sample_highdim();
void lambda_sample_highdim_infinite();
void common_transition_test();
void check_fluctuating_relation_for_common_vs_ind();
double test_cells_infinite();
double test_cells_infinite_common();
void estimate_y_z_distribution();
void check_fluctuating_relation_for_common_vs_base();
void lambda_curve_infinite();
template <typename T>
void draw_learning_lineage_tree_from_simulation(std::string input_dir_relative_path, std::function<double(T)> out_func);
void test_cells_learn_commom();
void test_learning_common();
void compare_common_and_individual_learning();
void test_lineage_push();
void test_calc_lambda();

template <typename T>
void draw_learning_lineage_tree_from_simulation(std::string input_dir_relative_path, std::string output_dir_relative_path, std::function<double(T)> out_func, bool is_full_lineage_used = false)
{
    std::ifstream in_other(input_dir_relative_path + "//other.dat");
    std::map<std::string, std::string> parameters = read_parameters(in_other);
    const int endtime = std::stoi(parameters["end_time"]);
    const int type_no = std::stoi(parameters["type_no"]);
    const int mem_no = std::stoi(parameters["mem_no"]);

    std::ifstream in_env(input_dir_relative_path + "//res//env.dat");
    std::vector<int> environments;
    readVec(endtime, environments, in_env);

    std::ifstream in_pop(input_dir_relative_path + "//res//pop.dat");
    std::ifstream in_pop_full(input_dir_relative_path + "//res//pop_full.dat");
    Lineage<Cell_Learn> lineage = read_learning_lineage(type_no, mem_no, in_pop);
    Lineage<Cell_Learn> lineage_full = read_learning_lineage(type_no, mem_no, in_pop_full);

    //output lineage graph
    std::ofstream outgraph(output_dir_relative_path + "//res//graph.dot");
    std::ofstream outgraph_full(output_dir_relative_path + "//res//graph_full.dot");
    std::ofstream out_graph_max_min(output_dir_relative_path + "//res//graph_max_min.dat");
    std::ofstream out_graph_max_min_full(output_dir_relative_path + "//res//graph_max_min_full.dat");
    auto output_func = [](Cell_Learn c) {
        return c.transition[0][0] / (c.transition[0][0] + c.transition[0][1]);
    };
    if (!is_full_lineage_used)
        lineage.graphic(out_func, outgraph, out_graph_max_min);
    else
        lineage_full.graphic(out_func, outgraph_full, out_graph_max_min_full);
}