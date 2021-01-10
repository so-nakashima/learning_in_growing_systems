#include <iostream>
#include <vector>
#include <string>
#include "headers/examples.h"
#include <stdlib.h>

using namespace std;

int main()
{
   //file_read_test();
   //test_analyze();

   //test_learning();
   //lambda_curve_infinite();

   //graphic_infinite();
   //lambda_sample();

   //lambda_sample_highdim_infinite();

   //common_transition_test();
   //test_analyze();

   //check_fluctuating_relation_for_common_vs_ind();

   //test_cells_infinite_common();
   //estimate_y_z_distribution();

   //check_fluctuating_relation_for_common_vs_base();

   //test_learning_common();
   //test_lineage_push();

   printf("hote");

   compare_common_and_individual_learning();
   //test_calc_lambda();
   system("dot -Tpng .//experiments//sim_2_no_growth_comp//res//common//graph.dat -o .//experiments//sim_2_no_growth_comp//res//common//graph.png");
   system("dot -Tpng .//experiments//sim_2_no_growth_comp//res//learning/graph.dat -o .//experiments//sim_2_no_growth_comp//res//learning//graph.png");
   system("dot -Tpng .//experiments//sim_2_no_growth_comp//res//whole/graph.dat -o .//experiments//sim_2_no_growth_comp//res//whole//graph.png");
}
