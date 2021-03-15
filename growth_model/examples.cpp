#include "headers/simulators.h"
#include <vector>
#include "headers/simulator_utility.h"
#include "headers/analyze.h"
#include "headers/linterp.h"
#include <iterator>
#include <algorithm>
#include <stdlib.h>

void test_env()
{
    MBPRE w;

    w.set_end_time(30);

    std::vector<double> init = {1.0, 1.0};
    std::vector<std::vector<double>> tran = {{1.0, 0}, {0.5, 0.5}};
    Markov_Environments env(tran, 2, init);

    std::ofstream out_env(".//res//environments.dat");
    env.set_env_record(&out_env);
    w.set_environments(&env);

    w.excecute();
}

void test_env_cells()
{
    MBPRE w;

    w.set_end_time(20);

    //env
    std::vector<double> init = {1.0, 1.0};
    std::vector<std::vector<double>> env_tran = {{0.5, 0.5}, {0.5, 0.5}};
    Markov_Environments env(env_tran, 2, init);

    std::ofstream out_env(".//res//environments.dat");
    env.set_env_record(&out_env);
    w.set_environments(&env);

    //cells
    std::vector<std::vector<double>> type_tran = {{0.8, 0.2}, {0.2, 0.8}};
    std::vector<Cell> init_pop = {Cell(0, "0")}; //, Cell(0, "1")};

    Cells cells(2, type_tran, init_pop);
    std::ofstream out_pop(".//res//population.dat");
    std::ofstream out_pop_full(".//res//population_full.dat");
    cells.set_maximum_population_size(5);

    cells.set_pop_record(&out_pop);
    cells.set_pop_full_record(&out_pop_full);

    w.set_population(&cells);

    //replication
    std::vector<std::vector<std::vector<double>>> replication = {{{0.0, 0.0, 1.0}, {0.0, 0.0, 1.0}}, {{0.0, 0.0, 1.0}, {0.0, 0.0, 1.0}}};
    w.set_offspring_distributions(replication);

    w.excecute();
}

void file_read_test()
{

    //initialize world setting
    std::ifstream in_other(".//experiments//sim_2//other.dat");
    std::map<std::string, std::string> parameters = read_parameters(in_other);
    MBPRE w;
    w.set_end_time(std::stoi(parameters["end_time"]));

    std::ifstream in_env(".//experiments//sim_2//env.dat");
    Markov_Environments env = read_env(in_env);
    std::ifstream in_cells(".//experiments//sim_2//for_lambda//initial_cells.dat");
    std::ifstream in_cell_tran(".//experiments//sim_2//for_lambda//cell_type_tran.dat");
    Cells cells = read_cells(in_cells, in_cell_tran);
    cells.set_maximum_population_size(10);

    //replication
    std::vector<std::vector<std::vector<double>>> replication;
    std::ifstream in_repl(".//experiments//sim_2//replication.dat");
    read3DTensor<double>(replication, in_repl);
    w.set_offspring_distributions(replication);

    //record
    std::ofstream out_env(".//experiments//sim_2//res//env.dat");
    env.set_env_record(&out_env);
    w.set_environments(&env);
    w.set_population(&cells);
    std::ofstream out_pop(".//experiments//sim_2//res//pop.dat");
    std::ofstream out_pop_full(".//experiments//sim_2//res//pop_full.dat");
    cells.set_pop_record(&out_pop);
    cells.set_pop_full_record(&out_pop_full);

    w.excecute();
}

void test_analyze()
{
    std::ifstream in_other(".//experiments//sim_2//other.dat");
    std::map<std::string, std::string> parameters = read_parameters(in_other);
    const int endtime = std::stoi(parameters["end_time"]);

    std::ifstream in_env(".//experiments//sim_2//res//env.dat");
    std::vector<int> environments;
    readVec(endtime, environments, in_env);

    std::ifstream in_pop(".//experiments//sim_2//res//pop.dat");
    std::ifstream in_pop_full(".//experiments//sim_2//res//pop_full.dat");
    Lineage<Cell> lineage = read_lineage(in_pop);
    Lineage<Cell> lienage_full = read_lineage(in_pop_full);

    std::ofstream out_lambda(".//experiments//sim_2//res//analyze//lambda.dat", std::ios_base::app);
    out_lambda << lienage_full.lambda(10) << std::endl;

    std::ofstream out_retro(".//experiments//sim_2//res//analyze//retro.dat");
    std::function<double(Cell, Cell, int, int, int)> func = [](Cell p, Cell c, int g, int p_type, int c_type) { //retrospective type transition
        if (p.type() == p_type && c.type() == c_type)
        {
            return 1.0;
        }
        else
        {
            return 0.0;
        }

    };

    std::function<double(Cell, Cell, int, int, int, const std::vector<int> &)> func2 = [](Cell p, Cell c, int g, int c_type, int env, const std::vector<int> &envs) {
        if (c.type() == c_type && envs[g] == env)
        {
            return 1.0;
        }
        else
        {
            return 0.0;
        }
    };

    std::vector<std::function<double(Cell, Cell, int)>> funcs;

    for (int i = 0; i != 2; i++)
    {
        for (int j = 0; j != 2; j++)
        {
            using namespace std::placeholders;
            std::function<double(Cell, Cell, int)> temp = std::bind(func, _1, _2, _3, i, j);
            funcs.push_back(temp);
        }
    }

    for (int i = 0; i != 2; i++)
    {
        for (int j = 0; j != 2; j++)
        {
            using namespace std::placeholders;
            std::function<double(Cell, Cell, int)> temp = std::bind(func2, _1, _2, _3, j, i, environments);
            funcs.push_back(temp);
        }
    }

    std::vector<double> res = lineage.backward_mean(funcs);
    for (auto d : res)
    {
        out_retro << d << std::endl;
    }
}

void lambda_curve()
{
    std::ofstream out_lambda_curve(".//experiments//sim_1//res//analyze//lambda_curve_infinite.dat");

    int mean_no = 200;
    for (int t = 0; t != 100; t++)
    {
        double mean = 0.0;
        for (int i = 0; i != mean_no; i++)
        {
            //initialize world setting
            std::ifstream in_other(".//experiments//sim_1//world_other.dat");
            MBPRE w = read_mbpre(in_other);

            std::ifstream in_env(".//experiments//sim_1//env.dat");
            Markov_Environments env = read_env(in_env);
            std::ifstream in_cells(".//experiments//sim_1//initial_cells.dat");
            std::ifstream in_cell_tran(".//experiments//sim_1//cell_type_tran.dat");
            Cells cells = read_cells(in_cells, in_cell_tran);
            cells.set_maximum_population_size(100);
            std::vector<std::vector<double>> cell_tran = {{0.5, 0.5}, {0.5, 0.5}};
            cell_tran[0][0] = t * 0.01;
            cell_tran[0][1] = 1.0 - t * 0.01;
            cell_tran[1][0] = 1.0 - t * 0.01;
            cell_tran[1][1] = t * 0.01;
            cells.set_type_transition(cell_tran);

            //replication
            std::vector<std::vector<std::vector<double>>> replication;
            std::ifstream in_repl(".//experiments//sim_1//replication.dat");
            read3DTensor<double>(replication, in_repl);

            //replication[0][0][4] = 0.0;
            //replication[0][0][3] = .5 - t / 60.0;
            //replication[0][0][9] = .5 + t / 60.0;
            //replication[1][1][4] = 0.0;
            //replication[1][1][1] = .5 - t / 60.0;
            //replication[1][1][7] = .5 + t / 60.0;

            w.set_offspring_distributions(replication);
            w.set_end_time(100);

            //record
            std::ofstream out_env(".//experiments//sim_1//res//env.dat");
            env.set_env_record(&out_env);
            w.set_environments(&env);
            w.set_population(&cells);
            std::ofstream out_pop(".//experiments//sim_1//res//pop.dat");
            std::ofstream out_pop_full(".//experiments//sim_1//res//pop_full.dat");
            cells.set_pop_record(&out_pop);
            cells.set_pop_full_record(&out_pop_full);

            w.excecute();

            std::ifstream in_pop_full(".//experiments//sim_1//res//pop_full.dat");

            Lineage<Cell> lienage_full = read_lineage(in_pop_full);

            mean += lienage_full.lambda(100) / mean_no;
        }
        std::cout << t << "-th calculation end" << std::endl;
        out_lambda_curve << mean << std::endl;
    }
}

void test_learning()
{
    //initialize world setting
    std::ifstream in_other(".//experiments//sim_2//other.dat");
    std::map<std::string, std::string> parameters = read_parameters(in_other);
    MBPRE w;
    w.set_end_time(std::stoi(parameters["end_time"]));

    std::ifstream in_env(".//experiments//sim_2//env.dat");
    Markov_Environments env = read_env(in_env);

    //initalize population
    std::ifstream in_cells(".//experiments//sim_2//initial_cells.dat");
    Cells_Learn cells = read_cells_learn(in_cells);

    //define learning rule
    int type_no = cells.cardinality();
    double learning_rate = std::stod(parameters["learning_rate"]);

    /* auto learning_rule
    = [=](int p_type, int d_type, int no_daughters, std::vector<std::vector<double>>& tran, std::vector<std::vector<double>>& jump_hist, std::vector<double>& rep_hist, std::vector<double>& mem, std::mt19937_64& mt){

        //update ancestral jump
        for(int i = 0; i != type_no; i++){
            for(int j = 0; j != type_no; j++){
                jump_hist[i][j] = (1 - learning_rate) * jump_hist[i][j] + learning_rate * ((i == p_type && j == d_type)? 1.0 : 0.0);
            }
        }

        //update transition using ancestral jump
        tran = jump_hist;
    }; */
    auto learning_rule = [=](int p_type, int d_type, int no_daughters, std::vector<std::vector<double>> &tran, std::vector<std::vector<double>> &jump_hist, std::vector<double> &rep_hist, std::vector<double> &mem, std::mt19937_64 &mt) {
        //update ancestral jump
        for (int i = 0; i != type_no; i++)
        {
            for (int j = 0; j != type_no; j++)
            {
                jump_hist[i][j] = (1 - learning_rate) * jump_hist[i][j] + learning_rate * ((i == p_type && j == d_type) ? 1.0 : 0.0);
            }
        }

        std::vector<double> pi(3, 0.0); //ancestral type_distribution at one-point
        for (int i = 0; i != type_no; i++)
        {
            for (int j = 0; j != type_no; j++)
            {
                pi[i] += jump_hist[i][j];
            }
        }

        for (int i = 0; i != type_no; i++)
        {
            for (int j = 0; j != type_no; j++)
            {
                tran[i][j] = pi[j];
            }
        }
    };

    cells.set_learning_rule(learning_rule);

    //replication
    std::vector<std::vector<std::vector<double>>> replication;
    std::ifstream in_repl(".//experiments//sim_2//replication.dat");
    read3DTensor<double>(replication, in_repl);
    w.set_offspring_distributions(replication);

    //record
    std::ofstream out_env(".//experiments//sim_2//res//env.dat");
    env.set_env_record(&out_env);
    w.set_environments(&env);
    w.set_population(&cells);
    std::ofstream out_pop(".//experiments//sim_2//res//pop.dat");
    std::ofstream out_pop_full(".//experiments//sim_2//res//pop_full.dat");
    cells.set_pop_record(&out_pop);
    cells.set_pop_full_record(&out_pop_full);

    w.excecute();
}

std::vector<double> linspace(double first, double last, int len)
{
    std::vector<double> result(len);
    double step = (last - first) / (len - 1);
    for (int i = 0; i < len; i++)
    {
        result[i] = first + i * step;
    }
    return result;
}

void graphic_test()
{
    const int endtime = 5 + 1;
    const int type_no = 2;
    const int mem_no = 0;

    std::ifstream in_env(".//experiments//sim_2//res//env.dat");
    std::vector<int> environments;
    readVec(endtime, environments, in_env);

    std::ifstream in_pop(".//experiments//sim_2//res//pop.dat");
    std::ifstream in_pop_full(".//experiments//sim_2//res//pop_full.dat");
    Lineage<Cell_Learn> lineage = read_learning_lineage(type_no, in_pop);
    Lineage<Cell_Learn> lineage_full = read_learning_lineage(type_no, in_pop_full);

    //lambda
    const int mesh_size = 50;
    const int length = mesh_size + 1;
    std::vector<double> grid1 = linspace(0.0, 1.0, length);
    std::vector<double> grid2 = linspace(0.0, 1.0, length);

    std::ifstream in_lambda_sample(R"(./experiments/sim_2/res/analyze/lambda_sample_point8_symmetric.dat)");
    std::vector<double> f_values;
    readVec(length * length, f_values, in_lambda_sample);

    std::vector<std::vector<double>::iterator> grid_iter_list;
    grid_iter_list.push_back(grid1.begin());
    grid_iter_list.push_back(grid2.begin());

    boost::array<int, 2> grid_sizes;
    grid_sizes[0] = length;
    grid_sizes[1] = length;

    int num_elements = length * length;

    InterpMultilinear<2, double> interp_ML(grid_iter_list.begin(), grid_sizes.begin(), f_values.data(), f_values.data() + num_elements);

    auto out_lambda = [&interp_ML](Cell_Learn c) {
        double t00 = c.transition[0][0] / (c.transition[0][0] + c.transition[0][1]);
        double t11 = c.transition[1][1] / (c.transition[1][0] + c.transition[1][1]);

        boost::array<double, 2> args = {t00, t11};

        return interp_ML.interp(args.begin());
    };

    //output lineage graph
    std::ofstream outgraph(".//experiments//sim_2//res//graph.dot");
    std::ofstream outgraph_full(".//experiments//sim_2//res//graph_full.dot");
    std::ofstream out_graph_max_min(".//experiments//sim_2//res//graph_max_min.dat");
    std::ofstream out_graph_max_min_full(".//experiments//sim_2//res//graph_max_min_full.dat");
    auto output_func = [](Cell_Learn c) {
        return c.transition[0][0] / (c.transition[0][0] + c.transition[0][1]);
    };
    lineage.graphic(out_lambda, outgraph, out_graph_max_min);
    lineage_full.graphic(output_func, outgraph_full, out_graph_max_min_full);
}

void output_lambda_helper_next(std::vector<int> &itrs, const std::vector<std::vector<double>> &coordinates, std::ofstream &out)
{

    for (int i = 0; i != itrs.size(); i++)
    {
        if (itrs[i] < coordinates[i].size() - 1)
        {
            itrs[i]++;
            break;
        }
        else
        {
            itrs[i] = 0;
            out << std::endl;
        }
    }
}

void output_lambdas_helper(const std::vector<std::vector<double>> &coordinates, std::ofstream &out,
                           std::function<double(const std::vector<double> &)> func)
{
    const int dim = coordinates.size();
    std::vector<int> itrs(dim, 0);

    int grid_no = 1;
    for (int i = 0; i != dim; i++)
    {
        grid_no *= coordinates[i].size();
    }

    for (int i = 0; i != grid_no; i++)
    {
        std::vector<double> args;
        for (int j = 0; j != dim; j++)
        {
            args.push_back(coordinates[j][itrs[j]]);
        }
        out << func(args) << " ";
        output_lambda_helper_next(itrs, coordinates, out);

        std::cout << std::to_string(i) + "-th calculation ends" << std::endl;
    }
}

void lambda_sample()
{
    std::ofstream out_lambda_sample(".//experiments//sim_2//res//analyze//lambda_sample_point8_symmetric.dat");

    int mean_no = 50;
    int mesh_size = 50;
    for (int t1 = 0; t1 != mesh_size + 1; t1++)
    {
        for (int t2 = 0; t2 != mesh_size; t2++)
        {
            double mean = 0.0;
            for (int i = 0; i != mean_no + 1; i++)
            {
                //initialize world setting
                std::ifstream in_other(".//experiments//sim_2//world_other.dat");
                MBPRE w = read_mbpre(in_other);

                std::ifstream in_env(".//experiments//sim_2//env.dat");
                Markov_Environments env = read_env(in_env);
                std::ifstream in_cells(".//experiments//sim_2//for_lambda//initial_cells.dat");
                std::ifstream in_cell_tran(".//experiments//sim_2//for_lambda//cell_type_tran.dat");
                Cells cells = read_cells(in_cells, in_cell_tran);
                std::vector<std::vector<double>> cell_tran = {{0.5, 0.5}, {0.5, 0.5}};
                cell_tran[0][0] = ((double)t1) / mesh_size;
                cell_tran[0][1] = 1 - ((double)t1) / mesh_size;
                cell_tran[1][0] = 1 - ((double)t2) / mesh_size;
                cell_tran[1][1] = ((double)t2) / mesh_size;
                cells.set_type_transition(cell_tran);

                //replication
                std::vector<std::vector<std::vector<double>>> replication;
                std::ifstream in_repl(".//experiments//sim_2//replication.dat");
                read3DTensor<double>(replication, in_repl);

                //replication[0][0][4] = 0.0;
                //replication[0][0][3] = .5 - t / 60.0;
                //replication[0][0][9] = .5 + t / 60.0;
                //replication[1][1][4] = 0.0;
                //replication[1][1][1] = .5 - t / 60.0;
                //replication[1][1][7] = .5 + t / 60.0;

                w.set_offspring_distributions(replication);
                w.set_end_time(100);

                //record
                std::ofstream out_env(".//experiments//sim_2//res//env.dat");
                env.set_env_record(&out_env);
                w.set_environments(&env);
                w.set_population(&cells);
                std::ofstream out_pop(".//experiments//sim_2//res//pop.dat");
                std::ofstream out_pop_full(".//experiments//sim_2//res//pop_full.dat");
                cells.set_pop_record(&out_pop);
                cells.set_pop_full_record(&out_pop_full);

                w.excecute();

                std::ifstream in_pop_full(".//experiments//sim_2//res//pop_full.dat");

                Lineage<Cell> lienage_full = read_lineage(in_pop_full);

                mean += lienage_full.lambda(100) / mean_no;
            }
            std::cout << (mesh_size + 1) * t1 + t2 << "-th calculation end" << std::endl;
            out_lambda_sample << mean << " ";
        }
        out_lambda_sample << std::endl;
    }
}

void lambda_sample_highdim()
{
    const int mesh_size = 6;
    const int dim = 6;
    const int mean_no = 10;
    std::ofstream out_lambda_sample(".//experiments//sim_2//res//analyze//lambda_sample_cyclic.dat");

    std::vector<std::vector<double>> coordinates;
    for (int i = 0; i != dim; i++)
    {
        std::vector<double> temp;
        for (int j = 0; j != mesh_size; j++)
            temp.push_back(1.0 / (mesh_size - 1) * j);

        coordinates.push_back(temp);
    }

    auto func = [](const std::vector<double> &args) {
        double mean = 0.0;

        for (int i = 0; i != mean_no; i++)
        {
            std::ifstream in_other(".//experiments//sim_2//world_other.dat");
            MBPRE w = read_mbpre(in_other);

            std::ifstream in_env(".//experiments//sim_2//env.dat");
            Markov_Environments env = read_env(in_env);
            std::ifstream in_cells(".//experiments//sim_2//for_lambda//initial_cells.dat");
            std::ifstream in_cell_tran(".//experiments//sim_2//for_lambda//cell_type_tran.dat");
            Cells cells = read_cells(in_cells, in_cell_tran);
            std::vector<std::vector<double>> cell_tran = {{0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}, {0.5, 0.5, 0.5}};

            cell_tran[0][0] = args[0];
            cell_tran[0][1] = args[1];
            cell_tran[0][2] = 1 - args[0] - args[1];
            cell_tran[1][0] = args[2];
            cell_tran[1][1] = args[3];
            cell_tran[1][2] = 1 - args[2] - args[3];
            cell_tran[2][0] = args[4];
            cell_tran[2][1] = args[5];
            cell_tran[2][2] = 1 - args[4] - args[5];
            cells.set_type_transition(cell_tran);

            //replication
            std::vector<std::vector<std::vector<double>>> replication;
            std::ifstream in_repl(".//experiments//sim_2//replication.dat");
            read3DTensor<double>(replication, in_repl);

            w.set_offspring_distributions(replication);
            w.set_end_time(100);

            //record
            std::ofstream out_env(".//experiments//sim_2//res//env.dat");
            env.set_env_record(&out_env);
            w.set_environments(&env);
            w.set_population(&cells);
            std::ofstream out_pop(".//experiments//sim_2//res//pop.dat");
            std::ofstream out_pop_full(".//experiments//sim_2//res//pop_full.dat");
            cells.set_pop_record(&out_pop);
            cells.set_pop_full_record(&out_pop_full);

            w.excecute();

            std::ifstream in_pop_full(".//experiments//sim_2//res//pop_full.dat");

            Lineage<Cell> lienage_full = read_lineage(in_pop_full);

            mean += lienage_full.lambda(100) / mean_no;
        }

        return mean;
    };

    output_lambdas_helper(coordinates, out_lambda_sample, func);
}

void common_transition_test()
{

    //initialize world setting
    std::ifstream in_other(".//experiments//sim_2//other.dat");
    std::map<std::string, std::string> parameters = read_parameters(in_other);
    MBPRE w;
    w.set_end_time(std::stoi(parameters["end_time"]));

    std::ifstream in_env(".//experiments//sim_2//env.dat");
    Markov_Environments env = read_env(in_env);
    std::ifstream in_cells(".//experiments//sim_2//for_lambda//initial_cells.dat");
    std::ifstream in_cell_tran(".//experiments//sim_2//for_lambda//cell_type_tran.dat");
    Cells_Common cells(read_cells(in_cells, in_cell_tran));
    cells.set_maximum_population_size(10);
    //cells.set_common_p_type(std::stoi(parameters["initial_p_type"]));

    //replication
    std::vector<std::vector<std::vector<double>>> replication;
    std::ifstream in_repl(".//experiments//sim_2//replication.dat");
    read3DTensor<double>(replication, in_repl);
    w.set_offspring_distributions(replication);

    //record
    std::ofstream out_env(".//experiments//sim_2//res//env.dat");
    env.set_env_record(&out_env);
    std::ofstream out_pop(".//experiments//sim_2//res//pop.dat");
    std::ofstream out_pop_full(".//experiments//sim_2//res//pop_full.dat");
    std::ofstream out_shared_p_type(".//experiments//sim_2//res//shared_p_type.dat");
    cells.set_pop_record(&out_pop);
    cells.set_pop_full_record(&out_pop_full);
    cells.set_out_shared_p_type(&out_shared_p_type);

    w.set_environments(&env);
    w.set_population(&cells);
    w.excecute();
}

/*void check_fluctuating_relation_for_common_vs_ind(){
    //initialize world setting
    std::ifstream in_other(".//experiments//sim_2//other.dat");
    std::map<std::string, std::string> parameters = read_parameters(in_other);
    const int avg_no = std::stoi(parameters["fluctuating_avg_no"]);

    const int endtime = std::stoi(parameters["end_time"]);


    std::vector<std::vector<std::vector<double>>> replication;
    std::ifstream in_repl(".//experiments//sim_2//replication.dat");
    read3DTensor<double>(replication, in_repl);


    std::ofstream out_lambda(".//experiments//sim_2//res//analyze//lambda.dat");

    std::ofstream out_res(".//experiments//sim_2//res//analyze//fluctuation.dat");
    
    double res = 0.0;
    for(int i = 0; i != avg_no; i++){
        //generate common env. sequence
        std::ifstream in_env(".//experiments//sim_2//env.dat");
        Markov_Environments env_temp = read_env(in_env);
        std::vector<int> env_seq = env_temp.generate(endtime);

        Environments_Sequence env(3, env_seq); 


        //compute lambda for individual transition
        //initialize world setting
        
        MBPRE w;
        w.set_end_time(std::stoi(parameters["end_time"]));
        std::ifstream in_other(".//experiments//sim_2//other.dat");
        std::map<std::string, std::string> parameters = read_parameters(in_other);
        w.set_end_time(std::stoi(parameters["end_time"]));



        std::ifstream in_cells(".//experiments//sim_2//for_lambda//initial_cells.dat");
        std::ifstream in_cell_tran(".//experiments//sim_2//for_lambda//cell_type_tran.dat");
        Cells cells = read_cells(in_cells, in_cell_tran);
        cells.set_maximum_population_size(1000);


        //replication
        w.set_offspring_distributions(replication);

        //record
        std::ofstream out_pop(".//experiments//sim_2//res//pop.dat");
        std::ofstream out_pop_full(".//experiments//sim_2//res//pop_full.dat");
        cells.set_pop_record(&out_pop);
        cells.set_pop_full_record(&out_pop_full);


        w.set_environments(&env);
        w.set_population(&cells);
        w.excecute();

        std::ifstream in_env_rec(".//experiments//sim_2//res//env.dat");
        std::vector<int> environments;
        readVec(endtime, environments, in_env_rec);

        std::ifstream in_pop(".//experiments//sim_2//res//pop.dat");
        std::ifstream in_pop_full(".//experiments//sim_2//res//pop_full.dat");
        Lineage<Cell> lineage = read_lineage(in_pop);
        Lineage<Cell> lienage_full = read_lineage(in_pop_full);

        double lambda_ind = lienage_full.lambda(1000);
        out_lambda << lambda_ind << std::endl;
        
        //compute lambda for common transition
        MBPRE w_;
        w_.set_end_time(std::stoi(parameters["end_time"]));


        std::ifstream in_cells_(".//experiments//sim_2//for_lambda//initial_cells.dat");
        std::ifstream in_cell_tran_(".//experiments//sim_2//for_lambda//cell_type_tran.dat");
        Cells_Common cells_(read_cells(in_cells_, in_cell_tran_));
        //cells.set_common_p_type(std::stoi(parameters["initial_p_type"]));


        //replication
        w_.set_offspring_distributions(replication);

        //record
        std::ofstream out_pop_(".//experiments//sim_2//res//pop.dat");
        std::ofstream out_pop_full_(".//experiments//sim_2//res//pop_full.dat");
        std::ofstream out_shared_p_type_(".//experiments//sim_2//res//shared_p_type.dat");
        cells_.set_pop_record(&out_pop_);
        cells_.set_pop_full_record(&out_pop_full_);
        cells_.set_out_shared_p_type(&out_shared_p_type_);


        w_.set_environments(&env);
        w_.set_population(&cells_);
        w_.excecute();


        std::ifstream in_env_rec_(".//experiments//sim_2//res//env.dat");
        std::vector<int> environments_;
        readVec(endtime, environments_, in_env_rec_);

        std::ifstream in_pop_full_(".//experiments//sim_2//res//pop_full.dat");
        Lineage<Cell> lienage_full_ = read_lineage(in_pop_full_);


        double lambda_common =  lienage_full_.lambda(1000);
        out_lambda << lambda_common << std::endl;

        //compute average
        res += (exp((- lambda_ind + lambda_common) * (endtime-1))) / avg_no;
    }

    out_res << "mean: " << res << std::endl;
}*/

double test_cells_infinite()
{
    //initialize world setting
    std::ifstream in_other(".//experiments//sim_3//other.dat");
    std::map<std::string, std::string> parameters = read_parameters(in_other);
    MBPRE w;
    w.set_end_time(std::stoi(parameters["end_time"]));

    std::ifstream in_env(".//experiments//sim_3//env.dat");
    Markov_Environments env = read_env(in_env);

    std::ifstream in_cells(".//experiments//sim_3//initial_cells.dat");
    std::ifstream in_cell_tran(".//experiments//sim_3//cell_type_tran.dat");
    Cells_Infinite cells(read_cells_inifinite(in_cells, in_cell_tran));

    //replication
    std::vector<std::vector<std::vector<double>>> replication;
    std::ifstream in_repl(".//experiments//sim_3//replication.dat");
    read3DTensor<double>(replication, in_repl);
    w.set_offspring_distributions(replication);

    //record
    std::ofstream out_env(".//experiments//sim_3//res//env.dat");
    env.set_env_record(&out_env);
    std::ofstream out_pop(".//experiments//sim_3//res//pop.dat");
    cells.set_output(&out_pop);

    w.set_environments(&env);
    w.set_population(&cells);
    w.excecute();

    double res = 0.0;
    for (auto l : cells.lambdas)
    {
        res += l;
    }
    return res;
}

double test_cells_infinite_common()
{
    //initialize world setting
    std::ifstream in_other(".//experiments//sim_3//other.dat");
    std::map<std::string, std::string> parameters = read_parameters(in_other);
    MBPRE w;
    w.set_end_time(std::stoi(parameters["end_time"]));

    std::ifstream in_env(".//experiments//sim_3//env.dat");
    Markov_Environments env = read_env(in_env);

    std::ifstream in_cells(".//experiments//sim_3//initial_cells.dat");
    std::ifstream in_cell_tran(".//experiments//sim_3//cell_type_tran.dat");
    Cells_Infinite_Common cells(read_cells_inifinite(in_cells, in_cell_tran));

    //replication
    std::vector<std::vector<std::vector<double>>> replication;
    std::ifstream in_repl(".//experiments//sim_3//replication.dat");
    read3DTensor<double>(replication, in_repl);
    w.set_offspring_distributions(replication);

    //record
    std::ofstream out_env(".//experiments//sim_3//res//env.dat");
    env.set_env_record(&out_env);
    std::ofstream out_pop(".//experiments//sim_3//res//pop.dat");
    cells.set_output(&out_pop);

    w.set_environments(&env);
    w.set_population(&cells);
    w.excecute();

    double res = 0.0;
    for (auto l : cells.lambdas)
    {
        res += l;
    }
    return res;
}

double lambda_cells_infinite(int env_no, const std::vector<int> envs)
{
    //initialize world setting
    std::ifstream in_other(".//experiments//sim_3//other.dat");
    std::map<std::string, std::string> parameters = read_parameters(in_other);
    MBPRE w;
    w.set_end_time(std::stoi(parameters["end_time"]));

    Environments_Sequence env(env_no, envs);

    std::ifstream in_cells(".//experiments//sim_3//initial_cells.dat");
    std::ifstream in_cell_tran(".//experiments//sim_3//cell_type_tran.dat");
    Cells_Infinite cells(read_cells_inifinite(in_cells, in_cell_tran));

    //replication
    std::vector<std::vector<std::vector<double>>> replication;
    std::ifstream in_repl(".//experiments//sim_3//replication.dat");
    read3DTensor<double>(replication, in_repl);
    w.set_offspring_distributions(replication);

    //record
    std::ofstream out_pop(".//experiments//sim_3//res//pop.dat");
    cells.set_output(&out_pop);

    w.set_environments(&env);
    w.set_population(&cells);
    w.excecute();

    double res = 0.0;
    for (auto l : cells.lambdas)
    {
        res += l;
    }
    return res;
}

double lambda_cells_infinite_common(int env_no, const std::vector<int> envs)
{
    //initialize world setting
    std::ifstream in_other(".//experiments//sim_3//other.dat");
    std::map<std::string, std::string> parameters = read_parameters(in_other);
    MBPRE w;
    w.set_end_time(std::stoi(parameters["end_time"]));

    Environments_Sequence env(env_no, envs);

    std::ifstream in_cells(".//experiments//sim_3//initial_cells.dat");
    std::ifstream in_cell_tran(".//experiments//sim_3//cell_type_tran.dat");
    Cells_Infinite_Common cells(read_cells_inifinite(in_cells, in_cell_tran));

    //replication
    std::vector<std::vector<std::vector<double>>> replication;
    std::ifstream in_repl(".//experiments//sim_3//replication.dat");
    read3DTensor<double>(replication, in_repl);
    w.set_offspring_distributions(replication);

    //record
    std::ofstream out_pop(".//experiments//sim_3//res//pop.dat");
    cells.set_output(&out_pop);

    w.set_environments(&env);
    w.set_population(&cells);
    w.excecute();

    double res = 0.0;
    for (auto l : cells.lambdas)
    {
        res += l;
    }
    return res;
}

double lambda_cells_infinite_common(int env_no, const std::vector<int> envs, std::vector<int> &vec)
{
    //initialize world setting
    std::ifstream in_other(".//experiments//sim_3//other.dat");
    std::map<std::string, std::string> parameters = read_parameters(in_other);
    MBPRE w;
    w.set_end_time(std::stoi(parameters["end_time"]));

    Environments_Sequence env(env_no, envs);

    std::ifstream in_cells(".//experiments//sim_3//initial_cells.dat");
    std::ifstream in_cell_tran(".//experiments//sim_3//cell_type_tran.dat");
    Cells_Infinite_Common cells(read_cells_inifinite(in_cells, in_cell_tran));

    cells.current_pop[0] = (0.22118 * 0.8 + 0.05415 * 0.2) * 3.0;
    cells.current_pop[1] = (0.05595 * 0.8 + 0.22118 * 0.2) * 3.0;
    cells.current_pop[2] = (0.05415 * 0.8 + 0.05595 * 0.2) * 3.0;

    //replication
    std::vector<std::vector<std::vector<double>>> replication;
    std::ifstream in_repl(".//experiments//sim_3//replication.dat");
    read3DTensor<double>(replication, in_repl);
    w.set_offspring_distributions(replication);

    //record
    std::ofstream out_pop(".//experiments//sim_3//res//pop.dat");
    cells.set_output(&out_pop);

    w.set_environments(&env);
    w.set_population(&cells);
    w.excecute();

    vec = cells.hist_common_p_type;

    double res = 0.0;
    for (auto l : cells.lambdas)
    {
        res += l;
    }
    return res;
}

void check_fluctuating_relation_for_common_vs_ind()
{
    std::ifstream in_other(".//experiments//sim_3//other.dat");
    std::map<std::string, std::string> parameters = read_parameters(in_other);
    const int avg_no = std::stoi(parameters["fluctuating_avg_no"]);
    const int end_time = std::stoi(parameters["end_time"]);

    std::ofstream out_lambda(".//experiments//sim_3//res//analyze//lambda.dat");
    std::ofstream out_res(".//experiments//sim_3//res//analyze//fluctuation.dat");

    double res = 0.0;
    for (int i = 0; i != avg_no; i++)
    {
        std::ifstream in_env(".//experiments//sim_3//env.dat");
        Markov_Environments env = read_env(in_env);
        std::vector<int> envs = env.generate(end_time + 1);

        double lambda_ind = lambda_cells_infinite(env.cardinality(), envs);
        double lambda_common = lambda_cells_infinite_common(env.cardinality(), envs);
        out_lambda << lambda_ind << " " << lambda_common << std::endl;
        res += exp(lambda_common - lambda_ind) / avg_no;
    }

    out_res << res << std::endl;
}

void estimate_y_z_distribution()
{
    std::ifstream in_other(".//experiments//sim_3//other.dat");
    std::map<std::string, std::string> parameters = read_parameters(in_other);
    const int end_time = std::stoi(parameters["end_time"]);
    const int env_type_no = std::stoi(parameters["env_type_no"]);
    const int z_type_no = std::stoi(parameters["z_type_no"]);

    std::ifstream in_pop(".//experiments//sim_3//res//pop.dat");
    std::vector<std::vector<double>> pop_data;
    readMat(end_time + 1, z_type_no + 2, pop_data, in_pop);

    std::vector<std::vector<double>> res(env_type_no, std::vector<double>(z_type_no));
    std::ofstream out_y_z(".//experiments//sim_3//res//analyze//y-z.dat");

    std::ifstream in_env(".//experiments//sim_3//res//env.dat");
    std::vector<int> env;
    readVec(end_time, env, in_env);

    const double inv_end_time = 1.0 / end_time;
    for (int i = 0; i != end_time; i++)
    {
        const int c_env = env[i];                //current_env
        const int c_p_type = pop_data[i + 1][1]; //current_common_p_type

        res[c_env][c_p_type] += inv_end_time;
    }

    out_mat(res, &out_y_z);
}

double lambda_base(int end_time, std::vector<int> envs, std::vector<int> z_seq, const std::vector<std::vector<double>> &y_z_dist, const std::vector<double> &y_dist, const std::vector<double> &z_dist)
{
    double res = 0.0;
    for (int t = 0; t != end_time; t++)
    {
        res += log(2) + log(y_z_dist[envs[t]][z_seq[t]]) - log(y_dist[envs[t]]) - log(z_dist[z_seq[t]]);
    }
    return res;
}

void check_fluctuating_relation_for_common_vs_base()
{
    std::ifstream in_other(".//experiments//sim_3//other.dat");
    std::map<std::string, std::string> parameters = read_parameters(in_other);
    const int avg_no = std::stoi(parameters["fluctuating_avg_no"]);
    const int end_time = std::stoi(parameters["end_time"]);
    const int env_type_no = std::stoi(parameters["env_type_no"]);
    const int z_type_no = std::stoi(parameters["z_type_no"]);

    std::ofstream out_lambda(".//experiments//sim_3//res//analyze//lambda.dat");
    std::ofstream out_res(".//experiments//sim_3//res//analyze//fluctuation_common_base.dat");

    std::ifstream in_y_z(".//experiments//sim_3//res//analyze//y-z_100000.dat");
    std::vector<std::vector<double>> y_z_dist;
    readMat(env_type_no, z_type_no, y_z_dist, in_y_z);
    std::vector<double> y_dist(env_type_no, 0.0), z_dist(z_type_no, 0.0);
    for (int i = 0; i != env_type_no; i++)
    {
        for (int j = 0; j != z_type_no; j++)
        {
            y_dist[i] += y_z_dist[i][j];
            z_dist[j] += y_z_dist[i][j];
        }
    }

    double res = 0.0;
    for (int i = 0; i != avg_no; i++)
    {
        std::ifstream in_env(".//experiments//sim_3//env.dat");
        Markov_Environments env = read_env(in_env);
        std::vector<int> envs = env.generate(end_time + 1);

        std::vector<int> z_seq;

        //double lambda_ind = lambda_cells_infinite(env.cardinality(), envs);
        double lambda_common = lambda_cells_infinite_common(env.cardinality(), envs, z_seq);
        double lambda_other = lambda_base(end_time, envs, z_seq, y_z_dist, y_dist, z_dist);
        out_lambda << lambda_common << " " << lambda_other << std::endl;
        res += exp(lambda_common - lambda_other) / avg_no;
    }

    out_res << res << std::endl;
}

void sim_learning(std::string setting_dir_rel_path, std::string output_dir_rel_path, std::function<void(int p_type, int d_type, int no_daughters, std::vector<std::vector<double>> &tran, std::vector<std::vector<double>> &jump_hist, std::vector<double> &rep_hist, std::vector<double> &mem, std::mt19937_64 &mt)> learning_rule, bool enable_common_learning = false, std::vector<std::vector<double>> transition = std::vector<std::vector<double>>())
{
    //initialize world setting
    std::ifstream in_other(setting_dir_rel_path + "//other.dat");
    std::map<std::string, std::string> parameters = read_parameters(in_other);
    MBPRE w;
    w.set_end_time(std::stoi(parameters["end_time"]));

    std::ifstream in_env(setting_dir_rel_path + "//env.dat");
    Markov_Environments env = read_env(in_env);

    //initalize population
    std::ifstream in_cells(setting_dir_rel_path + "//initial_cells.dat");
    Cells_Learn *cells;
    //std::ofstream out_spine(output_dir_rel_path + "//spine.dat"); //only used for common learning
    if (enable_common_learning)
    {
        Cells_Learn_Common *cells_temp = new_cells_learn_common_from_read(in_cells);
        // cells_temp->set_out_spine(&out_spine);
        cells = cells_temp;
    }
    else
    {
        cells = new_cells_learn_from_read(in_cells);
    }
    cells->set_learning_rule(learning_rule);

    //replication
    std::vector<std::vector<std::vector<double>>> replication;
    std::ifstream in_repl(setting_dir_rel_path + "//replication.dat");
    read3DTensor<double>(replication, in_repl);
    w.set_offspring_distributions(replication);

    //record
    std::ofstream out_env(output_dir_rel_path + "//env_res.dat");
    env.set_env_record(&out_env);
    w.set_environments(&env);
    w.set_population(cells);
    std::ofstream out_pop(output_dir_rel_path + "//pop.dat");
    std::ofstream out_pop_full(output_dir_rel_path + "//pop_full.dat");
    cells->set_pop_record(&out_pop);
    cells->set_pop_full_record(&out_pop_full);
    std::ofstream *Pout_spine;
    if (enable_common_learning)
    {
        Pout_spine = new std::ofstream(output_dir_rel_path + "//spine.dat");
        ((Cells_Learn_Common *)cells)->set_out_spine(Pout_spine);
    }

    w.excecute();

    delete cells;
    if (enable_common_learning)
    {
        delete Pout_spine;
    }
}

void test_learning_common()
{
    //set directories
    const std::string setting_rel_path = ".//experiments//sim_2_no_growth_comp";
    const std::string output_rel_path_4_common = ".//experiments//sim_2_no_growth_comp//res//common";
    const std::string output_rel_path_4_individual = ".//experiments//sim_2_no_growth_comp//res//learning";

    //read paramers and define const variables
    std::ifstream in_other(setting_rel_path + "//other.dat");
    std::map<std::string, std::string> parameters = read_parameters(in_other);
    const int type_no = std::stoi(parameters["type_no"]);
    const double learning_rate = std::stod(parameters["learning_rate"]);

    auto learning_rule = [=](int p_type, int d_type, int no_daughters, std::vector<std::vector<double>> &tran, std::vector<std::vector<double>> &jump_hist, std::vector<double> &rep_hist, std::vector<double> &mem, std::mt19937_64 &mt) {
        //update ancestral jump
        for (int i = 0; i != type_no; i++)
        {
            for (int j = 0; j != type_no; j++)
            {
                jump_hist[i][j] = (1 - learning_rate) * jump_hist[i][j] + learning_rate * ((i == p_type && j == d_type) ? 1.0 : 0.0);
            }
        }

        //update transition using ancestral jump
        tran = jump_hist;
    };

    sim_learning(setting_rel_path, output_rel_path_4_common, learning_rule, true);
    //sim_learning(setting_rel_path, output_rel_path_4_individual, learning_rule, false);
}

Lineage<Cell_Learn> generate_lineage_learning(std::string setting_dir_rel_path, std::string output_dir_rel_path, std::function<void(int p_type, int d_type, int no_daughters, std::vector<std::vector<double>> &tran, std::vector<std::vector<double>> &jump_hist, std::vector<double> &rep_hist, std::vector<double> &mem, std::mt19937_64 &mt)> learning_rule, bool enable_common_learning = false, std::string lineage_file_name = "pop.dat")
{
    //execute simulation
    sim_learning(setting_dir_rel_path, output_dir_rel_path, learning_rule, enable_common_learning);

    //construct lineage
    std::ifstream in_pop(output_dir_rel_path + "//" + lineage_file_name);

    std::ifstream in_other(setting_dir_rel_path + "//other.dat");
    std::map<std::string, std::string> parameters = read_parameters(in_other);
    const int mem_no = std::stoi(parameters["mem_no"]);
    const int type_no = std::stoi(parameters["type_no"]);

    return read_learning_lineage(type_no, in_pop);
}

std::vector<double> sim_no_learning_inf(const std::string &setting_dir_rel_path, const std::string &out_dir_rel_path, const std::vector<std::vector<double>> &transition = {}) //can change transition matrix by the last arg; return cells.lambdas
{
    //initialize world setting
    std::ifstream in_other(setting_dir_rel_path + "//other.dat");
    std::map<std::string, std::string> parameters = read_parameters(in_other);
    MBPRE w;
    w.set_end_time(std::stoi(parameters["end_time"]));

    std::ifstream in_env(setting_dir_rel_path + "//env.dat");
    Markov_Environments env = read_env(in_env);

    //cell setting
    std::ifstream in_cells(setting_dir_rel_path + "//initial_cells.dat");
    std::ifstream in_cell_tran(".//experiments//sim_3//cell_type_tran.dat");
    Cells_Infinite cells(read_cells_inifinite(in_cells, in_cell_tran));
    if (transition.size() != 0)
    {
        cells.set_type_transition(transition);
    }

    //replication
    std::vector<std::vector<std::vector<double>>> replication;
    std::ifstream in_repl(setting_dir_rel_path + "//replication.dat");
    read3DTensor<double>(replication, in_repl);
    w.set_offspring_distributions(replication);

    //record
    std::ofstream out_pop(out_dir_rel_path + " //pop.dat");
    cells.set_output(&out_pop);

    w.set_environments(&env);
    w.set_population(&cells);
    w.excecute();

    return cells.lambdas;
}

double calc_lambda(const Cell_Learn &cell, const std::string &setting_dir_rel_path, const std::string &out_dir_rel_path)
{
    //execute simulation and acquire lambdas
    std::vector<double> lambdas = sim_no_learning_inf(setting_dir_rel_path, out_dir_rel_path, cell.transition);

    double res = 0.0;
    for (auto l : lambdas)
    {
        res += l;
    }
    return res / lambdas.size();
}

void execute_graphviz(std::string input_rel_path, std::string output_rel_path, std::string extension = "png")
{
    std::string command = "dot -T" + extension + " " + input_rel_path + " -o " + output_rel_path;
    int is_failed = system(command.c_str());
    if (is_failed)
    {
        std::cout << "Graph drawing failed at execute_graphviz with command: " + command << std::endl;
    }
}

void compare_common_and_individual_learning()
{
    //set directories (modify if necessary)
    //for generating lineage
    const std::string setting_rel_path = ".//experiments//sim_2_no_growth_comp";
    const std::string output_rel_path_4_common = ".//experiments//sim_2_no_growth_comp//res//common";
    const std::string output_rel_path_4_individual = ".//experiments//sim_2_no_growth_comp//res//learning";

    //for calculating lambda (in order to pass lienage.graphic)
    const std::string setting_calc_lambda_dir_rel_path = ".//experiments//sim_2_no_growth_comp//calc_lambdas";
    const std::string out_calc_lambda_dir_rel_path = ".//experiments//sim_2_no_growth_comp//calc_lambdas//res";

    //for drawing graph
    std::ofstream out_ind_learning_lineage(".//experiments//sim_2_no_growth_comp//res//learning//graph.dat");
    std::ofstream out_common_learning(".//experiments//sim_2_no_growth_comp//res//common//graph.dat");
    std::ofstream out_whole_lineage(".//experiments//sim_2_no_growth_comp//res//whole//graph.dat");

    //not necessary, just complete arg. of lienage.grphic
    std::ofstream out_ind_learning_lineage_max_min(".//experiments//sim_2_no_growth_comp//res//learning//max_min.dat");
    std::ofstream out_common_learning_max_min(".//experiments//sim_2_no_growth_comp//res//common//max_min.dat");
    std::ofstream out_whole_max_min(".//experiments//sim_2_no_growth_comp//res//whole//max_min.dat");

    //read paramers to define const variables
    std::ifstream in_other(setting_rel_path + "//other.dat");
    std::map<std::string, std::string> parameters = read_parameters(in_other);

    //set learning rule to online EM
    const int type_no = std::stoi(parameters["type_no"]);
    const double learning_rate = std::stod(parameters["learning_rate"]);

    auto learning_rule = [=](int p_type, int d_type, int no_daughters, std::vector<std::vector<double>> &tran, std::vector<std::vector<double>> &jump_hist, std::vector<double> &rep_hist, std::vector<double> &mem, std::mt19937_64 &mt) {
        //update ancestral jump
        for (int i = 0; i != type_no; i++)
        {
            for (int j = 0; j != type_no; j++)
            {
                jump_hist[i][j] = (1 - learning_rate) * jump_hist[i][j] + learning_rate * ((i == p_type && j == d_type) ? 1.0 : 0.0);
            }
        }

        //update transition using ancestral jump
        tran = jump_hist;
    };

    //execute simulation and get lineage
    Lineage<Cell_Learn> lineage_ind = generate_lineage_learning(setting_rel_path, output_rel_path_4_individual, learning_rule, /* enable_common_learning */ false);

    //collect the same number of spines as max_pop_no
    Lineage<Cell_Learn> lineage_common;
    const int max_cell_no = std::stoi(parameters["max_cell_no"]);
    for (int i = 0; i != max_cell_no; i++)
    {
        Lineage<Cell_Learn> spine_lin = generate_lineage_learning(setting_rel_path, output_rel_path_4_common, learning_rule, /* enable_common_learning */ true, /*lineage_file_name*/ "spine.dat");
        lineage_common.push(spine_lin);
    }

    //whole lineage
    Lineage<Cell_Learn> lineage_whole;
    lineage_whole.push(lineage_common);
    lineage_whole.push(lineage_ind);

    //pass for lineage.graphic to calculate lambda of each cell
    auto cell2lambda = [&](const Cell_Learn &cell) {
        return calc_lambda(cell, setting_calc_lambda_dir_rel_path, out_calc_lambda_dir_rel_path);
    };

    //determine the range of plot
    double max_val = std::max(lineage_ind.max(cell2lambda), lineage_common.max(cell2lambda));
    double min_val = std::min(lineage_ind.min(cell2lambda), lineage_common.min(cell2lambda));

    lineage_ind.graphic(cell2lambda, min_val, max_val, out_ind_learning_lineage, out_ind_learning_lineage_max_min);
    lineage_common.graphic(cell2lambda, min_val, max_val, out_common_learning, out_common_learning_max_min);
    lineage_whole.graphic(cell2lambda, out_whole_lineage, out_whole_max_min);

    execute_graphviz(".//experiments//sim_2_no_growth_comp//res//learning//graph.dat", ".//experiments//sim_2_no_growth_comp//res//learning//graph.png");
    execute_graphviz(".//experiments//sim_2_no_growth_comp//res//common//graph.dat", ".//experiments//sim_2_no_growth_comp//res//common//graph.png");
    execute_graphviz(".//experiments//sim_2_no_growth_comp//res//whole//graph.dat", ".//experiments//sim_2_no_growth_comp//res//whole//graph.png");
}

void test_lineage_push()
{
    std::ifstream in_other(".//experiments//sim_2_no_growth_comp//other.dat");
    std::map<std::string, std::string> parameters = read_parameters(in_other);
    const int endtime = std::stoi(parameters["end_time"]);
    const int mem_no = std::stoi(parameters["mem_no"]);
    const int type_no = std::stoi(parameters["type_no"]);

    std::ifstream in_env(".//experiments//sim_2_no_growth_comp//res//learning//env_res.dat");
    std::vector<int> environments;
    readVec(endtime, environments, in_env);

    std::ifstream in_pop(".//experiments//sim_2_no_growth_comp//res//learning//pop.dat");
    std::ifstream in_pop_full(".//experiments//sim_2_no_growth_comp//res//learning//pop_full.dat");
    Lineage<Cell_Learn> lineage = read_learning_lineage(type_no, in_pop);
    Lineage<Cell_Learn> lienage_orginal = read_learning_lineage(type_no, in_pop);

    lineage.push(lineage);

    lineage.push(lineage);
}

void test_calc_lambda()
{
    //directories
    const std::string setting_calc_lambda_dir_rel_path = ".//experiments//sim_2_no_growth_comp//calc_lambdas";
    const std::string out_calc_lambda_dir_rel_path = ".//experiments//sim_2_no_growth_comp//calc_lambdas//res";
    //lineage location and parameters
    const std::string setting_dir_rel_path = ".//experiments//sim_2_no_growth_comp";
    const std::string output_dir_rel_path = ".//experiments//sim_2_no_growth_comp//res//learning";
    //for drawing graph
    std::ofstream out_ind_learning_lineage(".//experiments//sim_2_no_growth_comp//res//learning//graph.dat");
    std::ofstream out_ind_learning_lineage_max_min(".//experiments//sim_2_no_growth_comp//res//learning//max_min.dat");

    //pass for lineage.graphic to calculate lambda of each cell
    auto cell2lambda = [&](const Cell_Learn &cell) {
        return calc_lambda(cell, setting_calc_lambda_dir_rel_path, out_calc_lambda_dir_rel_path);
    };

    //construct lineage
    std::ifstream in_pop(output_dir_rel_path + "//pop.dat");

    std::ifstream in_other(setting_dir_rel_path + "//other.dat");
    std::map<std::string, std::string> parameters = read_parameters(in_other);
    const int mem_no = std::stoi(parameters["mem_no"]);
    const int type_no = std::stoi(parameters["type_no"]);

    Lineage<Cell_Learn> lineage = read_learning_lineage(type_no, in_pop);

    lineage.graphic(cell2lambda, out_ind_learning_lineage, out_ind_learning_lineage_max_min);
}

void no_spine_learning()
{
    //set directories (modify if necessary)
    //for generating lineage
    const std::string setting_rel_path = ".//experiments//sim_2_no_growth_comp";
    const std::string output_rel_path_4_common = ".//experiments//sim_2_no_growth_comp//res//common";

    //for calculating lambda (in order to pass lienage.graphic)
    const std::string setting_calc_lambda_dir_rel_path = ".//experiments//sim_2_no_growth_comp//calc_lambdas";
    const std::string out_calc_lambda_dir_rel_path = ".//experiments//sim_2_no_growth_comp//calc_lambdas//res";

    //for drawing graph
    std::ofstream out_common_learning(".//experiments//sim_2_no_growth_comp//res//common//graph.dat");

    //not necessary, just complete arg. of lienage.grphic
    std::ofstream out_common_learning_max_min(".//experiments//sim_2_no_growth_comp//res//common//max_min.dat");

    //read paramers to define const variables
    std::ifstream in_other(setting_rel_path + "//other.dat");
    std::map<std::string, std::string> parameters = read_parameters(in_other);

    //set learning rule to online EM
    const int type_no = std::stoi(parameters["type_no"]);
    const double learning_rate = std::stod(parameters["learning_rate"]);

    auto learning_rule = [=](int p_type, int d_type, int no_daughters, std::vector<std::vector<double>> &tran, std::vector<std::vector<double>> &jump_hist, std::vector<double> &rep_hist, std::vector<double> &mem, std::mt19937_64 &mt) {
        //update ancestral jump
        for (int i = 0; i != type_no; i++)
        {
            for (int j = 0; j != type_no; j++)
            {
                jump_hist[i][j] = (1 - learning_rate) * jump_hist[i][j] + learning_rate * ((i == p_type && j == d_type) ? 1.0 : 0.0);
            }
        }

        //update transition using ancestral jump
        tran = jump_hist;
    };

    //collect the same number of spines as max_pop_no (= 40)
    Lineage<Cell_Learn> lineage_common;
    const int max_cell_no = std::stoi(parameters["max_cell_no"]);
    for (int i = 0; i != max_cell_no; i++)
    {
        Lineage<Cell_Learn> spine_lin = generate_lineage_learning(setting_rel_path, output_rel_path_4_common, learning_rule, /* enable_common_learning */ true, /*lineage_file_name*/ "spine.dat");
        lineage_common.push(spine_lin);
    }

    //pass for lineage.graphic to calculate lambda of each cell
    auto cell2lambda = [&](const Cell_Learn &cell) {
        return calc_lambda(cell, setting_calc_lambda_dir_rel_path, out_calc_lambda_dir_rel_path);
    };

    lineage_common.graphic(cell2lambda, out_common_learning, out_common_learning_max_min);
}

void compare_common_and_individual_learning_iid()
{
    //set directories (modify if necessary)
    //for generating lineage
    const std::string setting_rel_path = ".//experiments//sim_4_no_growth_learning_artificial";
    const std::string setting_common_rel_path = ".//experiments//sim_4_no_growth_learning_artificial//setting_common";
    const std::string output_rel_path_4_common = ".//experiments//sim_4_no_growth_learning_artificial//res//common";
    const std::string output_rel_path_4_individual = ".//experiments//sim_4_no_growth_learning_artificial//res//learning";

    //for calculating lambda (in order to pass lienage.graphic)
    const std::string setting_calc_lambda_dir_rel_path = ".//experiments//sim_4_no_growth_learning_artificial//calc_lambdas";
    const std::string out_calc_lambda_dir_rel_path = ".//experiments//sim_4_no_growth_learning_artificial//calc_lambdas//res";

    //for drawing graph
    std::ofstream out_ind_learning_lineage(".//experiments//sim_4_no_growth_learning_artificial//res//learning//graph.dat");
    std::ofstream out_common_learning(".//experiments//sim_4_no_growth_learning_artificial//res//common//graph.dat");
    std::ofstream out_whole_lineage(".//experiments//sim_4_no_growth_learning_artificial//res//whole//graph.dat");

    //not necessary, just complete arg. of lienage.grphic
    std::ofstream out_ind_learning_lineage_max_min(".//experiments//sim_4_no_growth_learning_artificial//res//learning//max_min.dat");
    std::ofstream out_common_learning_max_min(".//experiments//sim_4_no_growth_learning_artificial//res//common//max_min.dat");
    std::ofstream out_whole_max_min(".//experiments//sim_4_no_growth_learning_artificial//res//whole//max_min.dat");

    //read paramers to define const variables
    std::ifstream in_other(setting_rel_path + "//other.dat");
    std::map<std::string, std::string> parameters = read_parameters(in_other);

    //set learning rule to online EM
    const int type_no = std::stoi(parameters["type_no"]);
    const double learning_rate = std::stod(parameters["learning_rate"]);
    const double noise_intensity = std::stod(parameters["noise_intensity"]);

    // auto learning_rule = [=](int p_type, int d_type, int no_daughters, std::vector<std::vector<double>> &tran, std::vector<std::vector<double>> &jump_hist, std::vector<double> &rep_hist, std::vector<double> &mem, std::mt19937_64 &mt) {
    //     //update ancestral jump
    //     for (int i = 0; i != type_no; i++)
    //     {
    //         for (int j = 0; j != type_no; j++)
    //         {
    //             jump_hist[i][j] = (1 - learning_rate) * jump_hist[i][j] + learning_rate * ((i == p_type && j == d_type) ? 1.0 : 0.0);
    //         }
    //     }

    //     std::vector<double> pi(type_no, 0.0); //ancestral type_distribution at one-point
    //     for (int i = 0; i != type_no; i++)
    //     {
    //         for (int j = 0; j != type_no; j++)
    //         {
    //             pi[i] += jump_hist[i][j];
    //         }
    //     }

    //     for (int i = 0; i != type_no; i++)
    //     {
    //         for (int j = 0; j != type_no; j++)
    //         {
    //             tran[i][j] = pi[j];
    //         }
    //     }
    // };

    auto learning_rule = [=](int p_type, int d_type, int no_daughters, std::vector<std::vector<double>> &tran, std::vector<std::vector<double>> &jump_hist, std::vector<double> &rep_hist, std::vector<double> &mem, std::mt19937_64 &mt) {
        //iid strategy
        std::vector<double> pi(type_no, 0.0);
        for (int i = 0; i != type_no; i++)
        {
            for (int j = 0; j != type_no; j++)
            {
                pi[i] += jump_hist[i][j];
            }
        }

        //gradient and normalize
        double sum = 0.0;
        for (int i = 0; i != type_no; i++)
        {
            double noise = std::normal_distribution<double>(1.0, noise_intensity)(mt);
            pi[i] = std::max(0.0001,
                             (learning_rate * (i + 1) * noise + (1 - learning_rate)) * pi[i]);
            sum += pi[i];
        }

        //normalize
        for (int i = 0; i != type_no; i++)
        {
            pi[i] /= sum;
        }

        for (int i = 0; i != type_no; i++)
        {
            for (int j = 0; j != type_no; j++)
            {
                jump_hist[i][j] = pi[i];
                tran[i][j] = pi[j];
            }
        }
    };

    //execute simulation and get lineage
    Lineage<Cell_Learn> lineage_ind = generate_lineage_learning(setting_rel_path, output_rel_path_4_individual, learning_rule, /* enable_common_learning */ false);

    //collect the same number of spines as max_pop_no (= 40)
    Lineage<Cell_Learn> lineage_common;
    const int max_cell_no = std::stoi(parameters["max_cell_no"]);
    for (int i = 0; i != max_cell_no; i++)
    {
        Lineage<Cell_Learn> spine_lin = generate_lineage_learning(setting_common_rel_path, output_rel_path_4_common, learning_rule, /* enable_common_learning */ true, /*lineage_file_name*/ "spine.dat");
        lineage_common.push(spine_lin);
    }

    //whole lineage
    Lineage<Cell_Learn> lineage_whole;
    lineage_whole.push(lineage_common);
    lineage_whole.push(lineage_ind);

    //pass for lineage.graphic to calculate lambda of each cell
    auto cell2lambda = [&](const Cell_Learn &cell) {
        return calc_lambda(cell, setting_calc_lambda_dir_rel_path, out_calc_lambda_dir_rel_path);
    };

    //determine the range of plot
    double max_val = std::max(lineage_ind.max(cell2lambda), lineage_common.max(cell2lambda));
    double min_val = std::min(lineage_ind.min(cell2lambda), lineage_common.min(cell2lambda));

    lineage_ind.graphic(cell2lambda, min_val, max_val, out_ind_learning_lineage, out_ind_learning_lineage_max_min);
    lineage_common.graphic(cell2lambda, min_val, max_val, out_common_learning, out_common_learning_max_min);
    lineage_whole.graphic(cell2lambda, out_whole_lineage, out_whole_max_min);

    execute_graphviz(".//experiments//sim_4_no_growth_learning_artificial//res//learning//graph.dat", ".//experiments//sim_4_no_growth_learning_artificial//res//learning//graph.png");
    execute_graphviz(".//experiments//sim_4_no_growth_learning_artificial//res//common//graph.dat", ".//experiments//sim_4_no_growth_learning_artificial//res//common//graph.png");
    execute_graphviz(".//experiments//sim_4_no_growth_learning_artificial//res//whole//graph.dat", ".//experiments//sim_4_no_growth_learning_artificial//res//whole//graph.png");
}

void sim5_random_search_and_growth()
{
    //set directories (modify if necessary)
    //for generating lineage
    const std::string setting_rel_path = ".//experiments//sim_5_random_search_and_growth";
    const std::string setting_common_rel_path = ".//experiments//sim_5_random_search_and_growth//setting_common";
    const std::string output_rel_path_4_random = ".//experiments//sim_5_random_search_and_growth//res//random";
    const std::string output_rel_path_4_adaptive = ".//experiments//sim_5_random_search_and_growth//res//adaptive";
    const std::string output_rel_path_4_learning = ".//experiments//sim_5_random_search_and_growth//res//learning";
    const std::string output_rel_path_4_whole = ".//experiments//sim_5_random_search_and_growth//res//whole";
    const std::string output_rel_path_4_random_plus_adaptive = ".//experiments//sim_5_random_search_and_growth//res//random_adaptive";

    //for calculating lambda (in order to pass lienage.graphic)
    const std::string setting_calc_lambda_dir_rel_path = ".//experiments//sim_5_random_search_and_growth//calc_lambdas";
    const std::string out_calc_lambda_dir_rel_path = ".//experiments//sim_5_random_search_and_growth//calc_lambdas//res";

    //for drawing graph
    std::ofstream out_random_lineage(output_rel_path_4_random + "//graph.dat");
    std::ofstream out_adaptive_learning(output_rel_path_4_adaptive + "//graph.dat");
    std::ofstream out_learning_lineage(output_rel_path_4_learning + "//graph.dat");
    std::ofstream out_whole_lineage(output_rel_path_4_whole + "//graph.dat");
    std::ofstream out_random_plus_adaptive_lineage(output_rel_path_4_random_plus_adaptive + "//graph.dat");

    //not necessary, just complete arg. of lienage.grphic
    std::ofstream out_random_lineage_min_max(output_rel_path_4_random + "//min-max.dat");
    std::ofstream out_adaptive_learning_min_max(output_rel_path_4_adaptive + "//min-max.dat");
    std::ofstream out_learning_lineage_min_max(output_rel_path_4_learning + "//min-max.dat");
    std::ofstream out_whole_lineage_min_max(output_rel_path_4_whole + "//min-max.dat");
    std::ofstream out_random_plus_adaptive_lineage_min_max(output_rel_path_4_random_plus_adaptive + "//min-max.dat");

    //read paramers to define const variables
    std::ifstream in_other(setting_rel_path + "//other.dat");
    std::map<std::string, std::string> parameters = read_parameters(in_other);

    //set learning rule to online EM
    const int type_no = std::stoi(parameters["type_no"]);
    const double learning_rate = std::stod(parameters["learning_rate"]);
    const double noise_intensity = std::stod(parameters["noise_intensity"]);

    auto ancestral_learning = [=](int p_type, int d_type, int no_daughters, std::vector<std::vector<double>> &tran, std::vector<std::vector<double>> &jump_hist, std::vector<double> &rep_hist, std::vector<double> &mem, std::mt19937_64 &mt) {
        //update ancestral jump
        jump_hist[p_type][d_type] += learning_rate;

        //iid strategy
        double sum = 0.0;
        std::vector<double> pi(type_no, 0.0);
        for (int i = 0; i != type_no; i++)
        {
            for (int j = 0; j != type_no; j++)
            {
                pi[i] += jump_hist[j][i];
            }
            sum += pi[i];
        }

        for (int i = 0; i != type_no; i++)
        {
            for (int j = 0; j != type_no; j++)
            {
                jump_hist[i][j] /= sum;
                tran[i][j] = pi[j] / sum;
            }
        }
    };

    auto random_search = [=](int p_type, int d_type, int no_daughters, std::vector<std::vector<double>> &tran, std::vector<std::vector<double>> &jump_hist, std::vector<double> &rep_hist, std::vector<double> &mem, std::mt19937_64 &mt) {
        double sum = 0.0;

        //iid strategy
        std::vector<double> pi(type_no, 0.0);
        for (int i = 0; i != type_no; i++)
        {
            for (int j = 0; j != type_no; j++)
            {
                pi[i] += jump_hist[i][j];
            }
            sum += pi[i];
        }

        //increase a selected component
        std::vector<double> unit(type_no, 1.0);
        std::discrete_distribution<int> dist(unit.begin(), unit.end());
        int selected = dist(mt);
        pi[selected] += noise_intensity;
        sum += noise_intensity;

        //normalize
        for (int i = 0; i != type_no; i++)
        {
            pi[i] /= sum;
        }

        for (int i = 0; i != type_no; i++)
        {
            for (int j = 0; j != type_no; j++)
            {
                jump_hist[i][j] = pi[i];
                tran[i][j] = pi[j];
            }
        }
    };

    auto adaptive_random_search = [=](int p_type, int d_type, int no_daughters, std::vector<std::vector<double>> &tran, std::vector<std::vector<double>> &jump_hist, std::vector<double> &rep_hist, std::vector<double> &mem, std::mt19937_64 &mt) {
        double sum = 0.0;

        //iid strategy
        std::vector<double> pi(type_no, 0.0);
        for (int i = 0; i != type_no; i++)
        {
            for (int j = 0; j != type_no; j++)
            {
                pi[i] += jump_hist[i][j];
            }
            sum += pi[i];
        }

        //increase a selected component
        std::discrete_distribution<int> dist(pi.begin(), pi.end());
        int selected = dist(mt);
        pi[selected] += noise_intensity;
        sum += noise_intensity;

        //normalize
        for (int i = 0; i != type_no; i++)
        {
            pi[i] /= sum;
        }

        for (int i = 0; i != type_no; i++)
        {
            for (int j = 0; j != type_no; j++)
            {
                jump_hist[i][j] = pi[i];
                tran[i][j] = pi[j];
            }
        }
    };

    //execute simulation and get lineage
    Lineage<Cell_Learn> lineage_learning = generate_lineage_learning(setting_rel_path, output_rel_path_4_learning, ancestral_learning, /* enable_common_learning */ false);

    //collect the same number of spines as max_pop_no (= 40)
    // Lineage<Cell_Learn> lineage_random_search;
    // const int max_cell_no = std::stoi(parameters["max_cell_no"]);
    // for (int i = 0; i != max_cell_no; i++)
    // {
    //     Lineage<Cell_Learn> spine_lin = generate_lineage_learning(setting_common_rel_path, output_rel_path_4_random, random_search, /* enable_common_learning */ true, /*lineage_file_name*/ "spine.dat");
    //     lineage_random_search.push(spine_lin);
    // }

    // Lineage<Cell_Learn> lineage_adaptive_random_search;
    // for (int i = 0; i != max_cell_no; i++)
    // {
    //     Lineage<Cell_Learn> spine_lin = generate_lineage_learning(setting_common_rel_path, output_rel_path_4_adaptive, adaptive_random_search, /* enable_common_learning */ true, /*lineage_file_name*/ "spine.dat");
    //     lineage_adaptive_random_search.push(spine_lin);
    // }
    Lineage<Cell_Learn> lineage_random_search = generate_lineage_learning(setting_rel_path, output_rel_path_4_random, random_search, /* enable_common_learning */ false);
    Lineage<Cell_Learn> lineage_adaptive_random_search = generate_lineage_learning(setting_rel_path, output_rel_path_4_adaptive, adaptive_random_search, /* enable_common_learning */ false);

    //whole lineage
    Lineage<Cell_Learn> lineage_whole;
    Lineage<Cell_Learn> lineage_random_plus_adaptive;
    lineage_random_plus_adaptive.push(lineage_random_search);
    lineage_random_plus_adaptive.push(lineage_adaptive_random_search);
    lineage_whole.push(lineage_random_plus_adaptive);
    lineage_whole.push(lineage_learning);

    //pass for lineage.graphic to calculate lambda of each cell
    auto cell2lambda = [&](const Cell_Learn &cell) {
        return calc_lambda(cell, setting_calc_lambda_dir_rel_path, out_calc_lambda_dir_rel_path);
    };

    //determine the range of plot
    const double max_val = std::max({lineage_learning.max(cell2lambda), lineage_random_search.max(cell2lambda), lineage_adaptive_random_search.max(cell2lambda)});
    const double min_val = std::min({lineage_learning.min(cell2lambda), lineage_random_search.min(cell2lambda), lineage_adaptive_random_search.min(cell2lambda)});

    lineage_learning.graphic(cell2lambda, min_val, max_val, out_learning_lineage, out_learning_lineage_min_max);
    lineage_random_search.graphic(cell2lambda, min_val, max_val, out_random_lineage, out_random_lineage_min_max);
    lineage_adaptive_random_search.graphic(cell2lambda, min_val, max_val, out_adaptive_learning, out_adaptive_learning_min_max);
    lineage_whole.graphic(cell2lambda, min_val, max_val, out_whole_lineage, out_whole_lineage_min_max);
    lineage_random_plus_adaptive.graphic(cell2lambda, out_random_plus_adaptive_lineage, out_random_plus_adaptive_lineage_min_max);

    //png
    execute_graphviz(".//experiments//sim_5_random_search_and_growth//res//learning//graph.dat", ".//experiments//sim_5_random_search_and_growth//res//learning//graph.png");
    execute_graphviz(".//experiments//sim_5_random_search_and_growth//res//random//graph.dat", ".//experiments//sim_5_random_search_and_growth//res//random//graph.png");
    execute_graphviz(".//experiments//sim_5_random_search_and_growth//res//adaptive//graph.dat", ".//experiments//sim_5_random_search_and_growth//res//adaptive//graph.png");
    execute_graphviz(".//experiments//sim_5_random_search_and_growth//res//whole//graph.dat", ".//experiments//sim_5_random_search_and_growth//res//whole//graph.png");
    execute_graphviz(output_rel_path_4_random_plus_adaptive + "//graph.dat", output_rel_path_4_random_plus_adaptive + "//graph.png");

    //pdf
    execute_graphviz(".//experiments//sim_5_random_search_and_growth//res//learning//graph.dat", ".//experiments//sim_5_random_search_and_growth//res//learning//graph.pdf", "pdf");
    execute_graphviz(".//experiments//sim_5_random_search_and_growth//res//random//graph.dat", ".//experiments//sim_5_random_search_and_growth//res//random//graph.pdf", "pdf");
    execute_graphviz(".//experiments//sim_5_random_search_and_growth//res//adaptive//graph.dat", ".//experiments//sim_5_random_search_and_growth//res//adaptive//graph.pdf", "pdf");
    execute_graphviz(".//experiments//sim_5_random_search_and_growth//res//whole//graph.dat", ".//experiments//sim_5_random_search_and_growth//res//whole//graph.pdf", "pdf");
    execute_graphviz(output_rel_path_4_random_plus_adaptive + "//graph.dat", output_rel_path_4_random_plus_adaptive + "//graph.pdf", "pdf");
}

double ff_thm_expected_gain(const Cell_Learn &cell, const std::vector<double> &Q_env, const std::vector<std::vector<double>> &mean_replication, double learning_rate)
{
    int type_no = cell.transition.size();
    int env_no = Q_env.size();

    //check the dimension of mean_replication is valid
    assert(mean_replication.size() == env_no);
    for (int i = 0; i != env_no; i++)
    {
        assert(mean_replication[i].size() == type_no);
    }

    //construct iid vector

    //calculate expected gain by definition
    double res = 0.0;
    for (int e = 0; e != env_no; e++)
    {
        for (int ef = 0; ef != env_no; ef++)
        {

            double first_moment1 = 0.0;
            double first_moment2 = 0.0;
            double second_moment = 0.0;
            for (int i = 0; i != type_no; i++)
            {
                first_moment1 += mean_replication[e][i] * cell.transition[0][i];
                first_moment2 += mean_replication[ef][i] * cell.transition[0][i];
                second_moment += mean_replication[e][i] * mean_replication[ef][i] * cell.transition[0][i];
            }
            res += log(second_moment / first_moment1 / first_moment1) * Q_env[e] * Q_env[ef];
        }
    }

    return res * learning_rate;
}

double ff_thm_KL(const Cell_Learn &cell, const std::vector<double> &Q_env, const std::vector<std::vector<double>> &mean_replication, double learning_rate)
{
    int env_no = Q_env.size();
    assert(env_no > 0);
    assert(mean_replication.size() == env_no);
    int type_no = mean_replication[0].size();

    //construct pi_b
    std::vector<std::vector<double>> pi_b;
    for (int y = 0; y != env_no; y++)
    {
        double sum = 0.0;
        std::vector<double> pi_b_y;
        for (int x = 0; x != type_no; x++)
        {
            double pi_b_x_y = mean_replication[y][x] * cell.transition[0][x];
            sum += pi_b_x_y;
            pi_b_y.push_back(pi_b_x_y);
        }
        for (int x = 0; x != type_no; x++)
        {
            pi_b_y[x] /= sum;
        }
        pi_b.push_back(pi_b_y);
    }

    //construct Q_bar
    std::vector<std::vector<double>> Q_bar;
    for (int y1 = 0; y1 != env_no; y1++)
    {
        double sum = 0.0;
        std::vector<double> Q_bar_y1; //Q(y2 | y1)
        for (int y2 = 0; y2 != env_no; y2++)
        {
            double Q_bar_y2_y1 = 0.0;
            for (int x = 0; x != type_no; x++)
            {
                double temp = mean_replication[y1][x] * pi_b[y2][x] * Q_env[y2];
                Q_bar_y2_y1 += temp;
                sum += temp;
            }
            Q_bar_y1.push_back(Q_bar_y2_y1);
        }

        for (int y2 = 0; y2 != env_no; y2++)
        {
            Q_bar_y1[y2] /= sum;
        }
        Q_bar.push_back(Q_bar_y1);
    }

    //calculate KL
    double res = 0.0;
    for (int y1 = 0; y1 != env_no; y1++)
    {
        for (int y2 = 0; y2 != env_no; y2++)
        {
            res += Q_env[y1] * Q_env[y2] * log(Q_env[y2] / Q_bar[y2][y1]);
        }
    }

    return res * learning_rate;
}

void check_ff_thm_from_lineage(Lineage<Cell_Learn> &lineage, std::function<double(Cell_Learn)> calc_lambda, std::function<double(Cell_Learn)> calc_expected_gain, std::function<double(Cell_Learn)> calc_KL, std::ofstream &out, std::ofstream &out_detail, const int pickup_no = std::numeric_limits<int>::max())
//cell.memory[0] = 0 if no learning occurs and > 0 otherwise
{
    //store fitness gain of the cells in lienage
    std::vector<double> expected_gains;
    std::vector<double> actual_gains;

    const size_t end_time = lineage.endtime();
    int counter = 1;
    for (int t = 0; t != end_time; t++)
    {
        for (int i = 0; i != lineage.pop_size(t); i++)
        {
            Cell_Learn cell = lineage.access(t, i);
            assert(cell.memory.size() > 0);
            if (cell.memory[0] > 0.01) //no learning, 0.01 is machine epsilon
                continue;

            //skip root cell
            if (cell.id().find('S') == std::string::npos)
                continue;

            Cell_Learn parent = lineage.parent(cell);

            double lambda_cell = calc_lambda(cell);
            double lambda_parent = calc_lambda(parent);
            double actual_gain = lambda_cell - lambda_parent;
            actual_gains.push_back(lambda_cell - lambda_parent);

            double expected_gain = calc_expected_gain(parent);
            expected_gains.push_back(expected_gain);

            double KL = calc_KL(parent);

            //output
            out << actual_gain << " " << expected_gain << " " << KL << std::endl;
            out_detail << counter << std::endl;
            parent.record(&out_detail);
            cell.record(&out_detail);
            counter++;
            if (counter > pickup_no)
                return;
        }
    }

    //drawing, to be implemented if necesarry
}

void check_ff_thm_from_path(const std::string &setting_rel_path, const std::string &output_rel_path, const std::string &setting_for_calc_lambda_rel_path, const std::string &output_for_calc_lambda_rel_path) //only for iid env. TODO? impletemnt Markov ver.
{

    //read paramers to define const variables
    std::ifstream in_other(setting_rel_path + "//other.dat");
    std::map<std::string, std::string> parameters = read_parameters(in_other);

    //set learning rule to online EM
    const int type_no = std::stoi(parameters["type_no"]);
    const double learning_rate = std::stod(parameters["learning_rate"]);
    const int time_estimate_retro = std::stoi(parameters["time_estimate_retro"]);

    //learning rule, to be updated
    auto ancestral_learning = [=](int p_type, int d_type, int no_daughters, std::vector<std::vector<double>> &tran, std::vector<std::vector<double>> &jump_hist, std::vector<double> &rep_hist, std::vector<double> &mem, std::mt19937_64 &mt) {
        //wait time_estimate_retro generations w/o learning and collect empirical distribution
        if (mem[0] + 1.0 < time_estimate_retro - 0.1) //to avoid numerical error, I inserted "-0.1"
        {
            mem[0] = floor(mem[0] + 1.1); //to avoid numerical error

            //update ancestral jump
            jump_hist[p_type][d_type] += 1.0 / time_estimate_retro;
            return;
        }

        //else learing
        mem[0] = 0.0;

        //iid strategy
        double sum_pi = 0.0;
        double sum_original = 0.0;
        std::vector<double> pi(type_no, 0.0);
        std::vector<double> original_pi(type_no, 0.0);
        for (int i = 0; i != type_no; i++)
        {
            for (int j = 0; j != type_no; j++)
            {
                pi[i] += jump_hist[j][i];
                sum_pi += jump_hist[j][i];
                original_pi[j] = tran[i][j];
                sum_original += tran[i][j];
            }
        }

        for (int i = 0; i != type_no; i++)
        {
            for (int j = 0; j != type_no; j++)
            {
                jump_hist[i][j] = 0.0;                                                                               //clear history for the next estimation of the empirical distribution
                tran[i][j] = learning_rate * pi[j] / sum_pi + (1.0 - learning_rate) * original_pi[j] / sum_original; //update
            }
        }
    };

    //execute simulation and get lineage
    Lineage<Cell_Learn> lineage = generate_lineage_learning(setting_rel_path, output_rel_path, ancestral_learning, /* enable_common_learning */ false);

    //pass for check_ff_thm_from_lineage to calculate ***actual*** fitness gain of each cell
    auto cell2lambda = [&](const Cell_Learn &cell) {
        return calc_lambda(cell, setting_for_calc_lambda_rel_path, output_for_calc_lambda_rel_path);
    };

    //pass for check_ff_thm_from_lineage to calculate ***expected*** fitness gain of each cell
    //first read parameters
    std::ifstream in_env(setting_rel_path + "//env.dat");
    Markov_Environments env = read_env(in_env);
    std::vector<double> Q_env;
    for (int y = 0; y != env.cardinality(); y++)
    {
        Q_env.push_back(env.get_transition(0, y));
    }

    //replication
    std::vector<std::vector<std::vector<double>>> replication;
    std::ifstream in_repl(setting_rel_path + "//replication.dat");
    read3DTensor<double>(replication, in_repl);
    std::vector<std::vector<double>> mean_replication;
    for (int y = 0; y != Q_env.size(); y++)
    {
        std::vector<double> mean_vec;
        for (int x = 0; x != type_no; x++)
        {
            double mean = 0.0;
            for (int i = 0; i != replication[y][x].size(); i++)
            {
                mean += replication[y][x][i] * i;
            }
            mean_vec.push_back(mean);
        }
        mean_replication.push_back(mean_vec);
    }

    auto cell2expected = [&](const Cell_Learn &cell) {
        return ff_thm_expected_gain(cell, Q_env, mean_replication, learning_rate);
    }; //check ff-thm

    auto cell2KL = [&](const Cell_Learn &cell) {
        return ff_thm_KL(cell, Q_env, mean_replication, learning_rate);
    }; //check ff-thm

    std::ofstream out_ff_thm(output_rel_path + "//ff_thm.dat");
    std::ofstream out_ff_thm_detail(output_rel_path + "//ff_thm_detail.dat");
    check_ff_thm_from_lineage(lineage, cell2lambda, cell2expected, cell2KL, out_ff_thm, out_ff_thm_detail);
}

void generate_random_initial_cell_learn(const std::string &setting_dir_rel_path, const int type_no, const int cell_no, const int max_cell_no) //memory structure is fixed
{
    std::ofstream initial_cell(setting_dir_rel_path + "//initial_cells.dat");

    initial_cell << type_no << " " << cell_no << " " << max_cell_no << std::endl
                 << std::endl;

    std::random_device rnd;
    std::mt19937_64 mt(rnd());
    std::uniform_int_distribution<> rand_type(0, type_no - 1);
    std::uniform_real_distribution<double> rand_transition(0, 1);

    for (int i = 0; i != cell_no; i++)
    {
        //type
        initial_cell << rand_type(mt) << std::endl
                     << std::endl;

        //jump history
        for (int i = 0; i != type_no; i++)
        {
            for (int j = 0; j != type_no; j++)
            {
                initial_cell << "0 ";
            }
            initial_cell << std::endl;
        }
        initial_cell << std::endl;

        //construct transition
        std::vector<double> pi(type_no);
        double sum = 0.0;
        for (int j = 0; j != type_no; j++)
        {
            pi[j] = rand_transition(mt);
            sum += pi[j];
        }
        //normalization and output
        for (int j = 0; j != type_no; j++)
        {
            for (int k = 0; k != type_no; k++)
            {
                initial_cell << pi[k] / sum << " ";
            }
            initial_cell << std::endl;
        }
        initial_cell << std::endl;

        //replication history
        initial_cell << 0 << std::endl
                     << std::endl;

        //memory
        initial_cell << "1 0\n"
                     << std::endl;
    }
}

void check_ff_thm_from_path_random_transition(const std::string &setting_rel_path, const std::string &output_rel_path, const std::string &setting_for_calc_lambda_rel_path, const std::string &output_for_calc_lambda_rel_path) //only for iid env. TODO? impletemnt Markov ver.
{

    //read paramers to define const variables
    std::ifstream in_other(setting_rel_path + "//other.dat");
    std::map<std::string, std::string> parameters = read_parameters(in_other);

    //set learning rule to online EM
    const int type_no = std::stoi(parameters["type_no"]);
    const double learning_rate = std::stod(parameters["learning_rate"]);
    const int time_estimate_retro = std::stoi(parameters["time_estimate_retro"]);
    const int max_cell_no = std::stoi(parameters["max_cell_no"]);
    const int init_cell_no = std::stoi(parameters["init_cell_no"]);
    const int pickup_no = std::stoi(parameters["pickup_no"]);

    //learning rule, to be updated
    auto ancestral_learning = [=](int p_type, int d_type, int no_daughters, std::vector<std::vector<double>> &tran, std::vector<std::vector<double>> &jump_hist, std::vector<double> &rep_hist, std::vector<double> &mem, std::mt19937_64 &mt) {
        //wait time_estimate_retro generations w/o learning and collect empirical distribution
        if (mem[0] + 1.0 < time_estimate_retro - 0.1) //to avoid numerical error, I inserted "-0.1"
        {
            mem[0] = floor(mem[0] + 1.1); //to avoid numerical error

            //update ancestral jump
            jump_hist[p_type][d_type] += 1.0 / time_estimate_retro;
            return;
        }

        //else learing
        mem[0] = 0.0;

        //iid strategy
        double sum_pi = 0.0;
        double sum_original = 0.0;
        std::vector<double> pi(type_no, 0.0);
        std::vector<double> original_pi(type_no, 0.0);
        for (int i = 0; i != type_no; i++)
        {
            for (int j = 0; j != type_no; j++)
            {
                pi[i] += jump_hist[j][i];
                sum_pi += jump_hist[j][i];
                original_pi[j] = tran[i][j];
                sum_original += tran[i][j];
            }
        }

        for (int i = 0; i != type_no; i++)
        {
            for (int j = 0; j != type_no; j++)
            {
                jump_hist[i][j] = 0.0;
                //clear history for the next estimation of the empirical distribution
                tran[i][j] = learning_rate * pi[j] / sum_pi + (1.0 - learning_rate) * original_pi[j] / sum_original; //update
            }
        }
    };

    //generate random transition
    generate_random_initial_cell_learn(setting_rel_path, type_no, init_cell_no, max_cell_no);

    //execute simulation and get lineage
    Lineage<Cell_Learn>
        lineage = generate_lineage_learning(setting_rel_path, output_rel_path, ancestral_learning, /* enable_common_learning */ false);

    //pass for check_ff_thm_from_lineage to calculate ***actual*** fitness gain of each cell
    auto cell2lambda = [&](const Cell_Learn &cell) {
        return calc_lambda(cell, setting_for_calc_lambda_rel_path, output_for_calc_lambda_rel_path);
    };

    //pass for check_ff_thm_from_lineage to calculate ***expected*** fitness gain of each cell
    //first read parameters
    std::ifstream in_env(setting_rel_path + "//env.dat");
    Markov_Environments env = read_env(in_env);
    std::vector<double> Q_env;
    for (int y = 0; y != env.cardinality(); y++)
    {
        Q_env.push_back(env.get_transition(0, y));
    }

    //replication
    std::vector<std::vector<std::vector<double>>> replication;
    std::ifstream in_repl(setting_rel_path + "//replication.dat");
    read3DTensor<double>(replication, in_repl);
    std::vector<std::vector<double>> mean_replication;
    for (int y = 0; y != Q_env.size(); y++)
    {
        std::vector<double> mean_vec;
        for (int x = 0; x != type_no; x++)
        {
            double mean = 0.0;
            for (int i = 0; i != replication[y][x].size(); i++)
            {
                mean += replication[y][x][i] * i;
            }
            mean_vec.push_back(mean);
        }
        mean_replication.push_back(mean_vec);
    }

    auto cell2expected = [&](const Cell_Learn &cell) {
        return ff_thm_expected_gain(cell, Q_env, mean_replication, learning_rate);
    };
    auto cell2KL = [&](const Cell_Learn &cell) {
        return ff_thm_KL(cell, Q_env, mean_replication, learning_rate);
    };

    //check ff-thm
    std::ofstream out_ff_thm(output_rel_path + "//ff_thm.dat", std::ios::app);
    std::ofstream out_ff_thm_detail(output_rel_path + "//ff_thm_detail.dat", std::ios::app);
    check_ff_thm_from_lineage(lineage, cell2lambda, cell2expected, cell2KL, out_ff_thm, out_ff_thm_detail, pickup_no);
}

bool starts_with(const std::string &s, const std::string &prefix)
{
    auto size = prefix.size();
    if (s.size() < size)
        return false;
    return std::equal(std::begin(prefix), std::end(prefix), std::begin(s));
}

std::vector<Cell_Learn> read_selected_generation(const std::string &file_rel_path, const int target_gen, const int type_no)
{
    std::ifstream in(file_rel_path);

    //for loading
    std::string id;
    int type;
    std::vector<std::vector<double>> ancestral_jump;
    std::vector<std::vector<double>> transition;
    std::vector<double> replication_history;
    std::vector<double> mem;

    std::vector<Cell_Learn> res;

    while (in >> id)
    {
        //load a cell from lineage
        in >> type;
        readMat(type_no, type_no, ancestral_jump, in);
        readMat(type_no, type_no, transition, in);
        //readVec(type_no, replication_history, in);
        int memory_no;
        in >> memory_no;
        readVec(memory_no, mem, in);

        //check generation
        const int gen = std::count(id.begin(), id.end(), 'S');
        if (gen == target_gen)
        {

            Cell_Learn cell(type, id, ancestral_jump, transition, replication_history, mem);
            res.push_back(cell);
        }
    }

    return res;
}
void check_ff_from_file_one_path(const int type_no, const int end_time, const std::string &file_rel_path, std::ofstream &out, std::function<double(Cell_Learn)> calc_lambda, std::function<double(Cell_Learn)> calc_expected_gain, std::function<double(Cell_Learn)> calc_KL)
//calculate gains if the cell is on the selected path
{

    //path information
    std::string path;

    for (int gen = end_time; gen > 0; gen--)
    {
        //load a cell from lineage
        std::vector<Cell_Learn> current_pop = read_selected_generation(file_rel_path, gen, type_no);
        std::vector<Cell_Learn> prev_pop = read_selected_generation(file_rel_path, gen - 1, type_no);

        //skip the generation w/o learning
        assert(current_pop.size() > 0);
        assert(current_pop[0].memory.size() > 0);
        if (current_pop[0].memory[0] > 0.01) //no learning, 0.01 is machine epsilon, also skip root cells
            continue;

        Cell_Learn cell;

        //serach a cell on the selected retrospective path
        for (auto cell_cur : current_pop)
        {
            if (path.empty() || starts_with(path, cell_cur.id())) //at the end of lineage, we can chose arbitrary sell (path.empty())
            {
                cell = cell_cur;
                break;
            }
        }
        //update path infomation
        path = cell.id();

        //search parent cell
        Cell_Learn parent;
        for (auto cell_prev : prev_pop)
        {
            if (cell_prev.id() == parent_id(cell.id()))
            {
                parent = cell_prev;
                break;
            }
        }

        //actual gain
        const double lambda_cell = calc_lambda(cell);
        const double lambda_parent = calc_lambda(parent);
        const double actual_gain = lambda_cell - lambda_parent;

        //expected gain (variance and KL)
        const double expected_gain = calc_expected_gain(parent);
        double KL = calc_KL(parent);

        //output
        out << actual_gain << " " << expected_gain << " " << KL << std::endl;
    }
}

void check_ff_thm_one_path(const std::string &setting_rel_path, const std::string &output_rel_path, const std::string &setting_for_calc_lambda_rel_path, const std::string &output_for_calc_lambda_rel_path) //only for iid env. TODO? impletemnt Markov ver.
{

    //read paramers to define const variables
    std::ifstream in_other(setting_rel_path + "//other.dat");
    std::map<std::string, std::string> parameters = read_parameters(in_other);

    //set learning rule to online EM
    const int type_no = std::stoi(parameters["type_no"]);
    const double learning_rate = std::stod(parameters["learning_rate"]);
    const int time_estimate_retro = std::stoi(parameters["time_estimate_retro"]);
    const int end_time = std::stoi(parameters["end_time"]);

    //learning rule, to be updated
    auto ancestral_learning = [=](int p_type, int d_type, int no_daughters, std::vector<std::vector<double>> &tran, std::vector<std::vector<double>> &jump_hist, std::vector<double> &rep_hist, std::vector<double> &mem, std::mt19937_64 &mt) {
        //wait time_estimate_retro generations w/o learning and collect empirical distribution
        if (mem[0] + 1.0 < time_estimate_retro - 0.1) //to avoid numerical error, I inserted "-0.1"
        {
            mem[0] = floor(mem[0] + 1.1); //to avoid numerical error

            //update ancestral jump
            jump_hist[p_type][d_type] += 1.0 / time_estimate_retro;
            return;
        }

        //else learing
        mem[0] = 0.0;

        //iid strategy
        double sum_pi = 0.0;
        double sum_original = 0.0;
        std::vector<double> pi(type_no, 0.0);
        std::vector<double> original_pi(type_no, 0.0);
        for (int i = 0; i != type_no; i++)
        {
            for (int j = 0; j != type_no; j++)
            {
                pi[i] += jump_hist[j][i];
                sum_pi += jump_hist[j][i];
                original_pi[j] = tran[i][j];
                sum_original += tran[i][j];
            }
        }

        for (int i = 0; i != type_no; i++)
        {
            for (int j = 0; j != type_no; j++)
            {
                jump_hist[i][j] = 0.0;                                                                               //clear history for the next estimation of the empirical distribution
                tran[i][j] = learning_rate * pi[j] / sum_pi + (1.0 - learning_rate) * original_pi[j] / sum_original; //update
            }
        }
    };

    //execute simulation and get lineage
    sim_learning(setting_rel_path, output_rel_path, ancestral_learning, /* enable_common_learning */ false);

    //pass for check_ff_thm_from_lineage to calculate ***actual*** fitness gain of each cell
    auto cell2lambda = [&](const Cell_Learn &cell) {
        return calc_lambda(cell, setting_for_calc_lambda_rel_path, output_for_calc_lambda_rel_path);
    };

    //pass for check_ff_thm_from_lineage to calculate ***expected*** fitness gain of each cell
    //first read parameters
    std::ifstream in_env(setting_rel_path + "//env.dat");
    Markov_Environments env = read_env(in_env);
    std::vector<double> Q_env;
    for (int y = 0; y != env.cardinality(); y++)
    {
        Q_env.push_back(env.get_transition(0, y));
    }

    //replication
    std::vector<std::vector<std::vector<double>>> replication;
    std::ifstream in_repl(setting_rel_path + "//replication.dat");
    read3DTensor<double>(replication, in_repl);
    std::vector<std::vector<double>> mean_replication;
    for (int y = 0; y != Q_env.size(); y++)
    {
        std::vector<double> mean_vec;
        for (int x = 0; x != type_no; x++)
        {
            double mean = 0.0;
            for (int i = 0; i != replication[y][x].size(); i++)
            {
                mean += replication[y][x][i] * i;
            }
            mean_vec.push_back(mean);
        }
        mean_replication.push_back(mean_vec);
    }

    auto cell2expected = [&](const Cell_Learn &cell) {
        return ff_thm_expected_gain(cell, Q_env, mean_replication, learning_rate);
    }; //check ff-thm

    auto cell2KL = [&](const Cell_Learn &cell) {
        return ff_thm_KL(cell, Q_env, mean_replication, learning_rate);
    }; //check ff-thm

    std::string pop_file_name = output_rel_path + "//pop.dat";
    std::ofstream out_ff_thm(output_rel_path + "//ff_thm.dat");
    check_ff_from_file_one_path(type_no, end_time, pop_file_name, out_ff_thm, cell2lambda, cell2expected, cell2KL);
}

void sim_6_check_ff_thm()
{
    //check const-env
    // check_ff_thm_from_path(
    //     ".//experiments//sim_6_ffthm//const_env",
    //     ".//experiments//sim_6_ffthm//const_env//res",
    //     ".//experiments//sim_6_ffthm//const_env//calc_lambdas",
    //     ".//experiments//sim_6_ffthm//const_env//calc_lambdas//res");

    //check const-end-random-gen
    // check_ff_thm_from_path_random_transition(
    //     ".//experiments//sim_6_ffthm//const_env_random_transition",
    //     ".//experiments//sim_6_ffthm//const_env_random_transition//res",
    //     ".//experiments//sim_6_ffthm//const_env_random_transition//calc_lambdas",
    //     ".//experiments//sim_6_ffthm//const_env_random_transition//calc_lambdas//res");

    // check_ff_thm_from_path_random_transition(
    //     ".//experiments//sim_6_ffthm//non_const_env",
    //     ".//experiments//sim_6_ffthm//non_const_env//res",
    //     ".//experiments//sim_6_ffthm//non_const_env//calc_lambdas",
    //     ".//experiments//sim_6_ffthm//non_const_env//calc_lambdas//res");

    // check_ff_thm_from_path_random_transition(
    //     ".//experiments//sim_6_ffthm//bet_hedge_plus_concentration",
    //     ".//experiments//sim_6_ffthm//bet_hedge_plus_concentration//res",
    //     ".//experiments//sim_6_ffthm//bet_hedge_plus_concentration//calc_lambdas",
    //     ".//experiments//sim_6_ffthm//bet_hedge_plus_concentration//calc_lambdas//res");

    check_ff_thm_one_path(
        ".//experiments//sim_6_ffthm//const_env_path",
        ".//experiments//sim_6_ffthm//const_env_path//res",
        ".//experiments//sim_6_ffthm//const_env_path//calc_lambdas",
        ".//experiments//sim_6_ffthm//const_env_path//calc_lambdas//res");
}
