#include "headers/simulators.h"
#include <vector>
#include "headers/simulator_utility.h"
#include "headers/analyze.h"
#include "headers/linterp.h"
#include <iterator>
#include <algorithm>

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
    Lineage<Cell_Learn> lineage = read_learning_lineage(type_no, mem_no, in_pop);
    Lineage<Cell_Learn> lineage_full = read_learning_lineage(type_no, mem_no, in_pop_full);

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

void sim_learning(std::string setting_dir_rel_path, std::string output_dir_rel_path, std::function<void(int p_type, int d_type, int no_daughters, std::vector<std::vector<double>> &tran, std::vector<std::vector<double>> &jump_hist, std::vector<double> &rep_hist, std::vector<double> &mem, std::mt19937_64 &mt)> learning_rule, bool enable_common_learning = false)
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
    if (enable_common_learning)
    {
        cells = new_cells_learn_common_from_read(in_cells);
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

    w.excecute();

    delete cells;
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
    sim_learning(setting_rel_path, output_rel_path_4_individual, learning_rule, false);
}

void compare_common_and_individual_learning()
{
    //set directories
    const std::string setting_rel_path = ".//experiments//sim_2_no_growth_comp";
    const std::string output_rel_path_4_common = ".//experiments//sim_2_no_growth_comp//res//common";
    const std::string output_rel_path_4_individual = ".//experiments//sim_2_no_growth_comp//res//learning";

    //read paramers and define const variables
    std::ifstream in_other(setting_rel_path + "//other.dat");
    std::map<std::string, std::string> parameters = read_parameters(in_other);

    //online EM
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

    //execute simulation
    sim_learning(setting_rel_path, output_rel_path_4_common, learning_rule, true);
    sim_learning(setting_rel_path, output_rel_path_4_individual, learning_rule, false);

    //draw graphs
    //construct lineage
    const int endtime = std::stoi(parameters["end_time"]);
    const int mem_no = std::stoi(parameters["mem_no"]);
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
    Lineage<Cell_Learn> lineage = read_learning_lineage(type_no, mem_no, in_pop);
    Lineage<Cell_Learn> lienage_orginal = read_learning_lineage(type_no, mem_no, in_pop);

    lineage.push(lineage);

    lineage.push(lineage);
}
