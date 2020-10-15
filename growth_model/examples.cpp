#include "headers/simulators.h"
#include <vector>
#include "headers/simulator_utility.h"
#include "headers/analyze.h"

void test_env(){
    MBPRE w;

    w.set_end_time(30);
    
    std::vector<double> init = {1.0, 1.0};
    std::vector<std::vector<double>> tran = {{1.0, 0}, {0.5, 0.5}};
    Markov_Environments env(tran, 2, init);

    std::ofstream out_env(".//res//environments.dat");
    w.set_env_record(&out_env);
    w.set_environments(&env);

    w.excecute();
}

void test_env_cells(){
        MBPRE w;

    w.set_end_time(20);
    
    //env
    std::vector<double> init = {1.0, 1.0};
    std::vector<std::vector<double>> env_tran = {{0.5, 0.5}, {0.5, 0.5}};
    Markov_Environments env(env_tran, 2, init);

    std::ofstream out_env(".//res//environments.dat");
    w.set_env_record(&out_env);
    w.set_environments(&env);


    //cells
    std::vector<std::vector<double>> type_tran = {{0.8, 0.2}, {0.2,0.8}};
    std::vector<Cell> init_pop = {Cell(0, "0")};//, Cell(0, "1")};

    Cells cells(2, type_tran, init_pop);
    std::ofstream out_pop(".//res//population.dat");
    std::ofstream out_pop_full(".//res//population_full.dat");
    cells.set_maximum_population_size(5);

    w.set_population(&cells);
    w.set_pop_record(&out_pop);
    w.set_pop_full_record(&out_pop_full);


    //replication
    std::vector<std::vector<std::vector<double>>> replication = {{{0.0, 0.0, 1.0}, {0.0, 0.0, 1.0}}, {{0.0, 0.0, 1.0},{0.0, 0.0, 1.0}}};
    w.set_offspring_distributions(replication);

    w.excecute();
}

void file_read_test(){
    //initialize world setting
    std::ifstream in_other(".//experiments//sim_1//world_other.dat");
    MBPRE w = read_mbpre(in_other);


    std::ifstream in_env(".//experiments//sim_1//env.dat");
    Markov_Environments env = read_env(in_env);
    std::ifstream in_cells(".//experiments//sim_1//initial_cells.dat");
    std::ifstream in_cell_tran(".//experiments//sim_1//cell_type_tran.dat");
    Cells cells = read_cells(in_cells, in_cell_tran);
    cells.set_maximum_population_size(10);


    //replication
    std::vector<std::vector<std::vector<double>>> replication;
    std::ifstream in_repl(".//experiments//sim_1//replication.dat");
    read3DTensor<double>(replication, in_repl);
    w.set_offspring_distributions(replication);

    //record
    std::ofstream out_env(".//experiments//sim_1//res//env.dat");
    w.set_env_record(&out_env);
    w.set_environments(&env);
    w.set_population(&cells);
    std::ofstream out_pop(".//experiments//sim_1//res//pop.dat");
    std::ofstream out_pop_full(".//experiments//sim_1//res//pop_full.dat");
    w.set_pop_record(&out_pop);
    w.set_pop_full_record(&out_pop_full);

    w.excecute();
}

void test_analyze(){
    std::ifstream in_env(".//experiments//sim_1//res//env.dat");
    std::vector<int> environments;
    readVec(5000, environments, in_env);

    std::ifstream in_pop(".//experiments//sim_1//res//pop.dat");
    std::ifstream in_pop_full(".//experiments//sim_1//res//pop_full.dat");
    Lineage<Cell> lineage = read_lineage(in_pop);
    Lineage<Cell> lienage_full = read_lineage(in_pop_full);

    std::ofstream out_lambda(".//experiments//sim_1//res//analyze//lambda.dat");
    out_lambda << lienage_full.lambda(10) << std::endl;

    std::ofstream out_retro(".//experiments//sim_1//res//analyze//retro.dat");
    std::function<double(Cell,Cell,int, int, int)> func = [](Cell p, Cell c,  int g, int p_type, int c_type){ //retrospective type transition
        if(p.type() == p_type && c.type() == c_type){
            return 1.0;
        }
        else
        {
            return 0.0;
        }
        
    };

    std::function<double(Cell,Cell,int, int, int, const std::vector<int>&)> func2 = [](Cell p,Cell c,int g, int c_type, int env, const std::vector<int>& envs){
        if(c.type() == c_type && envs[g] == env){
            return 1.0;
        }
        else{
            return 0.0;
        }
    };

    std::vector<std::function<double(Cell, Cell, int)>> funcs;

    for(int i = 0; i != 2; i++){
        for(int j = 0; j != 2; j++){
            using namespace std::placeholders;
            std::function<double(Cell,Cell,int)> temp = std::bind(func, _1, _2,_3, i, j);
            funcs.push_back(temp);
        }
    }

    for(int i = 0; i != 2; i++){
        for(int j = 0; j != 2; j++){
            using namespace std::placeholders;
            std::function<double(Cell,Cell,int)> temp = std::bind(func2, _1,_2,_3, j, i, environments);   
            funcs.push_back(temp);
        }
    }

    std::vector<double> res = lineage.backward_mean(funcs);
    for(auto d : res){
        out_retro << d << std::endl;
    }
}

void lambda_curve(){
    std::ofstream out_lambda_curve(".//experiments//sim_1//res//analyze//lambda_curve.dat");

    int mean_no = 200;
    for(int t = - 30; t != 31; t++){
        double mean = 0.0;
        for(int i = 0; i != mean_no; i++){
            //initialize world setting
            std::ifstream in_other(".//experiments//sim_1//world_other.dat");
            MBPRE w = read_mbpre(in_other);


            std::ifstream in_env(".//experiments//sim_1//env.dat");
            Markov_Environments env = read_env(in_env);
            std::ifstream in_cells(".//experiments//sim_1//initial_cells.dat");
            std::ifstream in_cell_tran(".//experiments//sim_1//cell_type_tran.dat");
            Cells cells = read_cells(in_cells, in_cell_tran);
            cells.set_maximum_population_size(100);
/*             std::vector<std::vector<double>> cell_tran = {{0.5,0.5}, {0.5,0.5}};
            cell_tran[0][0] = 0.7 + t * 0.01;
            cell_tran[0][1] = 0.3 - t * 0.01;
            cell_tran[1][0] = 0.3 - t * 0.01;
            cell_tran[1][1] = 0.7 + t * 0.01;
            cells.set_type_transition(cell_tran); */


            //replication
            std::vector<std::vector<std::vector<double>>> replication;
            std::ifstream in_repl(".//experiments//sim_1//replication.dat");
            read3DTensor<double>(replication, in_repl);

            replication[0][0][4] = 0.0;
            replication[0][0][3] = .5 - t / 60.0;
            replication[0][0][9] = .5 + t / 60.0;
            //replication[1][1][4] = 0.0;
            //replication[1][1][1] = .5 - t / 60.0;
            //replication[1][1][7] = .5 + t / 60.0;

            w.set_offspring_distributions(replication);
            w.set_end_time(100);


            //record
            std::ofstream out_env(".//experiments//sim_1//res//env.dat");
            w.set_env_record(&out_env);
            w.set_environments(&env);
            w.set_population(&cells);
            std::ofstream out_pop(".//experiments//sim_1//res//pop.dat");
            std::ofstream out_pop_full(".//experiments//sim_1//res//pop_full.dat");
            w.set_pop_record(&out_pop);
            w.set_pop_full_record(&out_pop_full);

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

    auto learning_rule
    = [=](int p_type, int d_type, int no_daughters, std::vector<std::vector<double>>& tran, std::vector<std::vector<double>>& jump_hist, std::vector<double>& rep_hist, std::vector<double>& mem, std::mt19937_64& mt){

        //update ancestral jump
        for(int i = 0; i != type_no; i++){
            for(int j = 0; j != type_no; j++){
                jump_hist[i][j] = (1 - learning_rate) * jump_hist[i][j] + learning_rate * ((i == p_type && j == d_type)? 1.0 : 0.0);
            }
        }

        //update transition using ancestral jump
        tran = jump_hist;
    };

    cells.set_learning_rule(learning_rule);

    //replication
    std::vector<std::vector<std::vector<double>>> replication;
    std::ifstream in_repl(".//experiments//sim_2//replication.dat");
    read3DTensor<double>(replication, in_repl);
    w.set_offspring_distributions(replication);

    //record
    std::ofstream out_env(".//experiments//sim_2//res//env.dat");
    w.set_env_record(&out_env);
    w.set_environments(&env);
    w.set_population(&cells);
    std::ofstream out_pop(".//experiments//sim_2//res//pop.dat");
    std::ofstream out_pop_full(".//experiments//sim_2//res//pop_full.dat");
    w.set_pop_record(&out_pop);
    w.set_pop_full_record(&out_pop_full);

    w.excecute();
}