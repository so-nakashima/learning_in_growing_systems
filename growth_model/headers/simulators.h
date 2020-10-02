#include<vector>
#include <random>

class World_base{
public:
    virtual void execute();
};

class World_MBPRE : public World_base
{
public:
    World_MBPRE(/* args */);
    ~World_MBPRE();

    
    void excecute();

private:
    
};

class Environments
{
private:
    //initial data
    int cardinality = 1;
    std::vector<std::vector<int>> transition; //transition[i][j] = prob. of the transition form i to j.
    int initial_state = 0;

    //states
    int current_state = 0;

    //random number geeerators
    std::mt19937_64 mt;
    std::uniform_int_distribution<> uniform;

public:
    Environments(/* args */);
    ~Environments();
    void set_initial(int init);
    void set_cardinality(int n);
    void set_transition(const std::vector<std::vector<int>>& tran_mat);
    int next_state();
};

Environments::Environments(/* args */)
{
    std::random_device rnd;
    mt.seed(rnd());

    uniform = std::uniform_int_distribution(0, cardinality);
}

Environments::~Environments()
{
    
}

