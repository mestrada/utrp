#ifndef SOLUTION_H
#define SOLUTION_H

#include <vector>
using namespace std;

class solution
{
    public:
        
        solution();
        solution(int sizeX, int sizeY);
        ~solution();
        void generate_solution(int **m_tt);
        void print_fact_routes(void);
        void print(void);
        bool is_in(int);
        std::vector<int> get_current_sol(void);
        int evaluate_time(int **tt);
        int evaluate_demand(int **td);
        int evaluate_cost(int**);
        bool change_sol(void);
        void reset(void);
        int route_lenght(int**, int);
        int get_n_routes(void);

    private:
        int **sol_m; 
        int n_nodes;
        int n_routes;
        vector<int> sol;

        //vector<coord> camino;
        //vector<coord> solution;
};
#endif
