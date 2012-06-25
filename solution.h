#ifndef SOLUTION_H
#define SOLUTION_H

#include <vector>
using namespace std;

class solution
{
    public:
        
        solution();
        solution(int, int);
        ~solution();
        void generate_solution(int **);
        void find_top_demand_nodes(int, int **);
        void reorder();
        void swap(int, int);
        bool is_distinct(int);
        void print_fact_routes(void);
        void print(void);
        bool is_in(int);
        std::vector<int> get_current_sol(void);
        int evaluate_time(int **);
        int evaluate_demand(int **);
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
        int *top_demand_nodes;
        int *top_demand;
        int top_tam;
        int last_index = 0;

        //vector<coord> camino;
        //vector<coord> solution;
};
#endif
