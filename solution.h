#ifndef SOLUTION_H
#define SOLUTION_H

#include <vector>
using namespace std;

typedef struct pair
{
    int start;
    int end;

}antigen;


typedef struct mut_pair{
        int ind;
        int ite;
}mut;


class solution
{
    public:
        
        solution();
        solution(int, int, int, int);
        ~solution();
        void generate_solution(int **);
        void generate_solution(int **, double);
        void generate_antigens(void);
        void generate_antigens(double);
        void print_antigens(void);
        int node_count(int);
        void find_top_demand_nodes(int, int **);
        void reorder();
        void swap(int, int);
        bool is_distinct(int);
        void print_fact_routes(void);
        void print(void);
        bool is_in(int);
        bool is_in(int, int);
        std::vector<int> get_current_sol(void);
        int evaluate_time(int **);
        long evaluate_time(int **, int **);
        long evaluate_demand(int **, int **);
        long evaluate_cost(int**);
        bool change_sol(void);
        void reset(void);
        int route_lenght(int**, int);
        int get_n_routes(void);
        int get_shortest_time(int, int, int **);
        int get_shortest_time(int, int, int, int **);
        int get_time(int, int, int, int **);
        double calculate_affinity(int, int, long, int);
        bool is_in_route(int, int);
        bool clonal_selection(void);
        bool mutation_process(double);
        bool mutate(int, int);

    private:
        int **sol_m; 
        int n_nodes;
        int n_routes;
        int n_pob;
        int max_routes;
        vector<int> sol;
        vector<int> *sol_group;
        int *top_demand_nodes;
        int *top_demand;
        int top_tam;
        int last_index;
        vector<antigen> antigens;

        //vector<coord> camino;
        //vector<coord> solution;
};
#endif
