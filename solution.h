#ifndef SOLUTION_H
#define SOLUTION_H

class solution
{
    public:
        
        solution();
        solution(int sizeX, int sizeY);
        void generate_solution(int **m_tt);
        void print(void);

    private: 
        int **sol_m; 
        int n_nodes;
        int n_routes;

        //vector<coord> camino;
        //vector<coord> solution;
};
#endif
