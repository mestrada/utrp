#ifndef SOLUTION_H
#define SOLUTION_H

#include <vector>

#include "matrix.h"

class solution
{
    public:
        
        solution(int, int, int, unsigned);
        ~solution();
        bool is_in(int, std::vector<int>);
        void print();
        void setDemandMatrix(int **);
        void setTimeMatrix(int **);
    private:
        int pop_size;
        int routes;
        int nodes;
        unsigned seed;

        std::vector< std::vector<int> > Q;
        std::vector< std::vector<int> > Ab;
        std::vector< std::vector<int> > Ag;
        std::vector< std::vector<int> > P;

        int** demand_matrix;
        int** time_matrix;
};
#endif
