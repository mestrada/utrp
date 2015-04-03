#ifndef SOLUTION_H
#define SOLUTION_H

#include <vector>

#include "matrix.h"

#define EMPTY -1
#define INF 100


class solution
{
    public:
        
        solution(int, int, int, unsigned);
        ~solution();
        bool is_in(int, std::vector<int>);
        void print();
        void setDemandMatrix(int **);
        void setTimeMatrix(int **);
        void calculate(void);
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

        bool is_feasible(std::vector<int>*);
        double OperatorCost(std::vector<int>*);
        double PassengerCost(std::vector<int>*);
        void evaluateCosts(void);
};
#endif
