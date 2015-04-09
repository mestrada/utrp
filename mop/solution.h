#ifndef SOLUTION_H
#define SOLUTION_H

#include <vector>

#include "matrix.h"

#define EMPTY -1
#define INF 100


typedef std::vector<int> Route;
typedef std::vector<Route> Routes;
typedef std::vector<Routes> Solutions;
typedef std::vector< std::vector<std::vector<int> > >::iterator SolIter;
typedef std::vector< std::vector<int> >::iterator RoutesIter;
typedef std::vector<int>::iterator RouteIter;


class solution
{
    public:
        
        solution(int, int, int, unsigned);
        ~solution();
        bool is_in(int, std::vector<int>);
        void print();
        void setDemandMatrix(int **);
        void setTimeMatrix(int **);
        void calculate(int);
    private:
        int pop_size;
        int routes;
        int nodes;
        unsigned seed;

        Solutions Q;
        Solutions Ab;
        Solutions Ag;
        Solutions P;

        int** demand_matrix;
        int** time_matrix;
        int** current_time_matrix;
        int** costMatrix;

        bool is_feasible(std::vector<int>*);
        double OperatorCost(Routes);
        double RouteOperatorCost(std::vector<int>*);
        double PassengerCost(std::vector<int>*);
        void setCurrentTimeMatrix(Routes);
        void evaluateCosts(Routes, int);
        void InitializeMatrix(int** &);
        void ResetMatrix(int** &);
        void InitializeCostMatrix(void);
        void ResetCostMatrix(void);
};
#endif

