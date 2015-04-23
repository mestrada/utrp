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

typedef struct individual
{
    double ocost; // Operator cost
    double pcost; // Passenger cost
}ind;

typedef std::vector<individual> Individuals;
typedef std::vector<individual>::iterator IndIter;

class solution
{
    public:
        
        solution(int, int, int, int, int, double, unsigned);
        ~solution();
        bool is_in(int, std::vector<int>);
        void print();
        void setDemandMatrix(int **);
        void setTimeMatrix(int **);
        void calculate(int);
        void generateRandomIndividuals(void);
        void printActualValues(void);
    private:
        int pop_size;
        int routes;
        int nodes;
        int minlength;
        int maxlength;
        double mutation_prob;
        unsigned seed;

        Individuals current_values;
        Individuals pool_values;

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
        double PassengerCost(Routes);
        void setCurrentTimeMatrix(Routes);
        void printAntigens();
        void mutateChange(double);
        void mutateResize(double, int, int);
        individual evaluateCosts(Routes, int);
        void InitializeMatrix(int** &);
        void ResetMatrix(int** &);
        void InitializeCostMatrix(void);
        void ResetCostMatrix(void);
        //void evaluateAllCosts(void);
        void evaluateAllCosts(Solutions, Individuals &);
        int getNonDominatedByOperatorCost(void);
        int getNonDominatedByPassengerCost(void);
        void clone(int, int);
};
#endif

