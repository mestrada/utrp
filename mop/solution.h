#ifndef SOLUTION_H
#define SOLUTION_H

#include <vector>

#include "matrix.h"

#define EMPTY -1
#define INF 100
#define MUTATION_TYPE_DIST 0.5
#define AFFINITY_PREFERENCE 0.7

typedef std::vector<int> Route;
typedef std::vector<Route> Solutions;
typedef std::vector<std::vector<int> >::iterator SolIter;
typedef std::vector< std::vector<int> >::iterator RoutesIter;
typedef std::vector<int>::iterator RouteIter;

typedef struct individual
{
    double ocost; // Operator cost
    double pcost; // Passenger cost
    int index;
}ind;

typedef std::vector<individual> Individuals;
typedef std::vector<individual>::iterator IndIter;

class solution
{
    public:
        
        solution(int, int, int, int, int, double, unsigned, double);
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
        double threshold;

        Individuals current_values;
        Individuals pool_values;
        Individuals ag_values;

        Solutions Q;
        Solutions Ab;
        Solutions Ag;
        Solutions P;

        int** demand_matrix;
        int** time_matrix;
        int** current_time_matrix;
        int** costMatrix;

        void fillAg(void);
        bool is_feasible(std::vector<int>*);
        double OperatorCost(Solutions);
        double RouteOperatorCost(std::vector<int>*);
        double RoutePassengerCost(Route);
        double PassengerCost(Solutions);
        void setCurrentTimeMatrix(Solutions);
        void printAntigens();
        void mutateChange(double);
        void mutateResize(double, int, int);
        individual evaluateCosts(Solutions, int);
        individual evaluateRouteCosts(Route, int);
        bool isDominated(double, double, int);
        std::vector<int> getNonDominated(void);
        void clone(std::vector<int>);
        void DestroyMatrix(int** &m);
        void InitializeMatrix(int** &);
        void ResetMatrix(int** &);
        void InitializeCostMatrix(void);
        void ResetCostMatrix(void);
        //void evaluateAllCosts(void);
        void evaluateAllCosts(Solutions, Individuals &);
        void evaluateAllRouteCosts(Solutions, Individuals &);
        int getNonDominatedByOperatorCost(void);
        int getNonDominatedByPassengerCost(void);
        void clone(int, int);
        void calculateCostMatrix(Solutions);
        std::vector<int> getBestRoutes(void);
};
#endif

