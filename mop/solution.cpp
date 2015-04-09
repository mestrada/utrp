
#include <iostream>
#include <stdlib.h>
#include <algorithm>

#include "Graph.h"
#include "solution.h"
#include "matrix.h"


solution::solution(int population, int n_routes, int n_nodes, unsigned seed):
pop_size(population), routes(n_routes), nodes(n_nodes), seed(seed){

    int rand_route;

    // Solution = Set of N routes

    /*  Q = set of solutions
        Ab = Antibodies (Worst solutions)
        Ag = Antigens (Best solutions)
        P = Pool of solutions to temporary store clones
    */

    // Q set of solutions
    Q = Solutions(pop_size, Routes(routes, Route(routes, EMPTY)));

    // Ab Antibodies set
    Ab = Solutions(pop_size, Routes(routes, Route(routes, EMPTY)));

    // Ag Antigens set
    Ag = Solutions(pop_size, Routes(routes, Route(routes, EMPTY)));
    // P pool of clones
    P = Solutions(2 * pop_size, Routes(routes, Route(routes, EMPTY)));

    // Steps MOAIS-HV Pierrard & Coello

    //  1.- Initialize population
    //      1a. Generate random individuals to fill Q

    std::cout << "Generate random individuals to fill Q" << std::endl;

    for(SolIter it=Q.begin(); it != Q.end(); ++it){
        for(RoutesIter jt=(*it).begin(); jt != (*it).end(); ++jt){
            for(RouteIter kt=(*jt).begin(); kt != (*jt).end(); ++kt){

                rand_route = (int) (rand() % nodes);
                // Constraint: No loops and no repetitions        
                if((*jt).empty()){
                    *kt = rand_route;
                }
                else{
                    if(!is_in(rand_route, *jt)){
                        *kt = rand_route;
                    }
                    else{
                        --kt;
                    }
                }
            }
        }
    }
}

solution::~solution(){

    Q.clear();
    Ab.clear();
    Ag.clear();
    P.clear();
}

bool solution::is_in(int node, std::vector<int> sol){
    std::vector<int>::iterator it;

    it = std::find(sol.begin(), sol.end(), node);

    if(it !=sol.end())
        return true;
    else
        return false;

}

void solution::print(){
    std::cout << "Printing Q" << std::endl;

    for(int i=0; i<pop_size; i++){
        for(int j=0; j<routes; j++){
            for(int k=0; k<routes; k++){
                std::cout << Q[i][j][k] << "\t";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

void solution::setDemandMatrix(int** dmatrix){
    demand_matrix = dmatrix;
}

void solution::setTimeMatrix(int** tmatrix){
    time_matrix = tmatrix;
}

void solution::calculate(void){
    //      1b. Store the best individuals in Ag
    Ag = Q;

    for(SolIter it=Q.begin(); it != Q.end(); ++it){
        for(RoutesIter jt=(*it).begin(); jt != (*it).end(); ++jt){
            if(is_feasible(&(*jt))) {

            }
        }
    }


    int ** costMatrix;

    costMatrix = (int**)malloc(nodes*sizeof(int*));
    for(int i =0; i < nodes; i++){
        costMatrix[i] = (int *) malloc(nodes *sizeof (int));
    }

    for(int i=0; i<nodes; i++){
        Graph G(nodes);
        G.read(time_matrix);
        G.set_source(i);
        G.dijkstra();
        
        G.SetPaths();
        
        //G.print_sol();

        G.fill_matrix(costMatrix, i);
    }

    evaluateCosts();



}

bool solution::is_feasible(std::vector<int>* s){

    for(std::vector<int>::iterator it=s->begin(); it != s->end(); ++it){
        if(it != (s->end() -1)){
            if (time_matrix[*it][*(it + 1)] < 0){
                /*std::cout << "Not feasible" << std::endl;*/
                return false;
            }
        }
    }
    /*std::cout << "Feasible" << std::endl;*/
    return true;
}


double solution::PassengerCost(std::vector<int>* s){

    /*
        

    */

    double fo_value = 0.0;


    for (int i = 0; i < nodes; ++i)
    {
        /* code */
        for (int j = 0; j < nodes; ++j)
        {
            /* code */

            // demand_matrix;
            // time_matrix;

        }
    }

    return fo_value;
}


double solution::OperatorCost(void){

    /*
        The operator cost is the sum of the time between nodes of routes.
    */

    double fo_cost;

     fo_cost = 0.0;


    for(SolIter it=Q.begin(); it != Q.end(); ++it){
        for(RoutesIter jt=(*it).begin(); jt != (*it).end(); ++jt){
            fo_cost += RouteOperatorCost(&(*jt));
        }
    }

    std::cout << "Total Operator Cost: " << fo_cost << std::endl;

    return fo_cost;
}


double solution::RouteOperatorCost(std::vector<int>* s){

    double fo_value;

    fo_value = 0.0;

    for(std::vector<int>::iterator it=s->begin(); it != s->end(); ++it){
            if(it != (s->end() -1)){
                if (time_matrix[*it][*(it + 1)] > 0){
                    fo_value += time_matrix[*it][*(it + 1)];
                }else{
                    fo_value += INF;
                }
            }
        }
    std::cout << "Route Operator Cost: " << fo_value << std::endl;

    return fo_value;
}

void solution::evaluateCosts(void){
    
    OperatorCost();  
}