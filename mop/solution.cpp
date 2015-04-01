#include "solution.h"
//#include "Graph.h"

#include <iostream>
#include <stdlib.h>
#include <algorithm>


/*using namespace std;*/

solution::solution(int population, int n_routes, int n_nodes, unsigned seed): pop_size(population),
routes(n_routes), nodes(n_nodes), seed(seed){

    int rand_route;

    // Solution = Set of N routes

    /*  Q = set of solutions
        Ab = Antibodies (Worst solutions)
        Ag = Antigens (Best solutions)
        P = Pool of solutions to temporary store clones
    */

    Q = std::vector< std::vector<int> >
        (pop_size, std::vector<int>(routes));
    // Ab Antibodies set
    Ab = std::vector< std::vector<int> >
        (pop_size, std::vector<int>(routes));
    // Ag Antigens set
    Ag = std::vector< std::vector<int> >
        (pop_size, std::vector<int>(routes));
    // P pool of clones
    P = std::vector< std::vector<int> >
        (2*pop_size, std::vector<int>(routes));

    // Steps MOAIS-HV Pierrard & Coello

    //  1.- Initialize population
    //      1a. Generate random individuals to fill Q

    std::cout << "Generate random individuals to fill Q" << std::endl;

    for(std::vector< std::vector<int> >::iterator it=Q.begin(); it != Q.end(); ++it){

        for(std::vector<int>::iterator jt=(*it).begin(); jt != (*it).end(); ++jt){

            rand_route = (int) (rand() % nodes);

            // Constraint: No loops and no repetitions        
            if((*it).empty()){
                *jt = rand_route;
            }
            else{
                if(!is_in(rand_route, *it))
                    *jt = rand_route;
            }
        }
    }

    //      1b. Store the best individuals in Ag

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
            std::cout << Q[i][j] << " ";
        }
        std::cout << std::endl;
    }
}