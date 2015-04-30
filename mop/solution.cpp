
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <iterator>     // std::vector::emplace

#include "Graph.h"
#include "solution.h"
#include "matrix.h"
#include <math.h> 


double HyperVolume(Individuals sorted_set, bool output){

    long hypervolume = 0;

    int prev_ocost = 3000;

    int refop = 3000;
    int refpas = 100;

    for(IndIter it=sorted_set.begin(); it != sorted_set.end(); ++it){
        hypervolume += (prev_ocost - it->ocost) * (refpas - it->pcost);
        prev_ocost = it->ocost;

    }
    if(output)
        std::cout << "HyperVolume = " << hypervolume << std::endl;

    return hypervolume;
}

bool SortbyOperator(individual i, individual j) { return i.ocost < j.ocost; };
bool SortbyOperatorReverse(individual i, individual j) { return i.ocost > j.ocost; };
bool SortbyPassenger(individual i, individual j) { return i.pcost < j.pcost; };
//bool SortByHyperVolume(individuals i, Individuals j) { return HyperVolume(i) > HyperVolume(j)}

solution::solution(int population, int n_routes, int n_nodes, int minlen,
    int maxlen, double mutation_prob, unsigned seed, double threshold):
    pop_size(population), routes(n_routes), nodes(n_nodes),
    minlength(minlen), maxlength(maxlen), mutation_prob(mutation_prob),
    seed(seed), threshold(threshold){

    int rand_route_node;

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
    // Auxiliar Antigens set
    // current_Ag = Solutions(pop_size, Routes(routes, Route(routes, EMPTY)));

    current_values = Individuals(pop_size);
    pool_values = Individuals(2 * pop_size);

    //current_antigens = Individuals(pop_size);
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

/*    for(int i=0; i<pop_size; i++){
        std::cout << "Solution # " << i << std::endl; 
        for(int j=0; j<routes; j++){
            for(int k=0; k<routes; k++){
                std::cout << Q[i][j][k] << "\t";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }*/
    int sol_n = 0;
    for(SolIter it=Q.begin(); it != Q.end(); ++it){
        std::cout << "Solution # " << sol_n << std::endl;
        for(RoutesIter jt=(*it).begin(); jt != (*it).end(); ++jt){
            for(RouteIter kt=(*jt).begin(); kt != (*jt).end(); ++kt){
                std::cout << *kt << "\t";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
        sol_n++;
    }
}

void solution::printAntigens(){
    std::cout << "Printing Ag" << std::endl;

    int sol_n = 0;
    for(SolIter it=Ag.begin(); it != Ag.end(); ++it){
        std::cout << "Solution # " << sol_n << std::endl;
        for(RoutesIter jt=(*it).begin(); jt != (*it).end(); ++jt){
            if(is_feasible(&(*jt))) {
                std::cout << "Feasible" << std::endl;
            }
            else{
                std::cout << "Unfeasible" << std::endl;
            }
            for(RouteIter kt=(*jt).begin(); kt != (*jt).end(); ++kt){
                std::cout << *kt << "\t";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
        sol_n++;
    }

/*
    for(int i=0; i<pop_size; i++){
        std::cout << "Solution # " << i << std::endl; 
        for(int j=0; j<routes; j++){
            for(int k=0; k<routes; k++){
                std::cout << Ag[i][j][k] << "\t";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }*/
}

void solution::generateRandomIndividuals(void){

    int rand_route_node;
    // Steps MOAIS-HV Pierrard & Coello

    //  1.- Initialize population
    //      1a. Generate random individuals to fill Q

    std::cout << "Generate random individuals to fill Q" << std::endl;

    for(SolIter it=Q.begin(); it != Q.end(); ++it){
        for(RoutesIter jt=(*it).begin(); jt != (*it).end(); ++jt){
            for(RouteIter kt=(*jt).begin(); kt != (*jt).end(); ++kt){

                rand_route_node = (int) (rand() % nodes);
                // Constraint: No loops and no repetitions        
                if((*jt).empty()){
                    *kt = rand_route_node;
                }
                else{
                    // Force feasible solutions
                    if(!is_in(rand_route_node, *jt)){
                        int i;
                        for(i=0; i<5; i++){
                            if(time_matrix[*(kt-1)][rand_route_node] != EMPTY
                                && !is_in(rand_route_node, *jt)){
                                *kt = rand_route_node;
                                break;
                            }
                            rand_route_node = (int) (rand() % nodes);
                        }

                        if(i == 5){
                            int j;
                            for(j=0; j<nodes; j++){
                                if(!is_in(j, *jt) && time_matrix[*(kt-1)][j] != EMPTY){
                                    *kt = j;
                                    break;
                                }
                            }
                            if(j == nodes){
                                (*jt).erase(kt);
                                kt--;
                            }
                        }
                        
                    }
                    else{
                        --kt;
                    }
                }
            }
        }
    }
}

void solution::setDemandMatrix(int** dmatrix){
    demand_matrix = dmatrix;
}

void solution::setTimeMatrix(int** tmatrix){
    time_matrix = tmatrix;
}

void solution::InitializeMatrix(int** &m){
    m = (int**)malloc(nodes*sizeof(int*));
    for(int i =0; i < nodes; i++){
        m[i] = (int *) malloc(nodes *sizeof (int));
    }
}

void solution::ResetMatrix(int** &m){
    for(int i=0; i<nodes; i++){
        for(int j=0; j<nodes; j++){
            m[i][j] = -1;
        }
    }
}

void solution::InitializeCostMatrix(void){
    costMatrix = (int**)malloc(nodes*sizeof(int*));
    for(int i =0; i < nodes; i++){
        costMatrix[i] = (int *) malloc(nodes *sizeof (int));
    }
}

void solution::ResetCostMatrix(void){
    for(int i=0; i<nodes; i++){
        for(int j=0; j<nodes; j++){
        costMatrix[i][j] = -1;
        }
    }
}

void solution::setCurrentTimeMatrix(Routes current_routes){

    InitializeMatrix(current_time_matrix);
    ResetMatrix(current_time_matrix);

    for(RoutesIter jt=current_routes.begin(); jt != current_routes.end(); ++jt){
        for(RouteIter kt=(*jt).begin(); kt != (*jt).end() - 1; ++kt){

            current_time_matrix[*kt][*(kt + 1)] = time_matrix[*kt][*(kt + 1)];
            current_time_matrix[*(kt + 1)][*kt] = time_matrix[*(kt + 1)][*kt];
        }
    }

}

void solution::mutateResize(double mutate_prob, int minLen, int maxLen){
    
    double p;
    int mutation_t;

    int rand_route_node;
    int rand_idx;

    for(SolIter it=P.begin(); it != P.end(); ++it){
        for(RoutesIter jt=(*it).begin(); jt != (*it).end(); ++jt){
            p = (double) rand() / RAND_MAX;

            if(p < mutate_prob){
                rand_idx = rand() % (*jt).size();
                mutation_t = (int) rand() % 2;
                if(mutation_t != 1){
                    // Add node
                    if((*jt).size() == maxLen)
                        continue;

                    rand_route_node = (int) rand() % nodes;
                    //jt->emplace(jt->begin() + rand_idx, rand_route_node);
                    //jt->push_back(rand_route_node);
                    jt->insert(jt->begin() + rand_idx, rand_route_node);
                }
                else{
                    // Remove node
                    if ((*jt).size() == minLen)
                        continue;

                    (*jt).erase((*jt).begin() + rand_idx);
                }
            }
        }

    }
}

void solution::mutateChange(double mutate_prob){

    double p;
    int rand_route_node;
        
    for(SolIter it=P.begin(); it != P.end(); ++it){
        for(RoutesIter jt=(*it).begin(); jt != (*it).end(); ++jt){
            for(RouteIter kt=(*jt).begin(); kt != (*jt).end(); ++kt){
                // Force feasibility
                if(kt != (jt->end()) - 1){
                    if(current_time_matrix[*(kt+1)][*kt] != -1){
                        //&& current_time_matrix[*(kt+1)][*kt] != -1)
                        /*cout << "ctt " << current_time_matrix[*(kt+1)][*kt] << endl;
                        cout << "kt & kt-1 " << *kt << *(kt+1) << endl;*/
                        continue;}
                }
                p = (double) rand() / RAND_MAX;
                if(p < mutate_prob){
                    while(true){
                        rand_route_node = (int) (rand() % nodes);
                        if(!is_in(rand_route_node, *jt)){
                            //std::cout << "Change " << *kt << "for " << rand_route_node << endl;
                            *kt = rand_route_node;
                            break;
                        }
                    }
                    
                    
                        
                }
            }
        }
    }
}


double solution::PassengerCost(Routes current_routes){

    /*
        

    */

    double total_cost = 0.0;
    double total_demand = 0.0;

    setCurrentTimeMatrix(current_routes);

    ResetCostMatrix();

    for(int i=0; i<nodes; i++){
        Graph G(nodes);
        G.read(current_time_matrix);
        G.set_source(i);
        G.dijkstra();
        G.SetPaths();
        G.fill_matrix(costMatrix, i);
    }

    for (int i = 0; i < nodes; ++i)
    {
        for (int j = 0; j < nodes; ++j)
        {
            // SUM of L COSTS * Demand / Total demand

            if (i == j)
                continue;

            if (costMatrix[i][j] == EMPTY){
                total_cost += demand_matrix[i][j] * INF;
            }
            else{
                total_cost += demand_matrix[i][j] * costMatrix[i][j];
            }
            total_demand += demand_matrix[i][j];
        }
    }

    return (double) total_cost / total_demand;
}


double solution::OperatorCost(Routes current_routes){

    /*
        The operator cost is the sum of the time between nodes of routes.
    */

    double fo_cost;

     fo_cost = 0.0;

    for(RoutesIter jt=current_routes.begin(); jt != current_routes.end(); ++jt){
        fo_cost += RouteOperatorCost(&(*jt));
    }
    

    //std::cout << "Total Operator Cost: " << fo_cost << std::endl;

    return fo_cost;
}


double solution::RouteOperatorCost(std::vector<int>* s){

    double fo_value;

    fo_value = 0.0;

    for(std::vector<int>::iterator it=s->begin(); it != s->end(); ++it){
            if(it != (s->end() -1)){
                if (current_time_matrix[*it][*(it + 1)] > 0){
                    fo_value += current_time_matrix[*it][*(it + 1)];
                }else{
                    fo_value += INF;
                }
            }
        }
    std::cout << "Route Operator Cost: " << fo_value << std::endl;

    return fo_value;
}

individual solution::evaluateCosts(Routes current_routes, int number){
    ind ind_eval;
    double op_cost, pass_cost;

    op_cost = OperatorCost(current_routes);
    pass_cost = PassengerCost(current_routes);

    //std::cout << "Solution # " << number << std::endl;
    //std::cout << "Total Operator cost: " << op_cost << std::endl;
    //std::cout << "Total Passenger cost: " << pass_cost << std::endl;

    ind_eval.ocost = op_cost;
    ind_eval.pcost = pass_cost;
    ind_eval.index = number;

    return ind_eval;
}

int solution::getNonDominatedByOperatorCost(void){
    int min_idx = 0;

    for(int i=1; i<pop_size; i++){
        if(current_values[i].ocost < current_values[min_idx].ocost)
            min_idx = i;
    }

    return min_idx;

}

int solution::getNonDominatedByPassengerCost(void){
    int min_idx = 0;

    for(int i=1; i<pop_size; i++){
        if(current_values[i].pcost < current_values[min_idx].pcost)
            min_idx = i;
    }

    return min_idx;

}

void solution::evaluateAllCosts(Solutions sol_ref, Individuals &ind_ref){
    std::cout << "Evaluate All Costs" << std::endl; 
    int sol_number_aux = 0;
    individual aux_cost;
    for(SolIter it=sol_ref.begin(); it != sol_ref.end(); ++it){
        aux_cost = evaluateCosts(*it, sol_number_aux);
        ind_ref[sol_number_aux] = aux_cost;
        std::cout << "Sol " << sol_number_aux << ": " << aux_cost.ocost << " | " << aux_cost.pcost << std::endl;
        sol_number_aux++;
    }
}

void solution::printActualValues(void){
    for(int i=0; i< pop_size; i++)
        std::cout << "Value for sol " << i << ":\t" << current_values[i].ocost << " | " << current_values[i].pcost << std::endl;
}

void solution::clone(int ndo, int ndp){
    
    for(int i=0; i<2 * pop_size; i++){
        if(i< pop_size){
            P[i] = Ag[i];
        }
        else{
            if(i % 2 == 0)
                P[i] = Ag[ndo];
            else
                P[i] = Ag[ndp];
        }
    }
}

void solution::calculate(int iter){
    
    int iteration = 0;
    double mutation_type;

    //      1b. Store the best individuals in Ag
    Ag = Q;


    for(SolIter it=Q.begin(); it != Q.end(); ++it){
        for(RoutesIter jt=(*it).begin(); jt != (*it).end(); ++jt){
            if(is_feasible(&(*jt))) {

            }
        }
    }

    //InitializeMatrix(current_time_matrix);

    InitializeCostMatrix();

    current_time_matrix = time_matrix;

    evaluateAllCosts(Ag, current_values);
    printActualValues();

    double initial_hv = 0.0;
    double ref_hv = 0.0;

    std::sort(current_values.begin(), current_values.end(), SortbyOperatorReverse);
    initial_hv = HyperVolume(current_values, true);

    while((iteration < iter) && !( (ref_hv > initial_hv) && (ref_hv - initial_hv) > threshold * initial_hv)) {
        ResetCostMatrix();

        int sol_number = 0;

        //printAntigens();

        for(SolIter it=Ag.begin(); it != Ag.end(); ++it){
            
            setCurrentTimeMatrix(*it);

            for(int i=0; i<nodes; i++){
                Graph G(nodes);
                G.read(current_time_matrix);
                G.set_source(i);
                G.dijkstra();
                G.SetPaths();
                G.fill_matrix(costMatrix, i);
            }
            //std::cout << "Total cost for sol: " << sol_number << std::endl;
            //evaluateAllCosts();
            //evaluateCosts(*it, sol_number);
            sol_number++;
        }

        evaluateAllCosts(Ag, current_values);

        int ndo_idx, ndp_idx;
        // Non-dominated Operator cost
        ndo_idx = getNonDominatedByOperatorCost();
        // Non-dominated Passenger cost
        ndp_idx = getNonDominatedByPassengerCost();

        //std::cout << "Non Dominated op & pass:\t" << ndo_idx << " | " << ndp_idx << std::endl;


        clone(ndo_idx, ndp_idx);

        //mutate

        //std::cout << "Mutation process executed with p: " << mutation_prob << std::endl;

        mutation_type = (double) rand() / RAND_MAX;
        if(mutation_type < MUTATION_TYPE_DIST)
            mutateChange(mutation_prob);
        else
            mutateResize(mutation_prob, minlength, maxlength);

        evaluateAllCosts(P, pool_values);

        double max_hv = 0.0;
        double current_hv = 0.0;
        int current_idx = 0;
        Individuals current_sol;
        Individuals best_sol;

        std::vector<int> best_ag;
        std::vector<int> current_ag;

        for(IndIter it=pool_values.begin(); it!=pool_values.end(); it++){

            /*for( std::vector<int>::const_iterator i = best_ag.begin(); i != best_ag.end(); ++i)
                std::cout << *i << ' ';
            std::cout << std::endl;*/

            if(current_ag.size() < pop_size){
                current_ag.push_back(current_idx);
                best_ag = current_ag;
            }
            else{
                best_sol.clear();
                current_sol.clear();

                for(int i=0; i<current_ag.size(); i++){
                    current_ag[i] = current_idx;

                    for(int j=0; j<best_ag.size(); j++){
                        best_sol.push_back(pool_values[best_ag[j]]);
                    }
                    for(int j=0; j<current_ag.size(); j++){
                        current_sol.push_back(pool_values[current_ag[j]]);
                    }

                    std::sort(best_sol.begin(), best_sol.end(), SortbyOperatorReverse);
                    std::sort(current_sol.begin(), current_sol.end(), SortbyOperatorReverse);
                    max_hv = HyperVolume(best_sol, false);
                    current_hv = HyperVolume(current_sol, false);

                    //std::cout << "Max HV: " << max_hv << " | Current HV: " << current_hv << std::endl;

                    if(max_hv > current_hv){
                        current_ag = best_ag;
                    }
                    else{
                        best_ag = current_ag;
                        break;
                    }

                }
            }

            current_idx++;
        }

        for(int i=0; i<best_ag.size(); i++){
            Ag[i] = P[pool_values[best_ag[i]].index];
        }

        /*std::sort(pool_values.begin(), pool_values.end(), SortbyOperator);

        for(int i=0; i< floor(pop_size / 2); i++){
            Ag[i] = P[pool_values[i].index];
        }

        std::sort(pool_values.begin(), pool_values.end(), SortbyPassenger);

        for(int i=0; i< floor(pop_size / 2); i++){
            Ag[i + floor(pop_size / 2)] = P[pool_values[i].index];
        }*/
        printAntigens();
        evaluateAllCosts(Ag, current_values);
        std::sort(current_values.begin(), current_values.end(), SortbyOperatorReverse);
        printActualValues();
        ref_hv = HyperVolume(current_values, true);
        iteration++;
    }


    std::cout << "\n--------\n\nFinal Evaluation\n\n--------\n" << std::endl;
    std::cout << "\n--------\n--------\n" << std::endl;
    std::cout << "N iterations: " << iteration << std::endl;
    std::cout << "Ref HV: " << ref_hv << std::endl;
    printAntigens();
    evaluateAllCosts(Ag, current_values);
    printActualValues();

    std::sort(current_values.begin(), current_values.end(), SortbyOperatorReverse);
    HyperVolume(current_values, true);

/*    for(SolIter it=Ag.begin(); it != Ag.end(); ++it){
        for(RoutesIter jt=(*it).begin(); jt != (*it).end(); ++jt){
            if(is_feasible(&(*jt))) {

            }
        }
    }*/
    

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