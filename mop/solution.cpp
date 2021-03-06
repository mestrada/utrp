
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

void PairCost(Individuals sorted_set){
    for(IndIter it=sorted_set.begin(); it != sorted_set.end(); ++it){
        std::cout << "Costs = " << it->ocost << "," << it->pcost << std::endl;
    }

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
    ag_values = Individuals(pop_size);

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

    //std::cout << "Generate random individuals to fill Q" << std::endl;

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

void solution::DestroyMatrix(int** &m){
    for(int i =0; i < nodes; i++){
        free(m[i]);
    }
    free(m);
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

    /*std::cout << "CT CHECK 1"<< std::endl;*/
    // InitializeMatrix(current_time_matrix);
    /*std::cout << "CT CHECK 2"<< std::endl;*/
    ResetMatrix(current_time_matrix);

    for(RoutesIter jt=current_routes.begin(); jt != current_routes.end(); ++jt){
        /*std::cout << "CT CHECK 3"<< std::endl;*/
        for(RouteIter kt=(*jt).begin(); (kt != (*jt).end() - 1 && kt != (*jt).end()); ++kt){
            /*std::cout << "CT CHECK 4 --"<< std::endl;*/
            /*std::cout << "values " << *kt << " | " << *(kt+1) << std::endl;
            std::cout << "CM values " << current_time_matrix[*kt][*(kt + 1)] << std::endl;
            std::cout << "TM values " << time_matrix[*kt][*(kt + 1)] << std::endl;
            current_time_matrix[*kt][*(kt + 1)] = time_matrix[*kt][*(kt + 1)];*/
            /*std::cout << "CT CHECK 5"<< std::endl;*/
            current_time_matrix[*(kt + 1)][*kt] = time_matrix[*(kt + 1)][*kt];
        }
    }

}

bool solution::isDominated(double ocost, double pcost, int index){
    for(int i=1; i<pop_size; i++){
        if(i == index)
            continue;
        if(current_values[i].ocost <=  ocost  && current_values[i].pcost < pcost){
                return true;
            }
        if(current_values[i].ocost <  ocost  && current_values[i].pcost <= pcost){
            return true;
        }
    }
    return false;
}

std::vector<int> solution::getNonDominated(void){
    int min_idx = 0;
    std::vector<int> non_dominated;
    for(int i=1; i<pop_size; i++){
        if(!isDominated(current_values[i].ocost, current_values[i].pcost, i))
            non_dominated.push_back(i);
    }
    return non_dominated;
}

void solution::mutateResize(double mutate_prob, int minLen, int maxLen){
    
    double p;
    int mutation_t;

    int rand_route_node;
    vector<int>::size_type rand_idx;
    vector<int>::size_type size_of_v;

    /*std::cout << "CHECK MR 1" << std::endl;*/
    for(SolIter it=P.begin(); it != P.end(); ++it){
        /*std::cout << "CHECK MR 2" << std::endl;*/
        for(RoutesIter jt=(*it).begin(); jt != (*it).end(); ++jt){
            /*std::cout << "CHECK MR 3" << std::endl;*/
            p = (double) rand() / RAND_MAX;
            size_of_v = (*jt).size();
            if(p < mutate_prob){
                /*std::cout << "CHECK MR 4" << std::endl;*/
                /*std::cout << "CHECK MR V " << (*jt).size() <<  std::endl;*/
                if (size_of_v == 0)
                    break;
                rand_idx = rand() % size_of_v;
                if (is_in(rand_idx, *jt)){
                    continue;
                }
                /*std::cout << "CHECK MR 4ii" << std::endl;*/
                mutation_t = (int) rand() % 2;
                if(mutation_t != 1){
                    /*std::cout << "CHECK MR 5a" << std::endl;*/
                    // Add node

                    if(size_of_v == maxLen)
                        continue;

                    rand_route_node = (int) rand() % nodes;
                    //jt->emplace(jt->begin() + rand_idx, rand_route_node);
                    //jt->push_back(rand_route_node);
                    if (rand_idx > 0){
                        if(time_matrix[*((*jt).begin() + rand_idx - 1)][rand_route_node] != EMPTY)
                        {
                            if(rand_idx > size_of_v -1){
                                jt->insert(jt->begin() + rand_idx, rand_route_node);
                            }
                            else{
                                if(time_matrix[*((*jt).begin() + rand_idx)][rand_route_node] != EMPTY){
                                jt->insert(jt->begin() + rand_idx, rand_route_node);
                                }
                            }
                        }
                    }
                    else
                    {
                        if(time_matrix[*((*jt).begin() + rand_idx)][rand_route_node] != EMPTY){
                            jt->insert(jt->begin() + rand_idx, rand_route_node);
                            // std::cout << "CHECK MR 5b rand :" << rand_route_node << std::endl;
                        }
                    }
                }
                else{
                    /*std::cout << "CHECK MR 5b" << std::endl;*/
                    // Remove node
                    if (size_of_v == minLen)
                        continue;

                    if((*jt).begin() + rand_idx == (*jt).begin()){
                                            // if(rand_idx == 0 or rand_idx == (*jt).size() - 1){
                        // std::cout << "idx: " << rand_idx << " size: " << (*jt).size() << " pos: " << *((*jt).begin() + rand_idx) << std::endl;
                        // std::cout << "begin: " << *((*jt).begin()) << std::endl;
                        (*jt).erase((*jt).begin() + rand_idx);

                    }else{
                        if( (*jt).begin() + rand_idx == (*jt).end() -1){
                                                // if(rand_idx == 0 or rand_idx == (*jt).size() - 1){
                        // std::cout << "idx: " << rand_idx << " size: " << " pos: " << *((*jt).begin() + rand_idx) << std::endl;
                        // std::cout << "begin: " << *((*jt).begin()) << std::endl;
                        // std::cout << "end: " << *((*jt).end()) << std::endl;
                        (*jt).erase((*jt).begin() + rand_idx);

                        }
                        else{
                        if(time_matrix[*((*jt).begin() + rand_idx - 1)][*((*jt).begin() + rand_idx + 1)] != EMPTY){
                            // std::cout << "V2 idx: " << rand_idx << " size: " << (*jt).size() << " pos: " << *((*jt).begin() + rand_idx) << std::endl;
                            (*jt).erase((*jt).begin() + rand_idx);
                            }
                        }
                    }
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
                int limit = nodes * 2;
                int idx = 0;
                if(p < mutate_prob){
                    while(idx < limit){
                        rand_route_node = (int) (rand() % nodes);
                        if(!is_in(rand_route_node, *jt))
                        {
                            //std::cout << "Change " << *kt << "for " << rand_route_node << endl;

                            if(kt != (*jt).begin() && kt != (*jt).end() - 1)
                            {
                                if(time_matrix[*(kt-1)][rand_route_node] != EMPTY
                                    && time_matrix[*(kt+1)][rand_route_node] != EMPTY)
                                {
                                    *kt = rand_route_node;
                                }
                            }
                            else
                            {
                                if(kt == (*jt).begin())
                                {
                                    if(time_matrix[*(kt+1)][rand_route_node] != EMPTY)
                                    {
                                        // std::cout << "Change " << *kt << "for " << rand_route_node << endl;
                                        *kt = rand_route_node;
                                    }
                                }
                                else{
                                    if(kt == (*jt).end()){
                                        if(time_matrix[*(kt-1)][rand_route_node] != EMPTY)
                                        {
                                            *kt = rand_route_node;   
                                        }

                                    }
                                }
                            }

                            
                            break;
                        }
                        idx++;
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
    

    // std::cout << "Total Operator Cost: " << fo_cost << std::endl;

    return fo_cost;
}


double solution::RouteOperatorCost(std::vector<int>* s){

    double fo_value;

    fo_value = 0.0;

    for(std::vector<int>::iterator it=s->begin(); it != s->end(); ++it){
            if(it != (s->end() -1)){
                // std::cout << "it- " << *it << " |it+1- " << *(it+1) << std::endl;
                // std::cout << "tt- " << current_time_matrix[*it][*(it + 1)] << std::endl;
                // std::cout << "tt- " << current_time_matrix[*(it + 1)][*it] << std::endl;
                 // if (current_time_matrix[*it][*(it + 1)] > 0 && current_time_matrix[*it][*(it + 1)] < 11){
                if (time_matrix[*it][*(it + 1)] > 0){
                    // if (time_matrix[*it][*(it + 1)] > 1000){
                    //     // std::cout << "ctm " << time_matrix[*it][*(it + 1)] << std::endl;
                    //     std::cout << "x,y: " << *it << "," << *(it + 1) << std::endl;
                    // }
                    
                    fo_value += time_matrix[*it][*(it + 1)];
                }else{

                    // if (time_matrix[*(it + 1)][*it] > 0 && time_matrix[*it][*(it + 1)] < 11){
                    if (time_matrix[*(it + 1)][*it] > 0){
                        // if (time_matrix[*(it + 1)][*it] > 1000){
                        //     // std::cout << "ctm " << time_matrix[*(it + 1)][*it] << std::endl;    
                        //     std::cout << "x,y: " << *(it + 1) << "," << *it << std::endl;
                        // }
                        
                        fo_value += time_matrix[*(it + 1)][*it];
                    }
                    else{
                    fo_value += INF;
                    }
                }
            }
        }

    //std::cout << "Route Operator Cost: " << fo_value << std::endl;

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
    //std::cout << "Evaluate All Costs" << std::endl; 
    int sol_number_aux = 0;
    individual aux_cost;
    for(SolIter it=sol_ref.begin(); it != sol_ref.end(); ++it){
        aux_cost = evaluateCosts(*it, sol_number_aux);
        ind_ref[sol_number_aux] = aux_cost;
        //std::cout << "Sol " << sol_number_aux << ": " << aux_cost.ocost << " | " << aux_cost.pcost << std::endl;
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


void solution::clone(std::vector<int> ndo){
    
    int counter = 0;
    double p;

    while(counter < 2 * pop_size){
        for(std::vector<int>::iterator it=ndo.begin(); it != ndo.end(); ++it){
            p = (double) rand() / RAND_MAX;
            if(p < AFFINITY_PREFERENCE){
                P[counter] = Ag[*it];
            }
            else{
                int alt_idx = rand() % pop_size;
                P[counter] = Ag[alt_idx];
            }


            counter++;
            if(counter >= 2 * pop_size){
                break;
            }
        }
    }
}

void solution::calculateCostMatrix(Solutions sol_set){
    int sol_number = 0;
    /*std::cout << "CM CHECK 1"<< std::endl;*/
    for(SolIter it=sol_set.begin(); it != sol_set.end(); ++it){
        /*std::cout << "CM CHECK 2"<< std::endl;*/
        setCurrentTimeMatrix(*it);
        /*std::cout << "CM CHECK 3"<< std::endl;*/
        for(int i=0; i<nodes; i++){
            Graph G(nodes);
            G.read(current_time_matrix);
            G.set_source(i);
            G.dijkstra();
            G.SetPaths();
            G.fill_matrix(costMatrix, i);
        }
        /*std::cout << "CM CHECK 4"<< std::endl;*/
        //std::cout << "Total cost for sol: " << sol_number << std::endl;
        //evaluateAllCosts();
        //evaluateCosts(*it, sol_number);
        sol_number++;
    }
    /*std::cout << "CM CHECK 5"<< std::endl;*/
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

    InitializeMatrix(current_time_matrix);

    InitializeCostMatrix();

    // current_time_matrix = time_matrix;

    for(int i=0;i<nodes;i++)
    {
        for(int j=0;j<nodes;j++)
        {
            current_time_matrix[i][j] = time_matrix[i][j];
        }
    }

    evaluateAllCosts(Ag, current_values);
    //printActualValues();

    double initial_hv = 0.0;
    double ref_hv = 0.0;

    std::sort(current_values.begin(), current_values.end(), SortbyOperatorReverse);
    initial_hv = HyperVolume(current_values, true);

    double ag_hv;

    // while((iteration < iter) && !( (ref_hv > initial_hv) && (ref_hv - initial_hv) > threshold * initial_hv)) {
    while(iteration < iter) {
        // std::cout << "iter " << iteration << std::endl; 
        ResetCostMatrix();

        //printAntigens();
        /*std::cout << "CHECKPOINT 1"<< std::endl;*/
        calculateCostMatrix(Ag);

        //evaluateAllCosts(Ag, ag_values);
        evaluateAllCosts(Ag, current_values);
        Individuals ag_sol;

        for(int j=0; j<Ag.size(); j++){
            ag_sol.push_back(current_values[j]);
        }

        std::sort(ag_sol.begin(), ag_sol.end(), SortbyOperatorReverse);
        ag_hv = HyperVolume(ag_sol, false);
        /*std::cout << "CHECKPOINT 2"<< std::endl;*/
        
        //DestroyMatrix(current_time_matrix);
        /*std::cout << "CHECKPOINT 3"<< std::endl;*/
        int ndo_idx, ndp_idx;
        std::vector<int> nondom;
        // Non-dominated Operator cost

        nondom = getNonDominated();
        // Non-dominated Passenger cost
        /*std::cout << "CHECKPOINT 4"<< std::endl;*/

        //std::cout << "Non Dominated op & pass:\t" << ndo_idx << " | " << ndp_idx << std::endl;

        if(nondom.size() > 0){
            /*std::cout << "CHECKPOINT 5a"<< std::endl;*/
            clone(nondom);
        }
        else{
            /*std::cout << "CHECKPOINT 5b"<< std::endl;*/
            ndp_idx = getNonDominatedByPassengerCost();
            ndo_idx = getNonDominatedByOperatorCost();
            clone(ndo_idx, ndp_idx);
        }
        /*std::cout << "CHECKPOINT 6"<< std::endl;*/
        //mutate

        //std::cout << "Mutation process executed with p: " << mutation_prob << std::endl;

        mutation_type = (double) rand() / RAND_MAX;
        if(mutation_type < MUTATION_TYPE_DIST){
            /*std::cout << "CHECKPOINT 7a"<< std::endl;*/
            mutateChange(mutation_prob);}
        else{
            /*std::cout << "CHECKPOINT 7b"<< std::endl;*/
            mutateResize(mutation_prob, minlength, maxlength);
        }

        /*std::cout << "CHECKPOINT 8"<< std::endl;*/
        // DestroyMatrix(current_time_matrix);
        /*std::cout << "CHECKPOINT 9"<< std::endl;*/
        calculateCostMatrix(P);
        /*std::cout << "CHECKPOINT 10"<< std::endl;*/
        evaluateAllCosts(P, pool_values);
        /*std::cout << "CHECKPOINT 11"<< std::endl;*/
        // DestroyMatrix(current_time_matrix);
        /*std::cout << "CHECKPOINT 12"<< std::endl;*/
        double max_hv = 0.0;
        double current_hv = 0.0;
        
        int current_idx = 0;
        bool betterThanAg;
        Individuals current_sol;
        Individuals best_sol;

        // ag_values = Individuals(pop_size);
        
        // std::cout << "Print values" << std::endl;
        // std::cout << "Print size: " << ag_values.size() << std::endl;
        // for(int ai = 0; ai < ag_values.size(); ai++){
        //     std::cout << "Pcost: " << (ag_values[ai]).pcost << " Ocost: " << (ag_values[ai]).ocost << std::endl;
        // }
        // std::cout << "Ag: " << ag_hv << std::endl;

        std::vector<int> best_ag;
        std::vector<int> current_ag;
        std::vector<int> temp_current_ag;

        Individuals tempValues = Individuals(pop_size);
        Solutions tempSol = Solutions(pop_size, Routes(routes, Route(routes, EMPTY)));

        tempValues = pool_values;

        int n_sucess = 0;

        for(IndIter it=pool_values.begin(); it!=pool_values.end(); it++){
            // max_hv = 0.0;
            // current_hv = 0.0;
            /*for( std::vector<int>::const_iterator i = best_ag.begin(); i != best_ag.end(); ++i)
                std::cout << *i << ' ';
            std::cout << std::endl;*/

            // betterThanAg = false;

            if(current_ag.size() < pop_size){
                current_ag.push_back(current_idx);
                best_ag = current_ag;
                temp_current_ag = current_ag;
            }
            else{

                best_sol.clear();
                current_sol.clear();

                for(int i=0; i<current_ag.size(); i++){
                    temp_current_ag[i] = current_idx;

                    for(int ia=0; ia<best_ag.size(); ia++){
                        tempSol[ia] = P[tempValues[best_ag[ia]].index];
                    }
                    calculateCostMatrix(tempSol);
                    evaluateAllCosts(tempSol, tempValues);

                    for(int ja=0; ja<best_ag.size(); ja++){
                        best_sol.push_back(tempValues[best_ag[ja]]);
                    }

                    std::sort(best_sol.begin(), best_sol.end(), SortbyOperatorReverse);
                    max_hv = HyperVolume(best_sol, false);
                    // DestroyMatrix(current_time_matrix);

                    for(int ib=0; ib<temp_current_ag.size(); ib++){
                        tempSol[ib] = P[tempValues[temp_current_ag[ib]].index];
                    }
                    calculateCostMatrix(tempSol);
                    evaluateAllCosts(tempSol, tempValues);

                    for(int jb=0; jb<temp_current_ag.size(); jb++){
                        current_sol.push_back(tempValues[temp_current_ag[jb]]);
                    }

                    std::sort(current_sol.begin(), current_sol.end(), SortbyOperatorReverse);
                    current_hv = HyperVolume(current_sol, false);

                    // std::cout << "Max HV: " << max_hv << " | Current HV: " << current_hv << std::endl;

                    if(max_hv > current_hv){
                        current_ag = best_ag;
                        n_sucess++;
                        // DestroyMatrix(current_time_matrix);
                    }
                    else{
                        best_ag = temp_current_ag;
                        // DestroyMatrix(current_time_matrix);
                        break;
                    }
                }

                
            }
            if(n_sucess > 5){
                break;
            }
            current_idx++;
        }
        /*std::cout << "CHECKPOINT 13"<< std::endl;*/
        
        if(max_hv > ag_hv){
            // std::cout << "Max hv: " << max_hv << " Ag hv: " << ag_hv << std::endl;    
            for(int i=0; i<best_ag.size(); i++){
                Ag[i] = P[pool_values[best_ag[i]].index];
            }
        }

        /*std::sort(pool_values.begin(), pool_values.end(), SortbyOperator);

        for(int i=0; i< floor(pop_size / 2); i++){
            Ag[i] = P[pool_values[i].index];
        }

        std::sort(pool_values.begin(), pool_values.end(), SortbyPassenger);

        for(int i=0; i< floor(pop_size / 2); i++){
            Ag[i + floor(pop_size / 2)] = P[pool_values[i].index];
        }*/
        //printAntigens();
        calculateCostMatrix(Ag);
        evaluateAllCosts(Ag, current_values);
        // DestroyMatrix(current_time_matrix);
        std::sort(current_values.begin(), current_values.end(), SortbyOperatorReverse);
        //printActualValues();
        ref_hv = HyperVolume(current_values, false);
        iteration++;

        Individuals().swap(current_sol);
        Individuals().swap(best_sol);
        std::vector<int>().swap(best_ag);
        std::vector<int>().swap(current_ag);
    }


    /*std::cout << "\n--------\n\nFinal Evaluation\n\n--------\n" << std::endl;
    std::cout << "\n--------\n--------\n" << std::endl;
    std::cout << "Ref HV: " << ref_hv << std::endl;*/
    // printAntigens();
    calculateCostMatrix(Ag);
    evaluateAllCosts(Ag, current_values);
    
    std::cout << "Executed iterations: " << iteration << std::endl;
    //printActualValues();

    std::sort(current_values.begin(), current_values.end(), SortbyOperatorReverse);
    HyperVolume(current_values, true);
    // DestroyMatrix(current_time_matrix);
    PairCost(current_values);

    // for(SolIter it=Ag.begin(); it != Ag.end(); ++it){
    //     for(RoutesIter jt=(*it).begin(); jt != (*it).end(); ++jt){
    //         if(is_feasible(&(*jt))) {
    //             std::cout << "Feasible" << std::endl;
    //         }
    //         else{
    //             std::cout << "Not feasible" << std::endl;
    //         }
    //     }
    // }
    

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