#include "solution.h"
#include "Graph.h"

#include <iostream>
#include <stdlib.h>

#define length(a) ( sizeof ( a ) / sizeof ( *a ) )

#define MIN_NODE 5

using namespace std;

solution::solution(){ 
}


solution::solution(int sizeX, int sizeY, int pob, int max):n_nodes(sizeX), n_routes(sizeY), n_pob(pob),
max_routes(max){

    last_index = 0;
    sol_m = new int*[n_nodes];
    sol_group = new vector<int>[n_pob];

        for(int i=0; i<n_nodes; i++){
            sol_m[i] = new int[n_routes];
//            sol_group = new vector<int>[n_pob];
        }

    for (int i = 0; i < n_nodes; i++) {
        for (int j = 0; j < n_routes; j++) {
            sol_m[i][j] = -1;
        }
    }
}

solution::~solution(){
    for (int i=0; i<n_nodes; i++)
        delete [] sol_m[i];
    delete [] sol_m; 
    delete [] sol_group;

    delete [] top_demand_nodes;
    delete [] top_demand;
}

void solution::generate_solution(int **m_tt){

    for(int index=0; index<top_tam; index++){
        Graph G(n_nodes);
        G.read(m_tt);
        G.set_source(top_demand_nodes[index]);
        G.dijkstra();

        //G.print_deb();
        G.output();
        last_index = G.fill_set(sol_m, last_index, n_routes, 3);
        //G.print_sol();
        //G.~Graph();
        //break;
    }

    for(int i=0; i<last_index; i++)
        if(node_count(i) > MIN_NODE)
            sol.push_back(i);

}

void solution::generate_solution(int **m_tt, double prob){
    double r;
    for(int index=0; index<top_tam; index++){
        if(last_index >= n_routes)
            break;
        Graph G(n_nodes);
        G.read(m_tt);
        G.set_source(top_demand_nodes[index]);
        G.dijkstra();

        //G.print_deb();
        G.output();
        last_index = G.fill_set(sol_m, last_index, n_routes, max_routes);
        G.print_sol();
        //G.~Graph();
        //break;
    }

    int current_routes;

    for(int i=0; i<n_pob; i++){
        current_routes = 0;
        for(int j=0; j<last_index; j++){
            r = (double) rand() / RAND_MAX;
            if(r < prob)
                if(node_count(j) > MIN_NODE){
                    sol_group[i].push_back(j);
                    current_routes++;
                }
            if( current_routes >= max_routes)
                    break;
        }
        
    }

}

void solution::generate_antigens(){

    for(int i=0; i<n_nodes; i++){
        for(int j=0; j<n_nodes; j++){
            antigen to_insert;
            to_insert.start = i; 
            to_insert.end =j;
            antigens.push_back(to_insert);
        }
    }
}

void solution::generate_antigens(double prob){
    double  r;
    for(int i=0; i<n_nodes; i++){
        for(int j=0; j<n_nodes; j++){
            r = (double) rand() / RAND_MAX;
            if(r < prob){
                antigen to_insert;
                to_insert.start = i; 
                to_insert.end =j;
                antigens.push_back(to_insert);
            }
        }
    }
}

void solution::print_antigens(){
    vector<antigen>::iterator it;
    int c = 0;
    cout << "Printing antigens" << endl;
    for(it=antigens.begin(); it < antigens.end(); it++){
        cout << it->start << ", " << it->end << " | ";
        c++;
    }
    cout << " TOTAL: " << c << endl;

}


int solution::node_count(int route){
    int route_lenght = 0;
    for(int i=0; i<n_nodes; i++){
        if(sol_m[i][route] >= 0){
            route_lenght += 1;
        }

    }
    return route_lenght;
}


void solution::find_top_demand_nodes(int tam, int **td){
    top_tam = tam;
    top_demand_nodes = new int[tam];
    top_demand = new int[tam];
    for (int i=0; i<tam; i++){
        top_demand[i] = 0;
        top_demand_nodes[i] = -1;
    }

    for(int i=0; i<n_nodes; i++){
        for(int j=0; j<n_nodes; j++){
            if(td[i][j] > top_demand[0] && is_distinct(i)){
                top_demand[0] = td[i][j];
                top_demand_nodes[0] = i;
                reorder();
            }
        }
    }
    /*
    cout << "Top Demand Nodes ..." << endl;
    for(int i=0; i<tam; i++)
        cout << top_demand_nodes[i] << ", ";
    cout << endl;
    */
}
//(int) length(top_demand_nodes)
void solution::reorder(){
    for(int i=0; i< top_tam - 1; i++){
        if(top_demand[i] > top_demand[i+1]){
            swap(i, i+1);
        }
    }
}

void solution::swap(int a, int b){
    int swap_index, swap_value;

    swap_index = top_demand_nodes[b];
    swap_value = top_demand[b];

    top_demand[b] = top_demand[a];
    top_demand_nodes[b] = top_demand_nodes[a];

    top_demand[a] = swap_value;
    top_demand_nodes[a] = swap_index;
}

bool solution::is_distinct(int index){
    for(int i=0; i<top_tam; i++){
        if(top_demand_nodes[i] == index)
            return false;
    }
    return true;
}

void solution::print_fact_routes(void){
    cout << "Printing Solutions\n";
    for(int j=0; j<n_routes; j++){
            for(int i=0; i<n_nodes; i++){
                cout << sol_m[i][j] << " ";
            }
            cout << endl;
        }
}

void solution::print(void){
    cout << "Printing actual vector solution" << endl;
    vector<int>::iterator it;
    for( int i=0; i < n_pob; i++){
        cout << "sol" << i << ": ";
        for ( it=sol_group[i].begin() ; it < sol_group[i].end(); it++ ){
            cout << *it << ", ";
        }
        cout << endl;
    }
    cout << endl;

}

bool solution::is_in(int value ){
    vector<int>::iterator it;
    for ( it=sol.begin() ; it < sol.end(); it++ ){
        if(value == *it)
            return true;
    }

    return false;
}

bool solution::is_in(int value, int vector_pos ){

    vector<int>::iterator it;

    for ( it=(sol_group[vector_pos]).begin() ; it < (sol_group[vector_pos]).end(); it++ ){
        if(value == *it)
            return true;
    }

    return false;
}


std::vector<int>  solution::get_current_sol(void){
    return sol;
}

int solution::evaluate_time(int **tt){
    /*evaluate_time()

        Retorna la evaluación respecto al tiempo de las rutas

        input: time matrix table.

        returns: double: sum of alls times of the routes.
    */
    int total_time = 0;
    int route_time = 0;
    cout << "Time by route" << endl;
    for(int k=0; k<n_pob; k++){


        for(int j=0; j<n_routes; j++){
            //route_time = 0;
            if(is_in(j, k)){
                for(int i=0; i<n_nodes; i++){
                    total_time += sol_m[i][j] > 0 && (i != sol_m[i][j]) && tt[i][sol_m[i][j]] > 0 ? tt[i][sol_m[i][j]] : 0;
                    //route_time += sol_m[i][j] > 0 && (i != sol_m[i][j]) && tt[i][sol_m[i][j]] > 0 ? tt[i][sol_m[i][j]] : 0;
                }
                cout << "Time for route " << j << ": " << route_time << endl;

            }
            
        }
    }

    return  total_time;
}

long solution::evaluate_time(int **tt, int **td){
    int total_demand = 0;
    long ptd = 0;
    int st = 0;
    for(int k=0; k<n_pob; k++)
    {
        for(int i=0; i<n_nodes; i++){
            if(!is_in(i, k))
                continue;

            for(int j=0; j<=i; j++){
                if(td[i][j] <= 0)
                    continue;
                st = get_shortest_time(i, j, k, tt);
                if(st != 0){
                    total_demand += td[i][j];
                    ptd += (long)  st*td[i][j];
                }
            }
        }
    }

    if(total_demand)
        return (long) ptd/ (long) total_demand;
    else
        return (long) ptd;
}

long solution::evaluate_demand(int **tt, int **td){
    int current_time;
    long total_demand = 0;
    for(int k=0; k<n_pob; k++){

        for(int i=0; i<n_nodes; i++){
            for(int j=0; j<=i; j++){
                current_time = 0;
                if(td[i][j] <= 0)
                    continue;

                current_time = get_shortest_time(i, j, k, tt);
                if(current_time > 0)
                    total_demand += td[i][j];

            }
        }
    }


    return total_demand*2;
}

int solution::get_shortest_time(int start, int end, int ssol, int **tt){
    int min_time = 999;
    int current_time = 999;
    //cout << "AAAA" << endl;
    for(int i=0; i<n_routes; i++){
        if(!is_in(i, ssol))
            continue;
        if(sol_m[end][i] >= 0){
            if(sol_m[start][i] >=0){

                current_time = get_time(start, end, i, tt);

                if(current_time < 0)
                    current_time = get_time(end, start, i, tt);

                if(current_time < 0)
                    cout << "[ERROR] no encuentra conexión entre 2 puntos ----------" << endl;
            }

            if(current_time < min_time){
                min_time = current_time;
            }
        }
    }
       
    if(min_time == 999 | min_time < 0){
        //cout << "MINT" << min_time << endl;
        return 0;
    }

    return min_time;
}


long solution::evaluate_cost(int **tt){
    vector<int>::iterator it;
    long total_cost = 0;
    for ( it=sol.begin() ; it < sol.end(); it++ ){
        total_cost += route_lenght(tt, *it);
    }
    return total_cost;
}

bool solution::change_sol(void){
    double r;
    vector<int>::iterator it;
    r = (double) rand() / RAND_MAX;

    if(r < 0.50){
        for ( it=sol.begin() ; it < sol.end(); it++ ){
            r = (double) rand() / RAND_MAX;

            if(r < 0.10){
                sol.erase(it);
                return true;
            }
        }
    }
    else{
        for(int i=0; i<n_routes; i++){
            r = (double) rand() / RAND_MAX;
            if(!is_in(i)){
                if(r < 0.10){
                    sol.push_back(i);
                    return true;
                }
            }
        }

    }

    return false;
}

void solution::reset(void){
    for(int i=0; i<n_routes; i++){
        if(!is_in(i))
            sol.push_back(i);
    }
}

int solution::route_lenght(int **tt, int route){
    int route_lenght = 0;
    for(int i=0; i<n_nodes; i++){
        if(sol_m[i][route] > 0){
            route_lenght += tt[i][sol_m[i][route]];
        }

    }
    return route_lenght;
}

int solution::get_n_routes(void){
    return sol.size();
}

int solution::get_shortest_time(int start, int end, int **tt){
    int min_time = 999;
    int current_time = 999;
    //cout << "AAAA" << endl;
    for(int i=0; i<n_routes; i++){
        if(!is_in(i))
            continue;
        if(sol_m[end][i] >= 0){
            if(sol_m[start][i] >=0){

                current_time = get_time(start, end, i, tt);

                if(current_time < 0)
                    current_time = get_time(end, start, i, tt);

                if(current_time < 0)
                    cout << "[ERROR] no encuentra conexión entre 2 puntos ----------" << endl;
            }

            if(current_time < min_time){
                min_time = current_time;
            }
        }
    }
       
    if(min_time == 999 | min_time < 0){
        //cout << "MINT" << min_time << endl;
        return 0;
    }

    return min_time;
}



int solution::get_time(int i, int j, int route, int **tt){
    if(i!=j && sol_m[j][route] != j && sol_m[j][route] >=0){
        return tt[j][sol_m[j][route]] + get_time(i, sol_m[j][route], route, tt );
    }
    else{
        if(i == j)
            return 0;
        if(sol_m[j][route] == j)
            return -1000000;
        if(sol_m[j][route] < 0 )
            return -1000000;
    }

}

bool solution::is_in_route(int index, int route){
    for(int i=0; i<n_nodes; i++){
        if (sol_m[i][route] >= 0)
            return true;
    }

    return false;
}


double solution::calculate_affinity(int x, int y, long dda, int route){
    int count;

    count = int(is_in_route(x, route)) + int(is_in_route(y, route));

    switch (count)
    {
        case 0: 
            return 0;
            break;
        case 1: 
            return dda/2;
            break;
        case 2: 
            return dda;
            break;

        default: 
            return 0;
            break;
}

}

bool solution::clonal_selection(){
    vector<antigen>::iterator it;   
    vector<int>::iterator it2;
    int c = 0;
    for(it=antigens.begin(); it<antigens.end(); it++){
        //cout << "Affinity A" << c << endl;
        for( int i=0; i < n_pob; i++){
            // For each solution calculate the affinity.
            for ( it2=sol_group[i].begin() ; it2 < sol_group[i].end(); it2++ ){
            //    cout << calculate_affinity(it->start, it->end, 1234, *it2 ) << " | ";
            }
    }
    c++;
    //cout << endl;

    }
    return true;
}

bool solution::mutation_process(double m_prob){
    double r;
    int counter;
    vector<int>::iterator it;
    std::vector<mut>::iterator it2;
    vector<mut> stm;
    mut aux;

    //cout << "Size of sol_group: " << sol_group->size() << endl;
    for( int i=0; i < n_pob; i++){
        counter = 0;
        for ( it=sol_group[i].begin() ; it < sol_group[i].end(); it++ ){
            //cout  << *it << ", " << counter << " | " ;

            r = (double) rand() / RAND_MAX;
            if (r < m_prob){
                /*
                aux.ind = i;
                aux.ite= counter;
                stm.push_back(aux);*/
                mutate(i, counter);
                break;
            }
            counter++;
        }
    }
/*
    for ( it2=stm.begin(); it2 < stm.end(); it2++ ){
        mutate(it2->ind, it2->ite);
    }
*/
}

bool solution::mutate(int ind, int pos){
    double r;
    vector<int>::iterator it;
    int op;
    int route;
    r = (int) (rand()%11);

    /*  3 ways of mutation:
            Add new route
            Delete an existent route
            Change an existent route
    */

    if(r<3)
        op = 0;
    else
        if(r<6)
            op = 1;
        else
            op = 2;
       
    switch(op){
        case 0:
            route = rand()%last_index;
            //cout << endl << "CASE 0 add";
            //cout << "ind, pos: " << ind << ", " << pos ;
            if(!is_in(route, ind))
                sol_group[ind].push_back(route);
            else
                return false;
            break;
        case 1:
            if(sol_group[ind].empty())
                return false;
            //sol_group[ind].erase( pos);
            sol_group[ind].pop_back();
            //cout << endl << "CASE 1 delete";
            break;
        case 2:
            route = rand()%last_index;
            if(sol_group[ind].empty() && pos >= sol_group[ind].size())
                return false;

            if(!is_in(route, ind)){
                //sol_group[ind].erase(sol_group[ind].begin() +pos);
                (sol_group[ind]).at(pos)= route;
                //cout << endl << "CASE 2 change";
                //sol_group[ind].push_back(route);
                return true;
            }
            return false;
            break;

    }

    return true;
}