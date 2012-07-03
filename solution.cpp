#include "solution.h"
#include "Graph.h"

#include <iostream>
#include <stdlib.h>

#define length(a) ( sizeof ( a ) / sizeof ( *a ) )

#define MIN_NODE 5

using namespace std;

solution::solution(){ 
}


solution::solution(int sizeX, int sizeY):n_nodes(sizeX), n_routes(sizeY){

    last_index = 0;
    sol_m = new int*[n_nodes];

        for(int i=0; i<n_nodes; i++){
            sol_m[i] = new int[n_routes];
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

    delete [] top_demand_nodes;
    delete [] top_demand;
}

void solution::generate_solution(int **m_tt){
    /*
    int initial_node;
    int current_node;

    for(int i=0; i<n_nodes; i++){
        initial_node = i;
        current_node = i;
        sol_m[i][initial_node] = 0;

        for(int j=i+1; j<n_nodes; j++){
            if(m_tt[current_node][j] > 0){
                sol_m[j][initial_node] = current_node;
                current_node = j;
            }
        
        }
    sol.push_back(i);
    }

    for(int i=n_nodes-1; i>=0; i--){
        initial_node = i + 15;
        current_node = i;
        sol_m[i][initial_node] = 0;

        for(int j=i-1; j>0; j--){
            if(m_tt[current_node][j] > 0){
                sol_m[j][initial_node] = current_node;
                current_node = j;
            }
        }
    sol.push_back(i+15);
    }*/

    for(int index=0; index<top_tam; index++){
        Graph G(n_nodes);
        G.read(m_tt);
        G.set_source(top_demand_nodes[index]);
        G.dijkstra();

        //G.print_deb();
        G.output();
        last_index = G.fill_set(sol_m, last_index, n_routes);
        //G.print_sol();
        //G.~Graph();
        //break;
    }

    for(int i=0; i<last_index; i++)
        if(node_count(i) > MIN_NODE)
            sol.push_back(i);

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
    for ( it=sol.begin() ; it < sol.end(); it++ ){
        cout << *it << ", ";
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
    //cout << "Time by route" << endl;
    for(int j=0; j<n_routes; j++){
        //route_time = 0;
        if(is_in(j)){
            for(int i=0; i<n_nodes; i++){
                total_time += sol_m[i][j] > 0 && (i != sol_m[i][j]) && tt[i][sol_m[i][j]] > 0 ? tt[i][sol_m[i][j]] : 0;
                //route_time += sol_m[i][j] > 0 && (i != sol_m[i][j]) && tt[i][sol_m[i][j]] > 0 ? tt[i][sol_m[i][j]] : 0;
            }
            //cout << "Time for route " << j << ": " << route_time << endl;

        }
        
    }

    return  total_time;
}

long solution::evaluate_time(int **tt, int **td){
    int total_demand = 0;
    long ptd = 0;
    int st = 0;
    for(int i=0; i<n_nodes; i++){
        //if(!is_in(i))
            //continue;

        for(int j=0; j<=i; j++){
            if(td[i][j] <= 0)
                continue;
            st = get_shortest_time(i, j, tt);
            if(st != 0){
                total_demand += td[i][j];
                ptd += (long)  st*td[i][j];
            }
        }
    }


    return (long) ptd/ (long) total_demand;
}

long solution::evaluate_demand(int **tt, int **td){
    /*evaluate_demand()

        Retorna la evaluación respecto a la demanda satisfecha

        input: demand matrix table.}

        returns: double: sum of alls satisfied demands.
    */
/*        
    int total_demand = 0;
    double route_demand = 0;
    //cout << "One-Way Demand by route" << endl;
    for(int j=0; j<n_routes; j++){
        //route_demand = 0;
        if(is_in(j)){
            for(int i=0; i<n_nodes; i++){
                for(int k=0; k<n_nodes; k++){
                    if(i==k)
                        continue;
                    else{
                        total_demand +=  sol_m[k][j] >= 0 && sol_m[k][j] == i ?  td[k][i] : 0; 
                        //route_demand +=  sol_m[k][j] >= 0 && sol_m[k][j] == i ?  td[k][i] : 0; 
                    }
                }
            }
            //cout << "Covered Demand for route " << j << ": " << route_demand << endl;
        }   
    }
    return  total_demand;
    */
    int current_time;
    long total_demand = 0;
    for(int i=0; i<n_nodes; i++){
        for(int j=0; j<=i; j++){
            current_time = 0;
            if(td[i][j] <= 0)
                continue;

            current_time = get_shortest_time(i, j, tt);
            if(current_time > 0)
                total_demand += td[i][j];

        }
    }


    return total_demand*2;
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