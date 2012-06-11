#include "solution.h"

#include <iostream>

using namespace std;

solution::solution(){ 
}


solution::solution(int sizeX, int sizeY):n_nodes(sizeX), n_routes(sizeY){

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



void solution::generate_solution(int **m_tt){
    int initial_node;
    int current_node;

    for(int i=0; i<n_nodes; i++){
        initial_node = i;
        current_node = i;
        sol_m[i][initial_node] = i;

        for(int j=i+1; j<n_nodes; j++){
            if(m_tt[current_node][j] > 0){
                sol_m[j][initial_node] = current_node;
                current_node = j;
            }
        
        }
    }
}

void solution::print(void){
    cout << "Printing Solutions\n";
    for(int j=0; j<n_routes; j++){
            for(int i=0; i<n_nodes; i++){
                cout << sol_m[i][j] << " ";
            }
            cout << endl;
        }
}


double solution::evaluate_time(int **tt){
    double total_time = 0;

    for(int j=0; j<n_routes; j++){
        for(int i=0; i<n_nodes; i++){
            total_time += sol_m[i][j] > 0 && (i != sol_m[i][j]) && tt[i][sol_m[i][j]] > 0 ? tt[i][sol_m[i][j]] : 0;
        }
    }

    return  total_time;
}

double solution::evaluate_demand(int **td){
    /*evaluate_demand()

        Retorna la evaluación respecto a la demanda satisfecha

        input: demand matrix table.}

        returns: double: sum of alls satisfied demands.
    */


    double total_demand = 0;

    for(int j=0; j<1; j++){
        for(int i=0; i<n_nodes; i++){
            for(int k=0; k<n_nodes; k++){
                if(i==k)
                    continue;
                else{
                    total_demand +=  sol_m[i][j] == k && sol_m[i][j] >= 0 ?  td[i][k] : 0; 
                }
            }
            
        }
    }

    return  total_demand;
}