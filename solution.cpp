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
