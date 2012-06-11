//#include "mapa.h"
#include "matrices.h"
#include "solution.h"

#include <iostream>
#include <vector>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

//Parámetros para Hill Climbing
#define N_ITER 100000

//Parámetro para la función objetivo.
#define N_ROUTES 20
#define PENALTY_RATE 2
#define N_NODES 15
//TODO: pasar por parámetros.

int main (int argc, char **argv)
{
	
	srand((unsigned)time(0));
    if(argc > 1){

    }
     else{
        cout << "Parámetros insuficientes\n";
        //return -1;
    }
        
        int fobj_ini;     
        int **m_td, **m_tt;

        matrices input(N_NODES);
        
        input.leer("./input/td1.txt","./input/tt1.txt");

/*
        int **sol_routes;

        sol_routes = new int*[N_NODES];

        for(int i=0; i<N_NODES; i++){
            sol_routes[i] = new int[N_ROUTES];
        }

        for(int i=0; i<N_NODES; i++){
            for(int j=0; j<N_ROUTES; j++){
                sol_routes[i][j] = -1;
            }
        }
*/      
        solution sol_set(N_NODES, N_ROUTES);

		cout << "\n<<<<<<Urban Routing Transit Problem >>>>>\n";
		cout << "Inicializando soluciones\n";
		
        m_td = input.getMatriz1();
        m_tt = input.getMatriz2();
        input.print();

        //sol_routes[i];
/*        
        int initial_node;
        int current_node;

        for(int i=0; i<N_NODES; i++){
            initial_node = i;
            current_node = i;
            sol_routes[i][initial_node] = i;

            for(int j=i+1; j<N_NODES; j++){
                if(m_tt[current_node][j] > 0){
                    sol_routes[j][initial_node] = current_node;
                    current_node = j;
                }
            
            }
        }
*/
        sol_set.generate_solution(m_tt);      

        sol_set.print();
/*        cout << "First Solution\n";

        for(int j=0; j<N_ROUTES; j++){
            for(int i=0; i<N_NODES; i++){
                cout << sol_routes[i][j] << " ";
            }
            cout << endl;
        }
*/
        /*
		while(  iter < N_ITER)
		{
			
			if(sol.make_movement(temp, &m1))
			{
				temp = enfriar(temp);
				esc = 0;
			}
			else
			{
				esc++;
			}
			
			if(esc > 10)
			{
				esc=0;
				temp = enfriar(temp);
			}
			
			iter++;
		}
		
		
		if(m1.is_fact(sol.get()))
		{
			
			sol.imprimir_archivo();
			cout << "\nF_OBJ de la solución final: " << sol.f_obj() << "\n";
			cout << "Diferencia: " << fobj_ini - sol.f_obj() << "\n";
		}
        */
  
    return 0;
}


