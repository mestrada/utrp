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
      
        solution sol_set(N_NODES, N_ROUTES);

		cout << "\n<<<<<<Urban Routing Transit Problem >>>>>\n";
		cout << "Inicializando soluciones\n";
		
        m_td = input.getMatriz1();
        m_tt = input.getMatriz2();
        input.print();

        sol_set.generate_solution(m_tt);      

        sol_set.print();

        cout << "Eval time: " << sol_set.evaluate_time(m_tt) << endl;
        cout << "Eval demand: " << sol_set.evaluate_demand(m_td) << endl;
  
    return 0;
}


