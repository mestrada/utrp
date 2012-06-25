//#include "mapa.h"
#include "matrices.h"
#include "solution.h"
#include "Graph.h"

#include <iostream>
#include <vector>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>

using namespace std;


#define N_NODES 15
#define N_ROUTES 30
#define PENALTY_RATE 2
#define MIN_ROUTES 5

int N_ITER;
//int N_ROUTES;
int N_RESTART;
unsigned SEED;
//int N_NODES;
//TODO: pasar por parámetros.

typedef struct cont
{
    int time;
    int cost;
    int demand;
}result;


/*
valgrind -v --tool=memcheck --leak-check=full ./urtp
Parámetros:  instancia_tiempo instanciat_demanda output N_ITER N_RESTART VECINDARIO SEMILLA
Ejemplo:  td1.txt tt1.txt output 1000 20 30 1234567890
*/

int main (int argc, char **argv)
{
	
	stringstream i1;
    stringstream i2;
    string f1;
    string f2;


    if(argc > 6){
        
        i1 << "./input/" << argv[1];
        i2 << "./input/" << argv[2];

        f1 = i1.str();
        f2 = i2.str();

        N_ITER = atoi(argv[4]);
        N_RESTART = atoi(argv[5]);
        SEED = (unsigned) atoi(argv[7]);

    }
     else{
        cout << "Parámetros insuficientes\n";
        return -1;
    }
        
        srand(SEED);

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

        sol_set.print_fact_routes();
        sol_set.print();
        cout << "Valores del conjunto total";
        cout << "Eval total time: " << sol_set.evaluate_time(m_tt) << endl;
        cout << "Eval total one-way demand: " << sol_set.evaluate_demand(m_td) << endl;
        cout << "N rutas en solución actual: " << sol_set.evaluate_cost(m_tt) << endl;


        Graph G(N_NODES);
        G.read(m_tt);
        G.set_source(1);
        G.dijkstra();
        G.output();

/*
        //HC
        int iter = 0, restart = 0;

        double current_time;
        double current_demand = 0;
        double current_cost;

        vector<int> br_tsol;
        vector<int> br_dsol;
        vector<int> br_osol;


        int best_tsol = 1000000;
        int best_dsol = 0;
        int best_osol = 1000000;

        result r_cost, r_time, r_demand;

        //r_cost = {1000000, 1000000, 0};
        r_cost.time = 1000000;
        r_cost.cost = 1000000;
        r_cost.demand = 0;
        r_time.time = 1000000;
        r_time.cost = 1000000;
        r_time.demand = 0;
        r_demand.time = 1000000;
        r_demand.cost = 1000000;
        r_demand.demand = 0;

        while(iter < N_ITER && restart < N_RESTART){

            current_cost = sol_set.evaluate_cost(m_tt);
            current_time = sol_set.evaluate_time(m_tt);
            current_demand = sol_set.evaluate_demand(m_td);

            if (current_cost < r_cost.cost){
                best_osol = current_cost;
                //br_osol = sol_set.get_current_sol();
                //r_cost = { current_time, current_cost, current_demand};
                r_cost.time = current_time;
                r_cost.cost = current_cost;
                r_cost.demand = current_demand;
            }
            else{
                if(current_cost == r_cost.cost && current_time < r_cost.time){
                    best_osol = current_cost;
                    //r_cost = { current_time, current_cost, current_demand};
                    r_cost.time = current_time;
                    r_cost.cost = current_cost;
                    r_cost.demand = current_demand;
                }

            }            

            if (current_time < r_time.time){
                best_tsol = current_time;
                //br_tsol = sol_set.get_current_sol();
                //r_time = { current_time, current_cost, current_demand};
                r_time.time = current_time;
                r_time.cost = current_cost;
                r_time.demand = current_demand;
            }


            if (current_demand > r_demand.demand){
                best_dsol = current_demand;
                //br_tsol = sol_set.get_current_sol();
                //r_time = { current_time, current_cost, current_demand};
                r_demand.time = current_time;
                r_demand.cost = current_cost;
                r_demand.demand = current_demand;
            }

            cout << "IT " << iter;
            cout << " Tiempo: " << current_time;
            cout << " Costo: " << current_cost;
            cout << " Demanda: " << current_demand << endl;
            //sol_set.print();

            if(sol_set.get_n_routes() == MIN_ROUTES){
                sol_set.reset();
                restart++;
                //best_osol = 1000000;
                //best_tsol = 1000000;
            }
            sol_set.change_sol();
            iter++;
        }

        cout << "Best Cost-Sol:";
        cout << " Tiempo: " << r_cost.time;
        cout << " Costo: " << r_cost.cost;
        cout << " Demanda: " << r_cost.demand;
        cout << endl;

        cout << "Best Time-Sol:";
        cout << " Tiempo: " << r_time.time;
        cout << " Costo: " << r_time.cost;
        cout << " Demanda: " << r_time.demand;
        cout << endl;

        cout << "Best Demand-Sol:";
        cout << " Tiempo: " << r_demand.time;
        cout << " Costo: " << r_demand.cost;
        cout << " Demanda: " << r_demand.demand;
        cout << endl;

*/
    return 0;
}


