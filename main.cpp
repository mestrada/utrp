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


//#define N_NODES 15
#define N_ROUTES 75
#define PENALTY_RATE 2
#define MIN_ROUTES 3
#define N_TOP_DEMAND_NODES 14
#define N_MAX_ROUTES 5

#define P_TIME 0.3
#define P_COST 0.2
#define P_DEMAND 0.5
#define CREATE_PROB 0.6
#define M_PROB 0.01

int N_ITER;
//int N_ROUTES;
int N_RESTART;
unsigned SEED;
int N_NODES;
int TRESHOLD;
int aff_set_size = 5;
int POBLACION;
//TODO: pasar por parámetros.

typedef struct cont
{
    int time;
    int cost;
    int demand;
}result;


/*
valgrind -v --tool=memcheck --leak-check=full ./urtp
Parámetros:  instancia_tiempo instanciat_demanda N_NODES N_ITER N_RESTART POBLACION SEMILLA
Ejemplo:  td1.txt tt1.txt 15 1000 20 30 1234567890
*/

long eval(long, long, long);

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
        N_NODES = atoi(argv[3]);
        N_ITER = atoi(argv[4]);
        N_RESTART = atoi(argv[5]);
        POBLACION = atoi(argv[6]);
        SEED = (unsigned) atoi(argv[7]);
        TRESHOLD = (int) (N_NODES / 4);

    }
     else{
        cout << "Parámetros insuficientes\n";
        return -1;
    }
        
        srand(SEED);

        int fobj_ini;     
        int **m_td, **m_tt;

        matrices input(N_NODES);

        stringstream sstd, sstt;
        
        input.leer(&f1[0], &f2[0]);
      
        solution sol_set(N_NODES, N_ROUTES, POBLACION, N_MAX_ROUTES);

		//cout << "\n<<<<<<Urban Routing Transit Problem >>>>>\n";
		//cout << "Inicializando soluciones\n";
		
        m_td = input.getMatriz1();
        m_tt = input.getMatriz2();
        input.print();

        sol_set.find_top_demand_nodes(N_TOP_DEMAND_NODES, m_td);
        //sol_set.generate_solution(m_tt);      
        sol_set.generate_solution(m_tt, CREATE_PROB);     
        sol_set.generate_antigens(0.05);
        sol_set.print_antigens(); 

        sol_set.print_fact_routes();
        cout << "--------------------0--------------------\n";
        sol_set.print();
        //cout << "Valores del conjunto total" << endl;
        //cout << "Eval total time: " << sol_set.evaluate_time(m_tt, m_td) << endl;
        //cout << "Eval total one-way demand: " << sol_set.evaluate_demand(m_tt, m_td) << endl;
        //cout << "Costo: " << sol_set.evaluate_cost(m_tt) << endl;
        //cout << "FO: " << eval( sol_set.evaluate_time(m_tt, m_td), sol_set.evaluate_cost(m_tt), 
         //   sol_set.evaluate_demand(m_tt, m_td) ) << endl;

        bool criteria = true;

        long current_time;
        long current_demand = 0;
        long current_cost;

        result r_best;
        r_best = {1000000, 1000000, 0};

        long best_fo = 9999999;

        long current_fo;

        int iter = 0, restart = 0;
        int t_routes = 0;

        
        cout << "Clonal Selection and mutation Alghotrithm " << endl;
        while(iter < N_ITER && criteria){

            sol_set.clonal_selection();
            sol_set.mutation_process(M_PROB);

            iter++;
        }



















/**
        cout << "START," << eval( sol_set.evaluate_time(m_tt, m_td), sol_set.evaluate_cost(m_tt), 
         sol_set.evaluate_demand(m_tt, m_td) ) << "," << sol_set.get_n_routes() << endl;

       int iter = 0, restart = 0;
       int t_routes = 0;

        long current_time;
        long current_demand = 0;
        long current_cost;

        result r_best;
        r_best = {1000000, 1000000, 0};

        long best_fo = 9999999;

        long current_fo;

        while(iter < N_ITER && restart < N_RESTART){

            current_cost = sol_set.evaluate_cost(m_tt);
            current_time = sol_set.evaluate_time(m_tt, m_td);
            current_demand = sol_set.evaluate_demand(m_tt, m_td);
            current_fo = eval(current_time, current_cost, current_demand);
            
            if(current_fo < best_fo  && current_fo != 0){
                best_fo = current_fo;

                r_best.time = current_time;
                r_best.cost = current_cost;
                r_best.demand = current_demand;
                t_routes = sol_set.get_n_routes();
            }
            
            if(sol_set.get_n_routes() == MIN_ROUTES){
                sol_set.reset();
                restart++;
            }
            
            while(!sol_set.change_sol());
            iter++;
        }

        //cout << "Best Cost/Time-Sol:";
        //cout << " Tiempo: " << r_best.time;
        //cout << " Costo: " << r_best.cost;
        //cout << " Demanda: " << r_best.demand;
        //cout << " N routes: " << t_routes;
        //cout << " FO: " << best_fo;
        //cout << endl;

        cout << "FINISH," << best_fo << "," << t_routes << endl;

    return 0;**/

    cout << "AAAAAAAAAAAAAAAAAAAAA" << endl;
}


long eval(long ttime, long cost, long demand){
    if (demand !=0)
        return (long) (100*ttime + 10*cost + 1000*(1/(demand)) )/3;
    else
        return 9999999;
    //return ttime*cost*(1/demand);

}