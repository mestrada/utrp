#ifndef __MAPA_H__
#define __MAPA_H__

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <ctime>
#include <vector>


#define UP 0
#define LEFT 2
#define DOWN 1
#define RIGHT 3

// Parámetros para la generación de la solución inicial aleatoria.
#define P_UP 0.4
#define P_LEFT 0.4
#define P_DOWN 0.1
#define P_RIGHT 0.1

//Parámetros para el Simulated Annealing
#define TEMP 1000
#define TASA_ENFRIAMIENTO 0.999   // 0.8 a 0.99
#define N_ITER 100000

//Parámetro para la función objetivo.
#define PENALTY_RATE 2


using namespace std;

typedef struct cont
{
    int x;
    int y;
}coord;

class mapa
{
    private: 
        int **matriz;
        int nfil, ncol;
        int lpx, lpy;
        vector<coord> camino;
		vector<coord> solution;

    public:
        
        mapa(void);
        mapa(int**, int , int);
        mapa(char *);
        mapa(int, int);

        ~mapa();
        int **get(void);

        int col(void);
        int fil(void);
        
        int get_coord(int , int );

        bool make_movement(double, mapa *); //Realiza los movimientos.

        void set_coord(int, int, int);

        bool moverse(int, mapa*);  //Función que es parte de la generación de la solución inicial.

        bool puede_moverse(int, int, int , mapa *); //Determina si desde un punto dado se puede realizar un cierto movimiento.
	
	    void forzar_movimiento(int, int , int );  //Realiza el movimiento dado. Se usa después de puede_moverse().
	
    	bool escapar(int);   //Determina si durante la generación de la sol. inicial se debe escapar o no cuando se está atrapado.

        bool atrapado(mapa *);  //Determina si la ruta está atrapada por obstáculos.

    	bool is_essential(int, int, mapa *); //Determina si se puede eliminar un punto.
	
        void imprimir(void);

        void imprimir_ruta(void);
		
		void imprimir_archivo(void);
		
		void limpiar_ruta(int, int); //Elimina todas las ocurrencias del punto en la lista debido a ciclos.

        bool is_fact(int**); //Determina si la solución es factible según el mapa de obstáculos.

        int get_last_pos_x(void);

        int get_last_pos_y(void);

        int obj_dist(void);

	    int f_obj(void);
 

};



double calcular_prob(double, double);
double enfriar(double);
void mover(mapa *, mapa *);
double funcion_obj(mapa *);
int escoger_mov(void);


#endif

