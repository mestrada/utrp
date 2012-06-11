#include "mapa.h"


using namespace std;

mapa::mapa(void)
{
    nfil = 0; 
    ncol = 0;
    matriz = NULL;   
    lpx=0;
    lpy=0; 
}

mapa::mapa(int **m, int nf , int nc)
{
    matriz = m;
    nfil = nf;
    ncol = nc;    
    lpx=0;
    lpy=0;
}

mapa::mapa(char *path)
{
    int i, j;
    char linea[1000];
    char *p;
    lpx=0;
    lpy=0;

    matriz = NULL;

    fstream arch;
    arch.open( path , ios::in);

    if(arch.getline(linea, 1000) != NULL)
    {
        p = strtok (linea," ");

        if(p != NULL)
        {       
            ncol = atoi(p);
            p = strtok (NULL, " ,");
        }

        if(p != NULL)
        {            
            nfil = atoi(p);

            p = strtok (NULL, " ,");
        }
    }
    
    matriz = new int*[nfil];
    
    for (i=0; i<nfil; i++)
        matriz[i] = new int[ncol];

    i=0;

    while(arch.getline(linea, 1000) != NULL)
    {
        char *p;
        j=0;
  
        p = strtok (linea," ");

        while (p != NULL)
        {        
            matriz[i][j] = atoi(p);
            p = strtok (NULL, " ,");
            j++;
        }
 
        i++;
    }
    arch.close();
}

mapa::mapa(int nf, int nc)
{
    int i, j;
    nfil = nf;
    ncol = nc;
    matriz = new int*[nf];
    
    for (i=0; i<nf; i++)
        matriz[i] = new int[nc];

    for(i=0; i<nf; i++)
        for(j=0; j<nc; j++)
            matriz[i][j] = 0;

    matriz[nf-1][nc-1] = 1;
    lpx = nc-1;
    lpy = nf-1;
    coord aux = {lpx, lpy};
    camino.push_back(aux);
}

mapa::~mapa()
{   
    int i;
    for (i=0; i<nfil; i++)
        delete [] matriz[i];

    delete [] matriz; 
}
int mapa::col(void)
{  
    return ncol;
}
int mapa::fil(void)
{
    return nfil;
}

int mapa::get_coord(int x, int y)
{
    return matriz[x][y];
}

void mapa::set_coord(int x, int y, int val)
{
    matriz[x][y] = val;
}

void mapa::imprimir(void)
{
    int i,j;
    for(i=0; i< nfil; i++)
    {
        for(j=0; j< ncol; j++)
             cout << matriz[i][j];
        cout << "\n";    
    }
}

void mapa::imprimir_ruta(void)
{
    vector<coord>::iterator it;
    for ( it= camino.begin(); it != camino.end(); it++)
        cout <<"(" << it->x << "," << it->y << ") ";
	
}

void mapa::imprimir_archivo(void)
{
	ofstream fileout;
	fileout.open("solucion.txt", ios::out);
	
	int i,j;
	
	fileout << ncol << " " << nfil << "\n"; 
	fileout << f_obj() << "\n";
	
    for(i=0; i< nfil; i++)
    {
        for(j=0; j< ncol; j++)
			fileout << matriz[i][j] << " ";
		
        fileout << "\n";    
    }
    
	fileout.close();

}

void mapa::limpiar_ruta(int x, int y)
{
	int i;
	bool find = true;
	
	vector<coord>::size_type sz;
	
	while(find)
	{
		sz = camino.size();
		find = false;
		for (i=0; i<sz; i++)
		{
			if(camino[i].x == x && camino[i].y == y)
			{
				camino.erase(camino.begin() + i);
				find = true;
				break;
			}
	
		}
	}
	
}

int **mapa::get(void)
{
    return matriz;
}


bool mapa::is_fact(int** sol)
{
    int i, j;
    for(i=0; i<nfil; i++)
        for(j=0; j<ncol; j++)
            if(matriz[i][j] + sol[i][j] > 1)
                return false;

    return true;
}

int mapa::get_last_pos_x(void)
{
    return lpx;
}

int mapa::get_last_pos_y(void)
{
    return lpy;
}



bool mapa::moverse(int dir, mapa *map)
{
    bool move = false;    
    switch(dir)
    {
        case UP:    if(lpy-1 >= 0)
                    {
                        move = map->matriz[lpy-1][lpx] == 0 && escapar(matriz[lpy-1][lpx]) ?  (matriz[lpy-1][lpx] = 1) : false;
                        if(move)
                            lpy = lpy-1;
                    }
                    break;              
      
        case LEFT:  if( (lpx-1) >= 0)
                    {
                        move = map->matriz[lpy][lpx-1] == 0 && escapar(matriz[lpy][lpx-1]) ?  (matriz[lpy][lpx-1] = 1) : false;
                        if(move)
                            lpx = lpx-1;
                    }
                    break;  

        case DOWN:  if(lpy+1 < nfil-1)
                    {
                        move = map->matriz[lpy+1][lpx] == 0 && escapar(matriz[lpy+1][lpx]) ?  (matriz[lpy+1][lpx] = 1) : false;
                        if(move)
                            lpy= lpy+1;
                    }                    
                    break;  

        case RIGHT: if(lpx+1 < ncol-1)
                    {
                        move = map->matriz[lpy][lpx+1] == 0 && escapar(matriz[lpy][lpx+1]) ?  (matriz[lpy][lpx+1] = 1) : false;
                        if(move)
                            lpx = lpx+1;
                    }
                    break;  
        default:    return false;        
    }    
    
    if(move)
    {
        coord aux = {lpx, lpy};
        camino.push_back(aux);
    }

    return move;
}

bool mapa::puede_moverse(int x, int y, int dir, mapa *map)
{
    bool move = false;        
    switch(dir)
    {
        case UP:    if(y-1 >= 0)
                    {
                        move = map->matriz[y-1][x] == 0 && !(matriz[y-1][x]) ?  true : false;    
                    }
                    break;              
      
        case LEFT:  if( (x-1) >= 0)
                    {
                        move = map->matriz[y][x-1] == 0 && !(matriz[y][x-1]) ?  true : false;
                    }
                    break;  

        case DOWN:  if(y+1 < nfil-1)
                    {
                        move = map->matriz[y+1][x] == 0 && !(matriz[y+1][x]) ?  true : false;
                    }                    
                    break;  

        case RIGHT: if(x+1 < ncol-1)
                    {
                        move = map->matriz[y][x+1] == 0 && !(matriz[y][x+1]) ?  true : false;
                    }
                    break;  
        default:    return false;        
    }    
    
    return move;
}

void mapa::forzar_movimiento(int dir, int x, int y)
{
	coord aux;
	switch(dir)
	{
		
		case UP:    matriz[y-1][x] = 1;
					aux.x =x;
					aux.y = y-1;
					camino.push_back(aux);
					solution.push_back(aux);
					break;              
		
		case LEFT:  matriz[y][x-1] = 1;
					aux.x =x-1;
					aux.y = y;
					camino.push_back(aux);
					solution.push_back(aux);
					break;  

		case DOWN:  matriz[y+1][x] = 1;   
					aux.x =x;
					aux.y = y+1;
					camino.push_back(aux);
					solution.push_back(aux);
					break;  

		case RIGHT: matriz[y][x+1] = 1;
					aux.x =x+1;
					aux.y = y;
					camino.push_back(aux);
					solution.push_back(aux);
						
	}
  
}
bool mapa::escapar(int a)
{
    switch(a)
    {
        case 0:     return true;
                    break;

        case 1:     double r = (double) rand() / RAND_MAX;
                    if(r < 0.1)
                        return true;
                    else
                        return false;
                   break;
           
    }
    return true;
}

bool mapa::atrapado(mapa *m)
{          

    int around = 0, lim = 0;
    bool up, down, left, right;

    right   =   (lpx+1 < ncol)  ?   (around += matriz[lpy][lpx+1] + m->matriz[lpy][lpx+1]) : false; 
    left    =   (lpx-1 >= 0)    ?   (around += matriz[lpy][lpx-1] + m->matriz[lpy][lpx-1]) : false;
    up      =   (lpy-1 >= 0)    ?   (around += matriz[lpy-1][lpx] + m->matriz[lpy-1][lpx]) : false;
    down    =   (lpy+1 < nfil)  ?   (around += matriz[lpy+1][lpx] + m->matriz[lpy+1][lpx]) : false;

    lim = (int)right + (int)left + (int)up + (int)down;
    
    if(around > lim-1)
    {
        return true;
    }
    else
        return false;

}

bool mapa::make_movement(double t , mapa *map)
{
    
    int r = (int) rand() % camino.size();
	
    if(is_essential(camino[r].x, camino[r].y, map))
    {
        int i =0;
		if(t > TEMP/2)
		{
			while(i<10)
			{
				int m = (int) rand() % 4;
				double p = (double) rand() / RAND_MAX;
				if(puede_moverse(camino[r].x, camino[r].y, m, map))
				{
					if(p < calcular_prob(1 ,t))
					{
						forzar_movimiento(m, camino[r].x, camino[r].y);
						return true;
						break;
					}
					else
					{
						return false;
						break;
					}
				}	
					
				i++;
			}
		}
		
		return false;
    }
    else
    {
        matriz[camino[r].y][camino[r].x] = 0;
		limpiar_ruta(camino[r].x, camino[r].y);
		return true;
    }
    
    
}


bool mapa::is_essential(int x, int y, mapa* map)
{
    bool up, down, left, right;
    int u, d , l , r ;
    u=0; d=0; l=0; r=0;
	
	
	if((x+1 < ncol) )
	{
		if(map->matriz[y][x+1] == 0)
		{
			r = matriz[y][x+1];
			right = true;
		}
		else 
			right = false;
	}
	else
		right = false;
	
	
	if(x-1 >= 0)
	{
		if(map->matriz[y][x-1] == 0 )
		{
			l = matriz[y][x-1];
			left = true;
		}
		else
			left = false;
	}
	else
		left = false;
	
	
	if(y-1 >= 0)
	{
		if(map->matriz[y-1][x] == 0 )
		{
			u = matriz[y-1][x];
			up= true;
		}
		else
			up = false;
	}
	else
		up = false;
	
	
	if((y+1 < nfil) )
	{
		if(map->matriz[y+1][x] == 0)
		{
			d = matriz[y+1][x];
			down= true;
		}
		else 
			down = false;
	}
	else
		down = false;
	
	
	
    if(  u+d+l+r == 0 )
    {
		return false;
    }
    else
    {
		if(u+d+l+r == 1)
		{
			if(x!=0 || y !=0  )
				if(x!=lpx || y!=lpy)
					if(x!=ncol-1 || y!=nfil-1)
					{						
						return false;
						
					}
		}
		else
		{
		
			if(up && right && left && down)
			{
				if( (((u && l) ? matriz[y-1][x-1] : 0) + ((u && r) ? matriz[y-1][x+1] : 0) + ((d && l) ? matriz[y+1][x-1] : 0) +((d && r) ? matriz[y+1][x+1] : 0) ) >= u+l+d+r-1 )
				{
					return false;
				}
				else
				{
					return true;
				}
				
			}
		}
    }
  
	 
  
	return true;
}


int mapa::obj_dist(void)
{
    return lpy + lpx;
}

int mapa::f_obj(void)
{
  return PENALTY_RATE*obj_dist() + camino.size();
}



/* *****************************************************************************

            Funciones y estructuras para la inicialización
       ***********************************************************
*/

int escoger_mov(void)
{
    double r = (double) rand() / RAND_MAX;
    
    if( r <= P_UP)
    {
        return UP;
    }
    else
    {
        if( r <= (P_UP + P_LEFT))
        {
            return LEFT;
        }
        else
        {
            if( r <= (P_UP + P_LEFT + P_DOWN))
            {
                return DOWN;
            }
            else
            {
                return RIGHT;
            }
        }
    }
}


void mover(mapa *map, mapa *sol)
{
    int i = 0;
    while( ((sol->get_last_pos_x() !=0 || sol->get_last_pos_y() != 0)) && i<10000 )
    {
        while(!(sol->moverse(escoger_mov(), map)))
        {
            //Intentar moverse hasta que lo haga
        }
        if(sol->atrapado(map))
            i++;
    }
}



/* *****************************************************************************

            Funciones y estructuras para el Modelo
       ***********************************************************
*/

double funcion_obj(mapa * sol)
{
    return sol->obj_dist();
}



/* ***************************************************************************** 

            Funciones y estructuras para Simulated Annealing. 
       ***********************************************************
*/



double calcular_prob( double delta_obj, double temp )
{   
    return  (double) exp((-1)*(double)delta_obj/temp);
}

double enfriar(double temp)
{
    return (double) temp*TASA_ENFRIAMIENTO;
}
