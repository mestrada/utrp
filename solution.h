#ifndef __SOLUTION_H__
#define __SOLUTION_H__

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <ctime>
#include <vector>

using namespace std;

typedef struct cont
{
    int x;
    int y;
}coord;

class solution
{
    private: 
        int **sol_m; 
        int n_nodes;
        int n_routes;

        //vector<coord> camino;
        //vector<coord> solution;

    public:
        
        solution(void);
        solution(int sizeX, int sizeY);
        void generate_solution(int **m_tt);
        void print(void);
};
#endif
