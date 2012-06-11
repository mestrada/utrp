/* 
 * File:   matrix.cpp
 * Author: mestrada
 * 
 * Created on June 10, 2012
 */

#include "matrix.h"
#include <cstdlib>
#include <cstdio>
#include <fstream>

using namespace std;

matrix::matrix(){
}

matrix::matrix(int n):n_nodes(n) {
}

matrix::matrix(int sizeX, int sizeY):dx(sizeX), dy(sizeY){

    allocArrays();
    for (int i = 0; i < dx; i++) {
        for (int j = 0; j < dy; j++) {
            m[i][j] = 0;
        }
    }
}


}
matrix::matrix(const matrix& orig) {
}

matrix::~matrix() {   
    for (int i = 0; i < dx; i++) {
        delete [] m[i];
    }
    delete [] m;
}

matrix::allocArrays()
{
    m = new int*[dx];
    for (int i = 0; i < dx; i++) {
        m[i] = new int[dy];
    }
}

matrix::matrix(const matrix& p)
: dx(p.dx), dy(p.dy) {
    allocArrays();
    for (int i=0; i<dx; ++i) {
        for (int j=0; j<dy; ++j) {
            p[i][j] = m.p[i][j];
        }
    }
}

matrix &matrix::operator=(const matrix& rhs){
    if (this == &m) {
    // avoid self-assignment
        return *this;
    } 
    else {
        if (dx != m.dx || dy != m.dy) {
            this->~Matrix();
            dx = m.dx; dy = m.dy;
            allocArrays();
        }
        for (int i = 0; i < dx; i++) {
            for (int j = 0; j < dy; j++) {
                p[i][j] = m.p[i][j];
            }
        }

    return *this;
}

int &Matrix::operator()(int i, int j) {
    return m[i][j];
}

}

void matrix::leer(char* arch_td, char *arch_tt){
    ifstream data;
    data.open(arch_td);
    matriz1 = (int**)malloc(n_nodes*sizeof(int*));
    matriz2 = (int**)malloc(n_nodes*sizeof(int*));
    for(int i =0; i < n_nodes; i++){
        matriz1[i] = (int *) malloc(n_nodes *sizeof (int));
        matriz2[i] = (int *) malloc(n_nodes *sizeof (int));
    }
    
    for(int i= 0; i < n_nodes ; i++){
        for(int j=0; j < n_nodes ; j++){
            data >> matriz1[i][j];
            printf(" %d", matriz1[i][j]);

        }
    }
    data.close();
    data.open(arch_tt);
    for(int i= 0; i < n_nodes ; i++){
        for(int j=0; j < n_nodes ; j++){
            data >> matriz2[i][j];
            printf(" %d", matriz2[i][j]);

        }
    }
    data.close();
}
void matrix::print(){
     for(int i=0;i<n_nodes;i++){
        for(int j=0;j<n_nodes;j++){
            printf(" %d", matriz1[i][j]);
        }
        printf("\n");
    }

}
int **matrix::getMatriz1(){
    return matriz1;
}

int **matrix::getMatriz2(){
    return matriz2;
}

