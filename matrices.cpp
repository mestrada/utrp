/* 
 * File:   matrices.cpp
 * Author: titool86
 * 
 * Created on October 10, 2010, 7:59 PM
 */

#include "matrices.h"
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

matrices::matrices(){
}

matrices::matrices(int n):n_nodes(n) {
}

matrices::matrices(const matrices& orig) {
}

matrices::~matrices() {   
    
    for (int i=0; i<n_nodes; i++){
        free(matriz1[i]);
        free(matriz2[i]);
    }
    free(matriz1); 
    free(matriz2);
    
}


void matrices::leer(char *arch_td, char *arch_tt){
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
        }
    }
    data.close();
    data.open(arch_tt);
    for(int i= 0; i < n_nodes ; i++){
        for(int j=0; j < n_nodes ; j++){
            data >> matriz2[i][j];
        }
    }
    data.close();
}
void matrices::print(){
    cout << "TD Matrix \n";
    for(int i=0;i<n_nodes;i++){
        for(int j=0;j<n_nodes;j++){
            printf(" %d", matriz1[i][j]);
        }
        printf("\n");

    }
    cout << "TT Matrix\n";
    for(int i=0;i<n_nodes;i++){
        for(int j=0;j<n_nodes;j++){
            printf(" %d", matriz2[i][j]);
        }
        printf("\n");
    }

}
int **matrices::getMatriz1(){
    return matriz1;
}

int **matrices::getMatriz2(){
    return matriz2;
}

