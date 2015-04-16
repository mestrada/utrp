/* 
 * File:   matrix.cpp
 * Author: titool86
 * 
 * Created on October 10, 2010, 7:59 PM
 */

#include "matrix.h"
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

matrix::matrix(){
}

matrix::matrix(int n):n_nodes(n) {
}

matrix::matrix(const matrix& orig) {
}

matrix::~matrix() {   
    
    for (int i=0; i<n_nodes; i++){
        free(m[i]);
    }
    free(m);
    
}


void matrix::load(char *input_file){
    ifstream data;
    cout << "Loading file: " << input_file << endl;
    data.open(input_file);
    m = (int**)malloc(n_nodes*sizeof(int*));
    for(int i =0; i < n_nodes; i++){
        m[i] = (int *) malloc(n_nodes *sizeof (int));
    }
    
    for(int i= 0; i < n_nodes ; i++){
        for(int j=0; j < n_nodes ; j++){
            data >> m[i][j];
        }
    }
    data.close();
}
void matrix::print(){
    cout << "Matrix \n";
    for(int i=0;i<n_nodes;i++){
        for(int j=0;j<n_nodes;j++){
            printf("%d\t", m[i][j]);
        }
        printf("\n");

    }
}
int **matrix::getMatrix(){
    return m;
}
