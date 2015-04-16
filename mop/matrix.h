/* 
 * File:   matrix.h
 * Author: titool86
 *
 * Created on October 10, 2010, 7:59 PM
 */

#ifndef MATRIX_H
#define MATRIX_H

#include <string>

using namespace std;

class matrix {
public:
    matrix();
    matrix(int);
    matrix(const matrix& orig);
    ~matrix();
    void load(char *);
    int **getMatrix();
    //virtual ~matrix();
    void print();
private:
    int **m;
    int n_nodes;

};

#endif  /* MATRIX_H */
