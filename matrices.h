/* 
 * File:   matrices.h
 * Author: titool86
 *
 * Created on October 10, 2010, 7:59 PM
 */

#ifndef MATRICES_H
#define MATRICES_H

class matrices {
public:
    matrices();
    matrices(int);
    matrices(const matrices& orig);
    ~matrices();
    void leer(char*, char*);
    int **getMatriz1();
    int **getMatriz2();
    //virtual ~matrices();
    void print();
private:
    int **matriz1;
    int **matriz2;
    int n_nodes;

};

#endif  /* MATRICES_H */
