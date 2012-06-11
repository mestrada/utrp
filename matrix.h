/* 
 * File:   matrix.h
 * Author: mestrada
 *
 * Created on June 10, 2012
 */

#ifndef MATRIX_H
#define MATRIX_H

class matrix {

public:
    matrix(int sizeX, int sizeY);
    matrix();
    ~matrix();
    matrix(const matrix& m);
    matrix &operator=(const matrix& rhs);
    int &operator()(int x, int y);



private:
    int dx, dy; // dimensions, dx × dy
    int **m; // pointer to a pointer to a long integer

    void allocArrays(); 

}
#endif  /* MATRIX_H */
