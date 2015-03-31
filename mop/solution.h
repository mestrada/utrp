#ifndef SOLUTION_H
#define SOLUTION_H

#include <vector>
using namespace std;

class solution
{
    public:
        
        solution(int, int, int, unsigned);
        ~solution();
        bool is_in(int, vector<int>);
        void print();
    private:
        int pop_size;
        int routes;
        int nodes;
        unsigned seed;

        //int **Q;
        //vector<int> *Q;
        std::vector< std::vector<int> > Q;
        int **Ag;
        int **Ab;
        int **P;
};
#endif
