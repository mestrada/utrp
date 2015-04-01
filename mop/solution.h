#ifndef SOLUTION_H
#define SOLUTION_H

#include <vector>

class solution
{
    public:
        
        solution(int, int, int, unsigned);
        ~solution();
        bool is_in(int, std::vector<int>);
        void print();
    private:
        int pop_size;
        int routes;
        int nodes;
        unsigned seed;

        std::vector< std::vector<int> > Q;
        std::vector< std::vector<int> > Ab;
        std::vector< std::vector<int> > Ag;
        std::vector< std::vector<int> > P;
};
#endif
