
#include <iostream>
#include <vector>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>

#include "matrix.h"
/*
#include "solution.h"
#include "Graph.h"
*/

int main (int argc, char **argv)
{
	
	std::stringstream input_path_1;
    std::stringstream input_path_2;
    std::string file_path_demand_matrix;
    std::string file_path_time_matrix;

    int number_routes;
    int n_nodes;
    int min_length;
    int max_length;
    int max_iter;
    int population;
    double mutation_prob;
    unsigned seed;

    /*
    Arguments parsing
        instance_name n_nodes n_routes min_len max_len mutation_prob max_iter population seed

        1 15 5 2 10 0.3 100 10 12345
    */
    if(argc > 8){

        input_path_1 << "./input/td" << argv[1] << ".txt";
        input_path_2 << "./input/tt" << argv[1] << ".txt";

        n_nodes = atoi(argv[2]);
        number_routes = atoi(argv[3]);
        min_length = atoi(argv[4]);
        max_length = atoi(argv[5]);
        mutation_prob = (double) atof(argv[6]);
        max_iter = atoi(argv[7]);
        population = atoi(argv[8]);
        seed = (unsigned) atoi(argv[9]);

        file_path_demand_matrix = input_path_1.str();
        file_path_time_matrix = input_path_2.str();

        std::cout << "================================" << std::endl;
        std::cout << "Using the following parameters:" << std::endl;
        std::cout << "--------------------------------" << std::endl;
        std::cout << "Instance Name:\t\t" << argv[1] << std::endl;
        std::cout << "Number of Nodes:\t" << argv[2] << std::endl;
        std::cout << "Number of Routes:\t" << argv[3] << std::endl;
        std::cout << "Minimum route lenght:\t" << argv[4] << std::endl;
        std::cout << "Maximum route lenght:\t" << argv[5] << std::endl;
        std::cout << "Mutation Probability:\t" << argv[6] << std::endl;
        std::cout << "Maximum iterations:\t" << argv[7] << std::endl;
        std::cout << "Popiulation size:\t" << argv[8] << std::endl;
        std::cout << "Seed:\t\t\t" << argv[9] << std::endl;

    }
     else{
        std::cout << "Parámetros insuficientes" << std::endl;
        return -1;
    }

    matrix demand_matrix(n_nodes);
    demand_matrix.load(&file_path_demand_matrix[0]);

    matrix time_matrix(n_nodes);
    time_matrix.load(&file_path_time_matrix[0]);

    demand_matrix.print();
    time_matrix.print();

}