/*

Title: Dijkstra’s single source shortest path algorithm
Author: Vinodh Kumar B

Description:
adjMatric[][] : This array is used to store the adjacency matrix of the input graph
predecessor[] : This array is used to store the predecessor elements. predecessor[i] is the previous element in the path from source to current node i.
ex: If 1->2->4->3 is the shortest path from node1(source) to node 3. predecessor[3] is 4.
distance[]      : This array will store the minimum distance from the source to corresponding node. It will be updated in each iteration in dijkstra’s                     algorithm.
mark[]          — This array is used to mark the nodes as visited. mark[i] = true, This means node i is marked as visited.
source          — source is the starting node from which we have to find out the shortest distance to all other nodes in the graph.
numOfVertices — The number of nodes in the graph.
Note: Node ranges will be 0 to numOfVertices-1

getClosestUnmarkedNode()    : This method will return the closest unmarked node from the current marked node.
printPath(i)                : This will print the path from source node to i in recursive manner
output()                    : This method will print paths to all nodes from source and their shortest distance
read()                        : This method will read the graph values from input
initialize()                : This method will set the initial values of predecessor[] to -1 and
dijkstra()                    : This method will implement the dijkstra’s algorithm

*/

#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <list>
#define MAXDIST 999
using namespace std;

class Graph
{
private:
int **adjMatrix;
int *predecessor, *distance;
bool *mark;
int source;
int numOfVertices;
list<int>  *short_routes;
int** cost_matrix;
public:
Graph();
Graph(int n);
~Graph();
void read(int **);
void set_source(int);
void initialize();
int getClosestUnmarkedNode();
void dijkstra();
void SetPaths();
void printPath(int, int);
void print_deb();
void insert_node(int, int);
int fill_set(int **, int, int, int);
void print_sol();
void fill_matrix(int** &, int);
};



#endif