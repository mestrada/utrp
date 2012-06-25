#include "Graph.h"

using namespace std;

Graph::Graph(){}

Graph::Graph(int n):numOfVertices(n){

    adjMatrix = new int*[numOfVertices];
        for(int i=0; i<numOfVertices; i++){
            adjMatrix[i] = new int[numOfVertices];
        }

    predecessor = new int[numOfVertices];
    distance = new int[numOfVertices];

    mark = new bool[numOfVertices];

}

void Graph::~Graph(){

    for(int i = 0; i<numOfVertices; i++){
        delete adjMatrix[i];
    }

    delete  adjMatrix;
    delete predecessor;
    delete  distance;
    delete mark;

}

//Read the adjacency matrix for the graph with +ve weights
void Graph::read(int **m_tt){

    for(int i=0;i<numOfVertices;i++)
            for(int j=0;j<numOfVertices;j++)
                adjMatrix[i][j] = m_tt[i][j];
}


void Graph::set_source(int s){
    source = s;
}


void Graph::initialize()
{
    for(int i=0;i<numOfVertices;i++) {
        mark[i] = false;
        predecessor[i] = -1;
        distance[i] = INFINITY;
    }
    distance[source] = 0;
}

int Graph::getClosestUnmarkedNode()
{
    int minDistance = INFINITY;
    int closestUnmarkedNode;
    for(int i=0;i<numOfVertices;i++) {
        if((!mark[i]) && ( minDistance >= distance[i])) {
            minDistance = distance[i];
            closestUnmarkedNode = i;
        }
    }
return closestUnmarkedNode;
}

void Graph::dijkstra()
{
    initialize();
    int minDistance = INFINITY;
    int closestUnmarkedNode;
    int count = 0;
    while(count < numOfVertices) {
        closestUnmarkedNode = getClosestUnmarkedNode();
        mark[closestUnmarkedNode] = true;
            for(int i=0;i<numOfVertices;i++) {
                if((!mark[i]) && (adjMatrix[closestUnmarkedNode][i]>0) ) {
                    if(distance[i] > distance[closestUnmarkedNode]+adjMatrix[closestUnmarkedNode][i]) {
                        distance[i] = distance[closestUnmarkedNode]+adjMatrix[closestUnmarkedNode][i];
                        predecessor[i] = closestUnmarkedNode;
                    }
                }
            }
        count++;
    }
}

void Graph::printPath(int node)
{
    if(node == source)
        cout<<node+1<<"..";
    else if(predecessor[node] == -1)
        cout<<"No path from "<<source<<"to "<<node<<endl;
        else {
            printPath(predecessor[node]);
            cout<<node+1<<"..";
        }
}

void Graph::output()
{
    for(int i=0;i<numOfVertices;i++) {
        if(i == source)
            cout<<source+1<<".."<<source;
        else
            printPath(i);

        cout<<"->"<<distance[i]<<endl;
    }
}

/*
int main()
{
Graph G;
G.read(m_tt);
G.set_source(1);
G.dijkstra();
G.output();
return 0;
}
*/
/*Sample output:

[vinod@thelegendbox cpp]$ g++ dijkstra.cpp
[vinod@thelegendbox cpp]$ ./a.out
Enter the number of vertices of the graph(should be > 0)
5
Enter the adjacency matrix for the graph
To enter infinity enter 999
Enter the (+ve)weights for the row 0
0
2
4
999
1
Enter the (+ve)weights for the row 1
2
0
7
3
999
Enter the (+ve)weights for the row 2
4
7
0
2
999
Enter the (+ve)weights for the row 3
999
3
2
0
6
Enter the (+ve)weights for the row 4
1
999
999
6
0
Enter the source vertex
0
0..0->0
0..1..->2
0..2..->4
0..1..3..->5
0..4..->1

*/
