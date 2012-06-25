#include "Graph.h"

using namespace std;

Graph::Graph(){}

Graph::Graph(int n):numOfVertices(n){

    short_routes = new list<int>[numOfVertices];
    adjMatrix = new int*[numOfVertices];
        for(int i=0; i<numOfVertices; i++){
            adjMatrix[i] = new int[numOfVertices];
        }

    predecessor = new int[numOfVertices];
    distance = new int[numOfVertices];

    mark = new bool[numOfVertices];

}

Graph::~Graph(){

    for(int i = 0; i<numOfVertices; i++){
        delete [] adjMatrix[i];
    }

    delete [] adjMatrix;
    delete [] predecessor;
    delete [] distance;
    delete [] mark;
    delete [] short_routes;

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

void Graph::printPath(int node, int parent)
{
    if(node == source){
        cout<<node+1<<"..";
        insert_node(parent, node);
    }
    else if(predecessor[node] == -1)
        cout<<"No path from "<<source<<"to "<<node<<endl;
    //TODO: Clean Route
        else {
            printPath(predecessor[node], parent);
            cout<<node+1<<"..";
            insert_node(parent, node);
        }
}

void Graph::output()
{
    for(int i=0;i<numOfVertices;i++) {
        if(i == source)
            cout<<source+1<<".."<<source+1;
        else
            printPath(i, i);

        cout<<"->"<<distance[i]<<endl;
    }
}

void Graph::insert_node(int index, int node){
    short_routes[index].push_front(node);
    //cout << "Index: " << index << " node: " << node << endl;
}


void Graph::print_deb(){
    for(int i=0; i<numOfVertices; i++){
        cout << "predecessor " << i << ": " << predecessor[i] << endl;
        cout << "distance " << i << ": " << distance[i] << endl;
        cout << "mark " << i << ": " << mark[i] << endl;
    }
}

int Graph::fill_set(int **sol, int l_ind, int max_set_routes){
    list<int>::iterator it; 
    int last = l_ind;
    int prev_node;
    for(int i=0; i<numOfVertices; i++){
        if(last >= max_set_routes)
            break;
        prev_node = 0;
        for(it = short_routes[i].begin(); it !=short_routes[i].end(); it++){
            sol[*it][last] = prev_node;
            prev_node = *it;
        }
        //cout << "LAST " << last << endl;
        last++;
    }
    return last;
}

void Graph::print_sol(){
    list<int>::iterator it; 
    cout << "SOL " << endl;
    for(int i=0; i<numOfVertices; i++){
        cout << "Sol " << i << endl;
        for(it = short_routes[i].begin(); it !=short_routes[i].end(); it++){
            cout << " " << *it << endl;
        }
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
