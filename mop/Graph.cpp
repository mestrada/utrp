#include "Graph.h"
#include <iterator>     // std::next

using namespace std;

// template <typename ForwardIt>
// ForwardIt next(ForwardIt it, 
//                typename std::iterator_traits<ForwardIt>::difference_type n = 1)
// {
//     std::advance(it, n);
//     return it;
// }


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
        distance[i] = MAXDIST;
    }
    distance[source] = 0;
}

int Graph::getClosestUnmarkedNode()
{
    int minDistance = MAXDIST;
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
    int minDistance = MAXDIST;
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
        //cout<<node+1<<"..";
        insert_node(parent, node);
    }
    else if(predecessor[node] == -1){
        //cout<<"No path from "<<source<<"to "<<node<<endl;
        //TODO: Clean Route
        }
        else {
            printPath(predecessor[node], parent);
            //cout<<node+1<<"..";
            insert_node(parent, node);
        }
}

void Graph::SetPaths()
{
    for(int i=0;i<numOfVertices;i++) {
        if(i == source){}
            //cout<<source+1<<".."<<source+1;
        else
            printPath(i, i);

        //cout<<"->"<<distance[i]<<endl;
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

int Graph::fill_set(int **sol, int l_ind, int max_set_routes, int max_routes){
    list<int>::iterator it; 
    int last = l_ind;
    int prev_node;
    for(int i=0; i<numOfVertices; i++){
        if(short_routes[i].size() < max_routes)
            continue;
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


void Graph::fill_matrix(int** &mat, int &source){
    for(int i=0; i<numOfVertices; i++){

        int cost = 0;
        int index = 0;

        for(list<int>::iterator it=short_routes[i].begin(); it!=short_routes[i].end(); it++){

            if(index != short_routes[i].size() - 1){
            cost += adjMatrix[*it][*(next(short_routes[i].begin(), index + 1))] ;
            }
            index++;
        }

        mat[source][i] = cost;
        mat[i][source] = cost;
    }
}
