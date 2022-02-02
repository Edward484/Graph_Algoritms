#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <bits/stdc++.h>

using namespace std;
int numberNodes;

void Print_graph( map<int,vector<int>> graph_list){
    for(auto it = graph_list.begin(); it != graph_list.end(); it++){
        cout << it->first << " => " ;
        for(auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
            cout << *it1<<" ";
        }
        cout << endl;
    }
}

void Print_graph_cost( map<int,vector<pair<int,int>>> graph_list){
    for(auto it = graph_list.begin(); it != graph_list.end(); it++){
        cout << it->first << " => " ;
        for(auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
            cout << it1->first<<"("<< it1->second << ")"<< "; ";
        }
        cout << endl;
    }
}

map<int,vector<int>> Read_graph(string file_name){
    map<int,vector<int>> graph_list;
    ifstream f ("graf.in");
    int nodes, edges;
    f >> nodes;
    f >> edges;
    numberNodes = nodes;
    for(int i = 0; i < edges; i++){
        int s_node, e_node;
        f >> s_node;
        f >> e_node;
        if(graph_list.count(s_node) == 0) {
            vector<int> temp;
            temp.push_back(e_node);
            graph_list.insert(pair<int, vector<int>>(s_node,temp ));
        }
        else{
            auto it = graph_list.find(s_node);
            it->second.push_back(e_node);
        }
    }
    return graph_list;
}

map<int,vector<pair<int,int>>> Read_graph_ordered_cost(string file_name){
    map<int,vector<pair<int,int>>> graph_list;
    ifstream f ("graf.in");
    int nodes, edges;
    f >> nodes;
    f >> edges;
    numberNodes = nodes;

    for(int i = 0; i < edges; i++){
        int s_node, e_node, e_cost;
        f >> s_node;
        f >> e_node;
        f >> e_cost;
        if(graph_list.count(s_node) == 0) {
            vector<pair<int,int>> temp;
            temp.push_back(make_pair(e_node,e_cost));
            graph_list.insert(pair<int, vector<pair<int,int>>>(s_node,temp));
        }
        else {
            auto it = graph_list.find(s_node);
            it->second.push_back(make_pair(e_node, e_cost));
        }
    }
    return graph_list;
}

map<int,vector<int>> insertNode(map<int,vector<int>> graph_list, int s_node, int e_node){
    if(graph_list.count(s_node) == 0) {
        vector<int> temp;
        temp.push_back(e_node);
        graph_list.insert(pair<int, vector<int>>(s_node,temp ));
    }
    else{
        auto it = graph_list.find(s_node);
        it->second.push_back(e_node);
    }
    return graph_list;
}

map<int,vector<pair<int,int>>> insertNodeCost(map<int,vector<pair<int,int>>> graph_list, int s_node, int e_node, int e_cost){
    if(graph_list.count(s_node) == 0) {
        vector<pair<int,int>> temp;
        temp.push_back(make_pair(e_node,e_cost));
        graph_list.insert(pair<int, vector<pair<int,int>>>(s_node,temp));
    }
    else {
        auto it = graph_list.find(s_node);
        it->second.push_back(make_pair(e_node, e_cost));
    }
    return graph_list;
};

map<int,vector<int>> Read_graph_unordered(string file_name){
    map<int,vector<int>> graph_list;
    ifstream f ("graf.in");
    int nodes, edges;
    f >> nodes;
    f >> edges;
    numberNodes = nodes;
    for(int i = 0; i < edges; i++){
        int s_node, e_node;
        f >> s_node;
        f >> e_node;
        graph_list = insertNode(graph_list,s_node,e_node);
        graph_list = insertNode(graph_list,e_node,s_node);
    }
    return graph_list;
}

map<int,vector<pair<int,int>>> Read_graph_unordered_cost(string file_name){
    map<int,vector<pair<int,int>>> graph_list;
    ifstream f ("graf.in");
    int nodes, edges;
    f >> nodes;
    f >> edges;
    numberNodes = nodes;

    for(int i = 0; i < edges; i++){
        int s_node, e_node, e_cost;
        f >> s_node;
        f >> e_node;
        f >> e_cost;
        if(graph_list.count(s_node) == 0) {
            vector<pair<int,int>> temp;
            temp.push_back(make_pair(e_node,e_cost));
            graph_list.insert(pair<int, vector<pair<int,int>>>(s_node,temp));
        }
        else {
            auto it = graph_list.find(s_node);
            it->second.push_back(make_pair(e_node, e_cost));
        }
            if(graph_list.count(e_node) == 0){
            vector<pair<int,int>> temp;
            temp.push_back(make_pair(s_node,e_cost));
            graph_list.insert(pair<int, vector<pair<int,int>>>(e_node,temp));
        }
            else{
                auto it = graph_list.find(e_node);
                it->second.push_back(make_pair(s_node,e_cost));
        }
    }
    return graph_list;
}

vector<int> BFS ( int start,map<int,vector<int>> graph_list, vector<int> previosNode){
    queue<int> queueBfs;
    set<int> visited;
    queueBfs.push(start);
    visited.insert(start);
    //init prev with -1
    for(int i = 0;i< 100;i++){
        previosNode.push_back(-1);
    }

    while(!queueBfs.empty()){
        int current = queueBfs.front();
        try {
            vector<int> currentNodeNeigbours;
            if(graph_list.count(current)) {
                auto x = graph_list.find(current)->second;
                for (auto nodeNeighbour: x) {
                    if (visited.count(nodeNeighbour) == 0) {
                        queueBfs.push(nodeNeighbour);
                        visited.insert(nodeNeighbour);
                        previosNode[nodeNeighbour] = current;
                        cout << nodeNeighbour << " ";
                    }
                }
            }
        }
        catch (exception e){
            throw e;
        }

        queueBfs.pop();
    }
    cout << endl;
    for(int i =0;i< 10;i++){
        cout << previosNode[i]<<" ";
    }
    return previosNode;
}
//for finding the shortest path from one node to another
vector<int> reconstructPath(int start, int endNode, vector<int> prev){
    vector<int> path,e;
    for(int i = endNode; i != -1; i = prev[i]){
        path.push_back(i);
    }
    reverse(path.begin(), path.end());
    if(path[0] == start){
        return path;
    } else
        return e ;
}


int CycleDfsRec(int i, map<int,vector<int>> graph,vector<int> fathers,vector<bool> visits  ) {
    visits[i] = true;
    for (auto j = graph.find(i)->second.begin(); j != graph.find(i)->second.end(); j++) {
        if(visits[*j] == false){
            fathers[*j] = i;
            int retur = CycleDfsRec(*j,graph,fathers,visits);
            if(retur == 0)
                return 0;
        }else{
            if(*j != fathers[i]){
                fathers[*j] = i;
                cout<< "Exista ciclu ";
                int x = fathers[*j];
                while(x != *j){
                    cout << x << " ";
                    x = fathers[x];
                }
                cout << *j;
                return 0;
                break;
            }
        }
    };
}

void CycleDFS(int i,map<int,vector<int>> graph){
    vector<int> fathers ;
    vector<bool> visits;
    for(int i =0;i< graph.size()+1;i++){
        visits.push_back(false);
        fathers.push_back(-1);
    }

    CycleDfsRec(i,graph,fathers,visits);
}

//works by starting from 0 edges, after sorting nodes by price,
// it adds them to the tree only if the 2 nodes of the edge don't belong to the same connected component
map<int,vector<pair<int,int>>> Kruskal(map<int,vector<pair<int,int>>> graph ){
    map<int,vector<pair<int,int>>> res;
    //sort
    vector<tuple<int,int,int>> sortedEdges;  //first vertice, second vertice, cost between them
    for(auto it = graph.begin(); it != graph.end(); it++){
        for(auto it1 = it->second.begin(); it1 != it->second.end(); it1++){
            auto t = make_tuple(it->first, it1->first, it1->second);
            if(std::find(sortedEdges.begin(), sortedEdges.end(),make_tuple(it1->first, it->first, it1->second) ) == sortedEdges.end())
                sortedEdges.push_back(t);
        }
    }

    sort(sortedEdges.begin(), sortedEdges.end(),
         [](tuple<int, int, int> const &t1, tuple<int, int , int> const &t2) {
             return get<2>(t1) < get<2>(t2); // lambda expression
         });

    vector<int> conexComp(numberNodes+1);
    int c = 0;
    for(auto node: conexComp){
        conexComp[c] = c;
        c++;
    }

    int nrMuchiiSel = 0;     //no of edges selected
    for(int i = 0 ; nrMuchiiSel < numberNodes -1;i++) {
        if (conexComp[get<0>(sortedEdges[i])] != conexComp[get<1>(sortedEdges[i])]) {           //it doesn't form cycles by being in the same connected component
            nrMuchiiSel++;

            res = insertNodeCost(res,get<0>(sortedEdges[i]),get<1>(sortedEdges[i]),get<2>(sortedEdges[i]));       //add node to new graph
            res = insertNodeCost(res,get<1>(sortedEdges[i]),get<0>(sortedEdges[i]),get<2>(sortedEdges[i]));

            //unify the connected components
            int min, max;
            if(conexComp[get<0>(sortedEdges[i])] < conexComp[get<1>(sortedEdges[i])]){
                min = conexComp[get<0>(sortedEdges[i])];
                max = conexComp[get<1>(sortedEdges[i])];
            }
            else{
                min = conexComp[get<1>(sortedEdges[i])];
                max = conexComp[get<0>(sortedEdges[i])];
            }
                           //keep the connected comp number as the lowest one.
            for(auto it = conexComp.begin(); it!= conexComp.end();it++){
                if(*it == max)
                    *it = min;
            }
        }
    }

    return res;
}

map<int,vector<pair<int,int>>> PrimLazy(map<int,vector<pair<int,int>>> graph ){
    struct Comparator {
            bool operator()(tuple<int, int, int>& t1, tuple<int, int, int>& t2) {
                return get<0>(t1) > get<0>(t2);
            }
    };

    priority_queue<tuple<int,int,int>, vector<tuple<int,int,int>>, Comparator> pq;
    map<int,vector<pair<int,int>>> res,emptyv;
    vector<bool> visited(numberNodes+1, false);
    int selectedEdgeCount = 0;
    int minCost = 0;

    for(auto nei : graph[0]){
        pq.push(make_tuple(nei.second, 0, nei.first));
    }
    visited[0] = true;

    while(!pq.empty() && selectedEdgeCount < numberNodes ){
        auto edge = pq.top();
        pq.pop();
        if(visited[get<2>(edge)])   continue;
        selectedEdgeCount++;
        minCost += get<0>(edge);
        res = insertNodeCost(res,get<1>(edge),get<2>(edge),get<0>(edge));
        res = insertNodeCost(res,get<2>(edge),get<1>(edge),get<0>(edge));

        for(const auto& nei : graph[get<2>(edge)]){
            pq.push(make_tuple(nei.second, get<2>(edge), nei.first));             //cost, origin, destination
        }
        visited[get<2>(edge)] = true;
    }
    if(selectedEdgeCount != numberNodes -1)
        return emptyv;
    Print_graph_cost(res);
    cout <<"The cost of the MSP is: " <<minCost;
    return res;
}

vector<int> DFS(int node, map<int,vector<int>> &graph, vector<int>&order, vector<int>&vis ){
    vis[node] = 1;
    for(auto nei: graph[node]){
        if(vis[nei] == 1){
            cout << "Eroare: exista ciclu";
            exit(1);
        }
        if(!vis[nei]){
            DFS(nei,graph,order,vis);
        }
    }
    vis[node] = 2;
    order.push_back(node);
    return vis;
}

void DFS_cost(int node, map<int,vector<pair<int,int>>> &graph, vector<int>&order, vector<int>&vis ){
    vis[node] = 1;
    for(auto nei: graph[node]){
        if(vis[nei.first] == 1){
            cout << "Eroare: exista ciclu";
            exit(1);
        }
        if(!vis[nei.first]){
            DFS_cost(nei.first,graph,order,vis);
        }
    }
    vis[node] = 2;
    order.push_back(node);
}

void printConnexComponents(map<int,vector<int>> graph){
    vector<int> visited(numberNodes+1,0);
    vector<int> order;
    int nrComp = 0;
    for(auto node : graph){
        if(visited[node.first] == 0){
            nrComp++;
            DFS(node.first,graph,order,visited);
        }
    }
    cout << "Numarulm de componente conexe ale acestui graf este: "<<nrComp<<endl;

}

vector<int> topSort1(map<int,vector<int>> graph){
    // order va mentine sortarea topologica
    vector<int> order;
    // vis va contoriza daca un nod e vizitat sau nu
    vector<int> vis(graph.size()+1, 0);

    for (int i = 1; i <= graph.size() ; ++i) {
        if (!vis[i]) {
            DFS(i, graph, order, vis);
        }
    }
    reverse(order.begin(), order.end());
    for (auto node : order)
        cout << node << " ";
    cout << '\n';
    return order;
}

vector<int> topSort1Cost(map<int,vector<pair<int,int>>> graph){
    // order va mentine sortarea topologica
    vector<int> order;
    // vis va contoriza daca un nod e vizitat sau nu
    vector<int> vis(graph.size()+1, 0);

    for (int i = 1; i <= graph.size() ; ++i) {
        if (!vis[i]) {
            DFS_cost(i, graph, order, vis);
        }
    }
    reverse(order.begin(), order.end());
    for (auto node : order)
        cout << node << " ";
    cout << '\n';
    return order;
}

vector<int> dagShortestPath(map<int,vector<pair<int,int>>> graph, int start){
    vector<int> topSort = topSort1Cost(graph);
    vector<int> distanceArray(graph.size()+10,9999);
    distanceArray[start] = 0;
    for(auto node : graph){
        int nodeIndex = topSort[node.first];
        vector<pair<int,int>> nodeNeighbours = node.second;
        for(auto neighbour: nodeNeighbours){
            int newDistance = distanceArray[node.first] + neighbour.second;
            distanceArray[neighbour.first] = min(newDistance, distanceArray[neighbour.first]);
        }
    }
     return distanceArray;
}

typedef pair<int, int> pi;

vector<int> DijkstraLazy(map<int,vector<pair<int,int>>> graph,int start){
    set<int> visitedArray;
    vector<int> previousNodeArray(graph.size()+4, -1);
    vector<int> distanceArray(graph.size()+4, 999999);
    distanceArray[start] = 0;
    priority_queue<pi, vector<pi>, greater<pi> > pq; //ordered by first element. So here we have pair<distance,node>
    pq.push(make_pair(0,start));  
    while(pq.size() != 0 ){
        auto possiblyBest = pq.top();                                   //iau cel mai bun nod posibil
        pq.pop();
        visitedArray.insert(possiblyBest.second);                       //il notez ca vizitat
        if(distanceArray[possiblyBest.second] >= possiblyBest.first){    //optimizare in caz ca iau o distanta mai proasta decat ce exista deja
            for(auto nodeNeighbour : graph[possiblyBest.second]){
                if(visitedArray.count(nodeNeighbour.first) == 0){       //daca nu am vizitat nodul
                    int newDist = distanceArray[possiblyBest.second] + nodeNeighbour.second;  //distanta pana la nod + distanta de la nod la vecin
                    if(newDist < distanceArray[nodeNeighbour.first]){
                        previousNodeArray[nodeNeighbour.first] = possiblyBest.second;
                        distanceArray[nodeNeighbour.first] = newDist;
                        pq.push(make_pair(newDist, nodeNeighbour.first ));
                    }
                }
            }
        }
    }
    return distanceArray;
}

vector<int> bellmanFord(map<int,vector<pair<int,int>>> graph, int start, int length = 100){
    vector<int> distanceArray(graph.size()+10, 999999);
    distanceArray[start] = 0;
    for(auto node: graph){
        for(auto nodeNeighbour : node.second){
            if(distanceArray[node.first] + nodeNeighbour.second < distanceArray[nodeNeighbour.first]){
                distanceArray[nodeNeighbour.first] = distanceArray[node.first]+ nodeNeighbour.second;
            }
        }
    }

    for(auto node: graph){
        for(auto nodeNeighbour : node.second){
            if(distanceArray[node.first] + nodeNeighbour.second < distanceArray[nodeNeighbour.first]){
                distanceArray[nodeNeighbour.first] = -99999;
            }
        }
    }

    return distanceArray;
}

int** fromListToMatrix(map<int,vector<pair<int,int>>> graph, int length){
    int **m = 0;
    m = new int*[length];
    for(int i =0;i < length; i++){
        for(int j =0; j < length; j++){
            m[i][j] = numeric_limits<int>::max();
        }
    }
    for(auto node: graph){
        for(auto nodeNeighbour : node.second){
            m[node.first][nodeNeighbour.first] = nodeNeighbour.second;
        }
    }

    return m;
}


vector<int> floydMarshalManager(map<int,vector<pair<int,int>>> graph,int start, int end){

    //move from adjacency list to adjacency matrix
    int length = numberNodes;
    int m[numberNodes+2][numberNodes+2];
    for(int i =0;i < length; i++){
        for(int j =0; j < length; j++){
            if(i == j)
                m[i][j] = 0;
            else
                m[i][j] = 99999;
        }
    }
    for(auto node: graph){
        for(auto nodeNeighbour : node.second){
            m[node.first][nodeNeighbour.first] = nodeNeighbour.second;
        }
    }



    //copy contents from m in dp and next
    int next[numberNodes+2][numberNodes+2];
    int dp[numberNodes+2][numberNodes+2];
    for(int i =0;i < length; i++){
        for(int j =0; j < length; j++){
            if(m[i][j] != 99999 ){
                next[i][j] = j;
            }
            else{
                next[i][j] = 0;
            }
            dp[i][j] = m[i][j];
        }
    }



    //solve
    for(int k =0;k< length;k++) {
        for (int i = 0; i < length; i++) {
            for (int j = 0; j < length; j++) {
                if(dp[i][k] + dp[k][j] < dp[i][j]){
                    dp[i][j] = dp[i][k] + dp[k][j];
                    next[i][j] = next[i][k];
                }
            }
        }
    }


    //identify negative cycles
    for(int k =0;k< length;k++) {
        for (int i = 0; i < length; i++) {
            for (int j = 0; j < length; j++) {
                if(dp[i][k] < 9999 &&  dp[k][j] < 9999 && dp[i][j] < 0){
                    dp[i][j] = -99999;
                    next[i][j] = -99999;
                }
            }
        }
    }

    //reconstruct shortestPath
    vector<int> path;
    vector<int> null;
    if(dp[start][end] > 9999) return path;          //starts from negative cycle
    int current = 0;
    for (current = start; current != end; current = next[current][end]){
        if(current > -9999) return null;            //reaches negative cycle
        path.push_back(current);
    }
    if(next[current][end] > -9999) return null;
    path.push_back(end);                //include last node in path
    return path;                    //good path


}

void dfsBridges(int node,int parent, map<int,vector<int>> &graph,  set<int>&vis,
                vector<pair<int,int>> &bridges, vector<int> &lowLinks, vector<int> &ids, int &id ){
    vis.insert(node);
    id++;
    ids[node] = id;
    lowLinks[node] = ids[node];
    for(auto nei: graph[node]){
        if(nei == parent) continue;
        if(vis.count(nei) == 0 ){
            dfsBridges(nei, node, graph, vis,bridges, lowLinks, ids, id);
            lowLinks[node] = min(lowLinks[node], lowLinks[nei]);
            if(ids[node] < lowLinks[nei]){
                bridges.push_back(make_pair(node,nei));
            }
        }
        else{
            lowLinks[node] = min(lowLinks[node], ids[nei]);
        }
    }
}

vector<pair<int,int>> findBridges(map<int,vector<int>> graph){
    vector<int> lowLinks(numberNodes);
    vector<int> ids(numberNodes,0);
    set<int> visited;
    vector<pair<int,int>> bridges;
    int id = 0;
    for(int i = 0; i < numberNodes; i++){
        if(visited.count(i) == 0){
            dfsBridges(i,-1, graph,visited,bridges,lowLinks, ids, id);
        }
    }
    return bridges;
}

void dfsArt(int root,int node, int parent,int &outEdgeCount, map<int,vector<int>> &graph,
                set<int>&vis, vector<bool> &artPts, vector<int> &lowLinks, vector<int> &ids, int &id ){
    vis.insert(node);
    id++;
    ids[node] = id;
    lowLinks[node] = ids[node];
    if (parent == root) outEdgeCount++;
    for(auto nei: graph[node]){
        if(nei == parent) continue;
        if(vis.count(nei) == 0 ){
            dfsArt(root,nei, node,outEdgeCount, graph, vis,artPts, lowLinks,ids,id);
            lowLinks[node] = min(lowLinks[node], lowLinks[nei]);
            if(ids[node] < lowLinks[nei]){
                artPts[node] = true;
            }
            if(ids[node] == lowLinks[nei]){
                artPts[node] = true;
            }
        }
        else{
            lowLinks[node] = min(lowLinks[node], ids[nei]);
        }
    }
}

vector<bool> findArticulationPoints(map<int,vector<int>> graph){
    vector<int> lowLinks(numberNodes);
    vector<int> ids(numberNodes,0);
    set<int> visited;
    vector<bool> artPts(numberNodes, false);
    int outEdgeCount = 0;
    int id = 0;
    for(int i = 0; i < numberNodes; i++){
        if(visited.count(i) == 0){
            outEdgeCount = 0;
            dfsArt(i,i,-1,outEdgeCount, graph,visited,artPts,lowLinks, ids,id);
            artPts[i] = (outEdgeCount > 1);
        }
    }
    return artPts;
}

void dfsTarjan(int node,map<int, vector<int>> &graph, vector<int> &lowLinks, vector<int>&ids,
               int &id, int &sccCount, stack<int> &stack, vector<bool> &onStack){
    id++;
    lowLinks[node] = id;
    ids[node] = lowLinks[node];
    stack.push(node);
    onStack[node] = true;
    //find strongly connected components
    for(auto nodeNeighbour: graph[node]){
        if(ids[nodeNeighbour] == -1 ){      //if this node is unvisited
            dfsTarjan(nodeNeighbour, graph, lowLinks, ids,id, sccCount, stack, onStack);
        }
        if(onStack[nodeNeighbour] == true){
            lowLinks[node] = min(lowLinks[node], lowLinks[nodeNeighbour]);      //update lowLinks just if the nei is on the stack
        }
    }
    //clear stack of scc if you arrived at a beginning of a scc
    if(ids[node] == lowLinks[node]){
        auto nodeOnStack = stack.top();
        onStack[nodeOnStack] = false;
        stack.pop();
        while(nodeOnStack != node){
            nodeOnStack = stack.top();
            onStack[nodeOnStack] = false;
            stack.pop();
            lowLinks[nodeOnStack] = ids[node];
        }
        sccCount++;
    }

}

//find strongly connected components
vector<int> Tarjan(map<int,vector<int>> graph){
    vector<int> lowLinks(numberNodes);
    vector<int> ids(numberNodes,-1);
    int id = 0;
    int sccCount = 0; //counts the number of strongly connected components
    stack<int> stack;
    vector<bool>onStack(numberNodes, false);

    for(int i =0;i< numberNodes;i++){
        if(ids[i] == -1){
            dfsTarjan(i, graph, lowLinks, ids,id, sccCount, stack, onStack);
        }
    }
    return lowLinks;
}


int main() {
    map<int,vector<pair<int,int>>> graph = Read_graph_unordered_cost("graf.in");
    Print_graph_cost(graph);
    PrimLazy(graph);
}
