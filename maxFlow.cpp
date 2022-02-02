#include <bits/stdc++.h>

using namespace std;

struct Graph {
    int numberNodes, src, dst;
    // edge is pair (neighbor, color)
    // color == 0 ==> BLACK (graph edge)
    // color == 1 ==> RED   (solution edge)
    vector<vector<long long>> cap;
    vector<vector<int>> adj;

    Graph(int n) : numberNodes(n), adj(n), cap(n, vector<long long>(n, 0)) {}

    void AddEdge(int a, int b, int c) {
        adj[a].push_back(b);
        adj[b].push_back(a);
        cap[a][b] += c;
        // cap[b][a] == 0
    }

    int bfs() {
        vector<bool> vis(numberNodes, false);
        vector<int> parent(numberNodes, -1);
        vector<int> q;

        auto push = [&](int node, int par) {
            if (vis[node]) return;
            vis[node] = true;
            parent[node] = par;
            q.push_back(node);
        };                              //pushes the node in queue

        push(src, -1);
        for (int i = 0; i < (int)q.size(); ++i) {
            int node = q[i];            //gets the node from queue
            for (auto nei : adj[node]) {
                if (cap[node][nei] > 0)
                    push(nei, node);
            }
        }

        if (parent[dst] == -1)
            return 0;

        long long flow = LONG_MAX;
        for (int node = dst; node != src; node = parent[node])
            flow = min(flow, cap[parent[node]][node]);
        assert(flow > 0);

        for (int node = dst; node != src; node = parent[node]) {
            cap[parent[node]][node] -= flow;
            cap[node][parent[node]] += flow;
        }
        return flow;
    }

    long long FordFulkerson(int src, int dst) {
        this->src = src; this->dst = dst;
        long long max_flow = 0;
        while (true) {
            long long curr_flow = bfs();
            if (curr_flow == 0) break;
            max_flow += curr_flow;
        }
        return max_flow;
    }
};

int main() {
    ifstream cin("graf.in");
    ofstream cout("graf.out");

    int n, m;
    cin >> n >> m;
    Graph G(n);
    for (int i = 0; i < m; ++i) {
        int a, b, c;
        cin >> a >> b >> c;
        G.AddEdge(a, b, c);
    }

    auto maxFlow = G.FordFulkerson(0, n - 1);
    cout << maxFlow << '\n';

    return 0;
}