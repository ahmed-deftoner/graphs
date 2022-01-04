#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <set>
using namespace std;

#define INF 10000

class DSU {
	int* parent;
	int* size;

public:
	DSU(int x) {
		parent = new int[x];
		size = new int[x];
		for (int i = 0; i < x; ++i) {
			parent[i] = -1;
			size[i] = 1;
		}
	}
	int find(int v) {
		if (parent[v]==-1)
			return v;
		return parent[v] = find(parent[v]);
	}
	void Union(int x, int y) {
		int a = find(x);
		int b = find(y);
		if (a != b) {
			if (size[a] < size[b])
				swap(a, b);
			parent[b] = a;
			size[a] += size[b];
		}
	}
};

class Graph {
	int V;
	vector<int>* adj;
	vector<pair<int,int>>* adjW;
	bool directed;

public:
	Graph() {
		directed = false;
		V = 0;
		adj = NULL;
		adjW = NULL;
	}
	Graph(int V) {
		this->V = V;
		adj = new vector<int>[V];
		adjW = new vector<pair<int, int>>[V];
		directed = false;
	}
	void isDirected(bool directed) {
		this->directed = directed;
	}
	void addEdge(int src, int des) {
		if (src > V || des > V) {
			cout << "index out bounds, vertex does not exist\n";
			return;
		}
		if (directed == true)
			adj[src].push_back(des);
		else {
			adj[src].push_back(des);
			adj[des].push_back(src);
		}
	}
	void addEdge(int u, int v, int w) {
		adjW[u].push_back({ w,v });
		adjW[v].push_back({ w,u });
	}
	void BFS(int s) {
		bool* visited = new bool[V];
		for (int i = 0; i < V; ++i)
			visited[i] = false;
		BFSUtil(s, visited);
		for (int i = 0; i < V; ++i)
			if (!visited[i])
				BFSUtil(i, visited);
	}
	void BFSUtil(int s,bool visited[]) {
		queue<int> q;
		q.push(s);
		visited[s] = true;
		while (!q.empty()) {
			int x = q.front();
			cout << x << " ";
			q.pop();
			vector<int>::iterator it;
			for (it = adj[x].begin(); it != adj[x].end(); it++) {
				if (!visited[*it]) {
					q.push(*it);
					visited[*it] = true;
				}
			}
		}
	}
	void DFS(int s) {
		bool* visited = new bool[V];
		for (int i = 0; i < V; ++i)
			visited[i] = false;
		DFSUtil(s, visited);
		for (int i = 0; i < V; ++i) 
			if (!visited[i])
				DFSUtil(i, visited);
	}
	void DFSUtil(int x, bool visited[]) {
		visited[x] = true;
		cout << x << " ";
		vector<int>::iterator it;
		for (it = adj[x].begin(); it != adj[x].end(); it++)
			if (!visited[*it])
				DFSUtil(*it, visited);
	}
	bool UndirectedCyclic(int src, bool visited[], int parent) {
		visited[src] = true;
		vector<int>::iterator it;
		for (it = adj[src].begin(); it != adj[src].end(); it++) {
			if (!visited[*it])
			{
				if (UndirectedCyclic(*it, visited, src))
					return true;
			}
			else if (*it != parent)
				return true;
		}
		return false;
	}
	bool DirectedCyclic(int src, bool visited[], bool stack[]) {
		visited[src] = true;
		stack[src] = true;
		vector<int>::iterator it;
		for (it = adj[src].begin(); it != adj[src].end(); it++) {
			if (!visited[*it])
			{
				if (DirectedCyclic(*it, visited, stack))
					return true;
			}
			if (stack[*it]==true)
				return true;
		}
		stack[src] = false;
		return false;
	}
	bool isCyclic() {
		bool* visited = new bool[V];
		for (int i = 0; i < V; ++i)
			visited[i] = false;
		if (directed == false) {
			for (int i = 0; i < V; ++i)
				if (!visited[i])
					return UndirectedCyclic(i, visited, -1);
		}
		else {
			bool* stack = new bool[V];
			for (int i = 0; i < V; ++i)
				stack[i] = false;
			for (int i = 0; i < V; ++i)
				if (!visited[i])
					return DirectedCyclic(i, visited, stack);
		}
		return false;
	}
	void TopologicalSortUtil(int x, bool visited[],stack<int> &s) {
		visited[x] = true;
		vector<int>::iterator it;
		for (it = adj[x].begin(); it != adj[x].end(); it++) 
			if (!visited[*it])
				TopologicalSortUtil(*it,visited,s);
		s.push(x);
	}
	void TopologicalSort() {
		if (directed == true && !this->isCyclic()) {
			bool* visited = new bool[V];
			for (int i = 0; i < V; ++i)
				visited[i] = false;
			stack<int> s;
			for (int i = 0; i < V; ++i)
				if (!visited[i])
					TopologicalSortUtil(i, visited, s);
			while (s.empty()==false) {
				cout << s.top() << " ";
				s.pop();
			} 
		}
	}
	Graph transpose() {
		Graph gr(V);
		gr.isDirected(true);
		for (int i = 0; i < V; ++i) {
			vector<int>::iterator it;
			for (it = adj[i].begin(); it != adj[i].end(); ++it) {
				gr.addEdge(*it, i);
			}
		}
		return gr;
	}
	void fillOrder(int x, bool visited[], stack<int>& s) {
		visited[x] = true;
		vector<int>::iterator it;
		for (it = adj[x].begin(); it != adj[x].end(); it++)
			if (!visited[*it])
				fillOrder(*it, visited, s);
		s.push(x);
	}
	void StronglyConnectedComponents() {
		stack<int> s;
		bool* visited = new bool[V];
		for (int i = 0; i < V; ++i)
			visited[i] = false;
		for (int i = 0; i < V; ++i)
			if (!visited[i])
				fillOrder(i, visited, s);
		Graph gr = transpose();
		for (int i = 0; i < V; ++i)
			visited[i] = false;
		while (!s.empty()) {
			int x = s.top();
			if(!visited[x])
				gr.DFSUtil(x, visited);		  
			s.pop();
		}	   
	}
	int KruskalMST() {
		vector<vector<int>> adjlist;
		for (int u = 0; u < V; ++u) {
			vector<pair<int, int>>::iterator it;
			for (it = adjW[u].begin(); it != adjW[u].end(); ++it)
				adjlist.push_back({ (*it).first,u,(*it).second });
		}
		sort(adjlist.begin(), adjlist.end());
		DSU s(V);
		int ans = 0;
		for (auto edge : adjlist) {
			int w = edge[0];
			int u = edge[1];
			int v = edge[2];
			if (s.find(u) != s.find(v)) {
				s.Union(u, v);
				ans += w;
			} 
		}
		return ans;
	}	
	void PrimMST() {
		priority_queue<pair<int,int>, vector<pair<int,int>>, greater<pair<int,int>>> pq;
		int src = 0; // Taking vertex 0 as source
		// Create a vector for keys and initialize all
		// keys as infinite (INF)
		vector<int> key(V, INF);

		// To store parent array which in turn store MST
		vector<int> parent(V, -1);

		// To keep track of vertices included in MST
		vector<bool> inMST(V+1, false);

		// Insert source itself in priority queue and initialize
		// its key as 0.
		pq.push({ 0, src });
		key[src] = 0;

		/* Looping till priority queue becomes empty   */
		while (!pq.empty())
		{
			// The first vertex in pair is the minimum key
			// vertex, extract it from priority queue.
			// vertex label is stored in second of pair (it
			// has to be done this way to keep the vertices
			// sorted key (key must be first item
			// in pair)
			int u = pq.top().second;
			pq.pop();

			//Different key values for same vertex may exist in the priority queue.
			//The one with the least key value is always processed first.
			//Therefore, ignore the rest.
			if (inMST[u] == true) {
				continue;
			}

			inMST[u] = true;  // Include vertex in MST

			vector<pair<int, int>>::iterator it;
			for (it=adjW[u].begin();it!=adjW[u].end();++it)
			{
				// Get vertex label and weight of current adjacent
				// of u.
				int v = (*it).second;
				int weight = (*it).first;

				//  If v is not in MST and weight of (u,v) is smaller
				// than current key of v
				if (inMST[v] == false && key[v] > weight)
				{
					// Updating key of v
					key[v] = weight;
					pq.push({ key[v], v });
					parent[v] = u;
				}
			}
		}

		// Print edges of MST using parent array
		for (int i = 1; i < V; ++i)
			cout << i << " -> " << parent[i] << endl;
	} 
	void Dijkstra(int src) {
		set<pair<int, int>> setd;
		vector<int> dist(V, INF);
		setd.insert({ 0,src });
		dist[src] = 0;
		while (!setd.empty()) {
			pair<int, int> temp = *(setd.begin());
			setd.erase(setd.begin());
			int u = temp.second;
			vector<pair<int, int>>::iterator it;
			for (it = adjW[u].begin(); it != adjW[u].end();++it) {
				int w = (*it).first;
				int v = (*it).second;
				if (dist[v] > dist[u] + w) {
					if (dist[v] != INF)
						setd.erase(setd.find({ dist[v], v }));
					dist[v] = dist[u] + w;
					setd.insert({ dist[v],v });
				}
			}	   
		}
		for (auto i : dist)
			cout << i << endl;
	}
	void FloydWarshal() {
		int** dist = new int* [V];
		for (int i = 0; i < V; ++i) {
			dist[i] = new int[V];
			for (int j = 0; j < V; ++j)
				dist[i][j] = INF;
		}
		vector<pair<int, int>>::iterator it;
		for (int i = 0; i < V; ++i) 
			for (it = adjW[i].begin(); it != adjW[i].end(); ++it) 
				dist[i][(*it).second] = (*it).first;
		for (int k = 0; k < V; ++k) {
			for (int i = 0; i < V; ++i) {
				for (int j = 0; j < V; ++j) {
					if (dist[i][j] > dist[i][k] + dist[k][j] && dist[i][k] != INF && dist[k][j] != INF)
						dist[i][j] = dist[i][k] + dist[k][j];
				}
			}
		}
		for (int i = 0; i < V; ++i) {
			for (int j = 0; j < V; ++j) {
				if (dist[i][j] == INF)
					cout << -1 << " ";
				else
					cout << dist[i][j]<<" ";
			}
			cout << endl;
		}
	}
	void BellmanFord(int src) {
		vector<int> dist(V, INF);
		dist[src] = 0;
		vector<pair<int, int>>::iterator it;
		for (int i = 0; i < V-1; ++i) {
			for (int j = 0; j < V; ++j) {
				for (it = adjW[j].begin(); it != adjW[j].end(); ++it) {
					int u = j;
					int w = (*it).first;
					int v = (*it).second;
					if (dist[u] != INF && dist[u] + w < dist[v])
						dist[v] = dist[u] + w;
				}
			}
		}
		for (auto i : dist)
			cout << i << " ";
		cout << endl;
	}
};

namespace Tests {
	void BFSTest() {
		cout << "BFS TEST\n";
		Graph g(7);
		g.addEdge(0, 1);
		g.addEdge(0, 2);
		g.addEdge(1, 3);
		g.addEdge(2, 3);
		g.addEdge(3, 4);
		g.addEdge(3, 5);
		g.addEdge(4, 6);
		g.addEdge(5, 6); 
		g.BFS(0);
		cout << endl;
	}
	void DFSTest() {
		cout << "DFS TEST\n";
		Graph g(7);
		g.addEdge(0, 1);
		g.addEdge(0, 2);
		g.addEdge(1, 3);
		g.addEdge(2, 3);
		g.addEdge(3, 4);
		g.addEdge(3, 5);
		g.addEdge(4, 6);
		g.addEdge(5, 6);
		g.DFS(0);
		cout << endl;
	}
	void CyclicTest() {
		cout << "CYCLE TEST\n";
		Graph g(7);
		g.addEdge(0, 1);
		g.addEdge(0, 2);
		g.addEdge(1, 3);
		g.addEdge(2, 3);
		g.addEdge(3, 4);
		g.addEdge(3, 5);
		g.addEdge(4, 6);
		g.addEdge(5, 6);
		if (g.isCyclic())
			cout << "is cyclic\n";
		else
			cout << "acyclic\n";
	}
	void TopologicalSortTest() {
		cout << "TOPOLOGICAL SORT TEST\n";
		Graph g(7);
		g.isDirected(true);
		g.addEdge(0, 1);
		g.addEdge(0, 2);
		g.addEdge(1, 3);
		g.addEdge(2, 3);
		g.addEdge(3, 4);
		g.addEdge(3, 5);
		g.addEdge(4, 6);
		g.addEdge(5, 6);
		g.TopologicalSort();
		cout << endl;
	}
	void SCCTest() {
		cout << "STRONGLY CONNECTED COMPONENTS TEST\n";
		Graph g(8);
		g.addEdge(0, 1);
		g.addEdge(1, 2);
		g.addEdge(2, 3);
		g.addEdge(2, 4);
		g.addEdge(3, 0);
		g.addEdge(4, 5);
		g.addEdge(5, 6);
		g.addEdge(6, 4);
		g.addEdge(6, 7);
		g.StronglyConnectedComponents();
		cout << endl;
	}
	void KruskalTest() {
		cout << "KRUSKALS MST TEST\n";
		Graph g(4);
		g.addEdge(0, 1, 1);
		g.addEdge(1, 3, 3);
		g.addEdge(3, 2, 4);
		g.addEdge(2, 0, 2);
		g.addEdge(0, 3, 2);
		g.addEdge(1, 2, 2);
		cout << g.KruskalMST()<<endl;
	}
	void PrimsTest() {
		cout << "PRIMS TEST\n";
		Graph g(9);
		g.addEdge(0, 1, 4);
		g.addEdge(0, 7, 8);
		g.addEdge(1, 2, 8);
		g.addEdge(1, 7, 11);
		g.addEdge(2, 3, 7);
		g.addEdge(2, 8, 2);
		g.addEdge(2, 5, 4);
		g.addEdge(3, 4, 9);
		g.addEdge(3, 5, 14);
		g.addEdge(4, 5, 10);
		g.addEdge(5, 6, 2);
		g.addEdge(6, 7, 1);
		g.addEdge(6, 8, 6);
		g.addEdge(7, 8, 7);
		g.PrimMST();
	}
	void DijkstraTest() {
		cout << "DIJKSTRA TEST\n";
		Graph g(9);
		g.addEdge(0, 1, 4);
		g.addEdge(0, 7, 8);
		g.addEdge(1, 2, 8);
		g.addEdge(1, 7, 11);
		g.addEdge(2, 3, 7);
		g.addEdge(2, 8, 2);
		g.addEdge(2, 5, 4);
		g.addEdge(3, 4, 9);
		g.addEdge(3, 5, 14);
		g.addEdge(4, 5, 10);
		g.addEdge(5, 6, 2);
		g.addEdge(6, 7, 1);
		g.addEdge(6, 8, 6);
		g.addEdge(7, 8, 7);
		g.Dijkstra(0);
	}
	void BellmanFordTest() {
		cout << "BELLMAN FORD TEST\n";
		Graph g(9);
		g.addEdge(0, 1, 4);
		g.addEdge(0, 7, 8);
		g.addEdge(1, 2, 8);
		g.addEdge(1, 7, 11);
		g.addEdge(2, 3, 7);
		g.addEdge(2, 8, 2);
		g.addEdge(2, 5, 4);
		g.addEdge(3, 4, 9);
		g.addEdge(3, 5, 14);
		g.addEdge(4, 5, 10);
		g.addEdge(5, 6, 2);
		g.addEdge(6, 7, 1);
		g.addEdge(6, 8, 6);
		g.addEdge(7, 8, 7);
		g.BellmanFord(0);
	}
	void FloydWarshallTest() {
		cout << "FLOYD WARHSALL TEST\n";
		Graph g(9);
		g.addEdge(0, 1, 4);
		g.addEdge(0, 7, 8);
		g.addEdge(1, 2, 8);
		g.addEdge(1, 7, 11);
		g.addEdge(2, 3, 7);
		g.addEdge(2, 8, 2);
		g.addEdge(2, 5, 4);
		g.addEdge(3, 4, 9);
		g.addEdge(3, 5, 14);
		g.addEdge(4, 5, 10);
		g.addEdge(5, 6, 2);
		g.addEdge(6, 7, 1);
		g.addEdge(6, 8, 6);
		g.addEdge(7, 8, 7);
		g.FloydWarshal();
	}
}

int main() {
	Tests::BFSTest();
	Tests::DFSTest();
	Tests::CyclicTest();
	Tests::TopologicalSortTest();
	Tests::SCCTest();
	Tests::PrimsTest();
	Tests::KruskalTest();
	Tests::DijkstraTest();
	Tests::BellmanFordTest();
	Tests::FloydWarshallTest();
	return 0;
}