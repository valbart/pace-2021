#pragma once

#include <vector>
#include <map>
#include <set>
#include <cassert>
#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <random>

using namespace std;

// Graph of cluster
// Supports efficient moving of vertices
// and cost recomputation ?
class ClusterGraph
{
private:

	typedef struct Edge
	{
		int u, v;
		Edge(int a, int b): u(a), v(b) { }
	} Edge;

	vector<int> cluster_of;
	vector<set<int>> adjs; // Adj lists
	vector<vector<int>> cluster_adjs; // cluster_adjs[c][u]: number of vertices adjacent to $u$ in $c$
	vector<set<int>> clusters; // Adj lists
	vector<Edge> edges;
	int _cost;
	int _n;
public:
	ClusterGraph(int n) : cluster_of(n), adjs(n), cluster_adjs(n, vector<int>(n, 0)), clusters(n), _cost(0)
	{
		_n = n;
		for (int i = 0; i < n; ++i)
		{
			cluster_of[i] = i;
			clusters[i].insert(i);
		}
	}
	
	inline int cost() const { return _cost; }
	inline int n() const { return _n; }
	inline int m() const { return edges.size(); }

	void debug() const
	{
		for (int i = 0; i < _n; ++i)
		{
			cout << "Cluster " << i << ": "; 
			for (int v: clusters[i])
			{
				cout << v << " ";
			}
			cout <<  "--\t";
			for (int j = 0; j < _n; ++j)
				cout << cluster_adjs[i][j] << " ";

			cout << endl;
		}
	}

	int cluster_size(int c) const
	{
		return clusters[c].size();
	}

	// O(log(n))
	void add_edge(int u, int v)
	{
		assert(u != v);
		edges.emplace_back(u, v);
		adjs[u].insert(v);
		adjs[v].insert(u);
		int cu = cluster_of[u];
		int cv = cluster_of[v];
		cluster_adjs[cu][v]++;
		cluster_adjs[cv][u]++;

		// Update cost: reduce if edge inside a cluster, otherwise increase
		if (cu == cv) --_cost;
		else ++_cost;
	}

	// Cost variation when moving v to c
	// O(1)
	int move_delta_cost(int v, int c) const
	{
		int cv = cluster_of[v];
		if (cv == c) return 0;
		return 
		(cluster_size(c) - 2*cluster_adjs[c][v]) // cost to link to new cluster minus 2*existing edges
		- (cluster_size(cv) - 1 - 2*cluster_adjs[cv][v]); // minus the same for current cluster, but do not count itself !
	}

	// Move v to c, recompute cost
	// O(deg(v))
	void move(int v, int c)
	{
		int cv = cluster_of[v];
		_cost += move_delta_cost(v, c);
		cluster_of[v] = c;
		clusters[cv].erase(v);
		clusters[c].insert(v);
		// need to update adjacency of v
		for (int u: adjs[v])
		{
			--cluster_adjs[cv][u];
			++cluster_adjs[c][u];
		}
	}

	// Merge clusters c and d into c, recompute cost
	// O(n)
	void merge(int c, int d)
	{
		assert(c != d);
		for (auto v: clusters[d])
			_cost += (cluster_size(c) - cluster_adjs[c][v]);

		for (auto v: clusters[d])
		{
			clusters[c].insert(v);
			cluster_of[v] = c;
		}

		clusters[d].clear();

		for (int i = 0; i < _n; ++i)
		{
			cluster_adjs[c][i] += cluster_adjs[d][i];
			cluster_adjs[d][i] = 0;
		}
	}
};


// Does not handle comments
static ClusterGraph load(string fname)
{
	ifstream graph_file;
	graph_file.open(fname);
	int n, m, u, v;
	string s;
	graph_file >> s >> s; // ignore header
	graph_file >> n >> m;

	// cout << n << "," << m << endl;
	ClusterGraph g(n);
	for (int i = 0; i < m; ++i)
	{
		graph_file >> u >> v;
		// cout << u << v << endl;
		g.add_edge(u-1, v-1); // vertices are 0-indexed
	}

	graph_file.close();
	return g;
}



int ClusterLocalSearch(ClusterGraph g)
{
	int best_cost, best_cluster, tmp;
	
    random_device rd;
    mt19937 rng(rd());
	vector<int> order(g.n());
	for (int i = 0; i < g.n(); ++i)
		order[i] = i;

	bool improve;
	do
	{
		improve = false;
		shuffle(order.begin(), order.end(), rng);
		for (int j: order)
		{
			best_cost = g.m();
			best_cluster = 0;
			for (int c = 0; c < g.n(); ++c)
			{
				tmp = g.move_delta_cost(j,c);
				if (tmp < best_cost)
				{
					best_cost = tmp;
					best_cluster = c;
				}
			}
			if (best_cost < 0)
			{
				g.move(j, best_cluster);
				improve = true;
			}

		}
	} while (improve);
	return g.cost();
}

int restart_cluster_local_search(ClusterGraph &g, int it)
{
	int mi = g.m();
	for (int i = 0; i < it; ++i)
	{
		mi = min(mi, ClusterLocalSearch(g));
	}
	return mi;
}
