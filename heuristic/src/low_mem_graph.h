#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;
using Edge = pair<int,int>;

/* Graph class with kernelization
 */
class LowMemGraph
{
private:
	vector<vector<int>> adjs; // Adj lists
	// vector<Edge> initial_edges;
	int64_t _n;
	int64_t _m;

	void add_edge(int u, int v);
	void sort_adjs();
public:
	// LowMemGraph(int n);
	LowMemGraph() {};
	LowMemGraph(int n, const vector<Edge> &edges);

	inline int64_t n() const { return _n; }
	inline int64_t m() const { return _m; }
	inline const vector<int> &neighbors(int u) const { return adjs[u]; }
	inline int64_t deg(int u) const { return adjs[u].size(); }
	inline bool adjacent(int u, int v) const { return binary_search(adjs[u].begin(), adjs[u].end(), v); }
	// inline bool is_initial(int u, int v) const { return binary_search(initial_edges.begin(), 
																		// initial_edges.end(), 
																		// make_pair(min(u, v), max(u, v))); }
	// inline const vector<Edge> &initial_edges_vec() const { return initial_edges; }

	// Kernelization
	bool remove_excess_degree_one();
	bool disjoint_neighborhoods(int u, int v);
	bool remove_edge_disjoint_neighbors();
	bool remove_C4();
	int neighborhoods_intersection_size(int u, int v);
	void remove_edge(int u, int v);
	bool remove_deg3_triangles();
	int64_t kernelize();

	vector<vector<int>> connected_components();


	// debug functions
	void debug() const;
	bool check_adj_lists() const;
	
	// Does not handle comments
	static LowMemGraph from_istream(istream &is);
	static LowMemGraph from_cin();
	static LowMemGraph from_file(const string &fname);
};