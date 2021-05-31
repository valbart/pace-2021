#pragma once

#include "low_mem_graph.h"
#include "instance.h"

#include <vector>
#include <functional>
#include <string>

using namespace std;
using Node = int; 
using Edge = pair<int, int>;
using Cluster = vector<Node>;
using Labeling = vector<int>;
using ConnectedComponent = vector<Node>;

constexpr bool DEBUG_MULTI = true;

// Handle an instance divided in multiple CCs
// Output a solution in O(m) instead of O(n²)
class KernelizedMultiCCInstance
{
private:
	// vector<LowMemGraph> cc_graphs;
	vector<Edge> initial_edges;
	vector<Instance> cc_instances;
	vector<ConnectedComponent> ccs; // ccs[i]: list of the vertices in the i-th connected component.
	vector<pair<int,int>> vertex_to_cc;
	vector<Labeling> solutions;		// labeling of each CC
	vector<int64_t> costs;
	int64_t total_cost;
public:
	KernelizedMultiCCInstance(int n, const vector<Edge> &edges);

	inline int64_t cost() 	const { return total_cost; }
	inline int64_t m() 		const { return initial_edges.size(); }
	inline int64_t _n_cc() 	const { return cc_instances.size(); }
	inline const vector<ConnectedComponent> &_ccs() const { return ccs; }
	inline vector<Instance> &i_ccs() { return cc_instances; }



	void solve(function<bool(int)> exit_condition = [](int i){ (void) i; return false; });
	void sa_multi_cc(
	function<bool(int)> exit_condition,
	function<double(double, int64_t, int64_t)> temp_iter,
	int initial_it,
	int ndestroy, double Tinit,
	int small_cc_threshold);
	
	void print_sol();
	vector<Cluster> get_sol();
	int64_t count_sol();

	static KernelizedMultiCCInstance from_istream(istream &is);
	static KernelizedMultiCCInstance from_cin();
	static KernelizedMultiCCInstance from_file(const string &fname);

};
