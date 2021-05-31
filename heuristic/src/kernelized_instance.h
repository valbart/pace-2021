#pragma once

#include "low_mem_graph.h"
#include "instance.h"

#include <iostream>
#include <vector>
#include <functional>
#include <string>

using namespace std;
using Node = int; 
using Edge = pair<int, int>;
using Cluster = vector<Node>;
using Labeling = vector<int>;
using ConnectedComponent = vector<Node>;


class KernelizedInstance
{
private:
	vector<Edge> initial_edges;
	Instance instance;
	vector<int> vertex_to_instance;
	vector<int> instance_to_vertex;

public:
	KernelizedInstance(int n, const vector<Edge> &edges);

	int64_t bfs_destroy_repair_sa(function<bool(int)> exit_condition,
				int ndestroy,
				double Tinit,
				double _decay_rate);

	void print_sol();
	int64_t count_sol();

	static KernelizedInstance from_istream(istream &is);
	static KernelizedInstance from_cin();
	static KernelizedInstance from_file(const string &fname);
};