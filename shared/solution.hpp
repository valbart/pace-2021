#ifndef SOLUTION
#define SOLUTION

#include <vector>
#include <iostream>
#include <cassert>
#include "graph.hpp"


using Node = int; 
using Edge = std::pair<int, int>;
using Cluster = std::vector<Node>;

class Solution {

private:
	int64_t cost;
	std::vector<Cluster> clusters;
	std::vector<int> cluster_of;
public:
  Solution(int64_t _cost, int n, const std::vector<Cluster> &_clusters);
  Solution(int64_t _cost, const std::vector<int> &_cluster_of);
  int get_cost();
	/**
	 * Outputs the solution sol for the graph g, according to the format specified by PACE.
	 * Takes as input the list of existing edges, must be sorted
	 */
  void print(Graph* G);

};

Solution NO_SOLUTION();


#endif
