#include "solution.hpp"

#include <algorithm>
#include <iostream>

using namespace std;

Solution::Solution(int64_t _cost, int n, const vector<Cluster> &_clusters): 
	cost(_cost), 
	clusters(_clusters),
	cluster_of(n)
{
	for (uint c = 0; c < clusters.size(); ++c)
		for (int v: clusters[c])
			cluster_of[v] = c;
}


Solution::Solution(int64_t _cost, const vector<int> &_cluster_of): 
	cost(_cost), 
	clusters(_cluster_of.size()), 
	cluster_of(_cluster_of)
{
	for (uint i = 0; i < _cluster_of.size(); ++i)
		clusters[_cluster_of[i]].push_back(i);
}

/**
 * Outputs the solution sol for the graph g, according to the format specified by PACE.
 * Takes as input the list of existing edges, must be sorted
 */
void Solution::print(Graph* G) {
  for (int i = 0; i < G->nr_vertices(); i++){
    for (int j = i+1; j < G->nr_vertices(); j++){
      if ((cluster_of[i] == cluster_of[j]) != G->has_edge(i,j)){
	cout << G->labels[i] << " " << G->labels[j] << "\n";
      }
    }
  }
  cout.flush();
}
// 	  for (const auto &c: clusters)
// 	{
// 		for (int u : c)
// 		{
// 			for (int v : c)
// 			{
// 				if (u < v)
// 				{
// 					Edge e = {u,v};
// 					if (!binary_search(sorted_edges.begin(), sorted_edges.end(), e))
// 					{
// 						cout << u + 1 << " " << v + 1 << "\n";
// 					}
// 				}
// 			}
// 		}
// 	}

// 	for (auto &[u,v]: sorted_edges)
// 	{
// 		if (cluster_of[u] != cluster_of[v])
// 			cout << u + 1 << " " << v + 1 << "\n";
// 	}
// 	cout.flush();
// }

int Solution::get_cost(){
  return this->cost;
}


Solution NO_SOLUTION(){
  std::vector<Cluster> empty;
  return Solution(std::numeric_limits<int>::max(),0,empty);
}


