#include "lp_solve.hpp"

#include "../shared/graph.hpp"
#include "../shared/solution.hpp"

#include "../coin/ClpSimplex.hpp"
#include "../coin/CoinPackedVector.hpp"
#include "../coin/CoinPackedMatrix.hpp"

#include <unordered_map> 
#include <vector>
#include <set>
#include <functional>
#include <iostream>
#include <cmath>


// A hash function used to hash a pair of any kind
struct hash_pair {
    template <class T1, class T2>
    size_t operator()(const std::pair<T1, T2>& p) const
    {
        auto hash1 = std::hash<T1>{}(p.first);
        auto hash2 = std::hash<T2>{}(p.second);
        return hash1 ^ hash2;
    }
};

class NullHandler : public CoinMessageHandler {
	
	
	public:
	virtual int print() {
		return 0 ; 
	}
};

Solution lp_solve(Graph * g){
	ClpSimplex model ; 
	NullHandler handler ; //Message handler to disable output by LP solver.
	model.passInMessageHandler(&handler) ;
	CoinPackedMatrix matrix(false, 0, 0) ; 
	
	int n = g->labels.size() ;
	if(n > 100) {
		return NO_SOLUTION() ;
	}
	int nb_pairs = (n * (n-1)) / 2 ; 
	std::unordered_map<Edge, int, hash_pair> pair_index ; 
	std::vector<double> upper_bounds ; 
	std::vector<double> lower_bounds ; 
	
	matrix.setDimensions(0, nb_pairs) ;
	int last_index = -1 ; 
	for(int i = 0 ; i < n ; i++) {
		for(int j = i +1 ; j < n ; j++) {
			last_index++ ;
			Edge e(i, j) ;
			Edge einv(j, i) ; 
			pair_index[e] = last_index ; 
			pair_index[einv] = last_index ;
		}
	}
	
	double* variable_lb = new double[nb_pairs]{0} ; 
	double* variable_ub = new double[nb_pairs]{0} ; 
	
	//add all P3 as constraints
	for(int u = 0 ; u < n ; u++) {
		for(const int v : g->neighbours(u)) {
			for(const int w : g->neighbours(u)) {
				if(w > v and not g->has_edge(v, w)) {
					int e1 = pair_index[Edge(u, v)] ; 
					int e2 = pair_index[Edge(v, w)] ; 
					int e3 = pair_index[Edge(u, w)] ; 
					
					CoinPackedVector row ; 
					row.insert(e1, 1) ; 
					row.insert(e2, 1) ; 
					row.insert(e3, 1) ; 
					matrix.appendRow(row) ;
					variable_ub[e1]++ ; 
					variable_ub[e2]++ ; 
					variable_ub[e3]++ ; 
					upper_bounds.push_back(3) ; 
					lower_bounds.push_back(1) ;
				}
			}
		}
	}
	
	//add all claw as constraints
	for(int u = 0 ; u < n ; u++) {
		for(const int v : g->neighbours(u)) {
			for(const int w : g->neighbours(u)) {
				if(w > v and not g->has_edge(v, w)) {
					for(const int z : g->neighbours(u)) {
						if(z > w and not g->has_edge(v, z) and not g->has_edge(w, z)) {
							std::vector<int> leaves{v, w, z} ;
							Star s(u, leaves) ;
							CoinPackedVector row ; 
							for(Edge e : s.all_pairs()) {
								row.insert(pair_index[e], 1) ;
								variable_ub[pair_index[e]]++ ; 
							}
							matrix.appendRow(row) ;
							upper_bounds.push_back(6) ; 
							lower_bounds.push_back(2) ;							
						}
					}
				}
			}
		}
	}
	double * objective = new double[nb_pairs] ; 
	for(int i = 0 ; i < nb_pairs ; i++) {
		objective[i] = 1 ;
	}
	
	model.loadProblem(matrix, variable_lb, variable_ub, objective, lower_bounds.data(), upper_bounds.data()) ;
	//model.writeMps("example.mps");
	model.initialSolve() ;
	
	int lower_bound = std::round(model.getObjValue()) ;
	double* solution = model.primalColumnSolution();
	
	Graph edited_graph(*g) ;
	int nb_edits = 0 ;
	for(int i = 0 ; i < n ; i++) {
		for(int j = i+1 ; j < n ; j++) {
			if(solution[pair_index[Edge(i, j)]] > 0.5) {
				nb_edits ++ ;
				if(edited_graph.has_edge(i, j)) 
					edited_graph.remove_edge(i, j) ;
				else 
					edited_graph.add_edge(i, j) ;
			}
		}
	}
	
	if(nb_edits != lower_bound) 
		return NO_SOLUTION() ;
	
	//check if the edited graph is a union of disjoint clusters
	for(int i = 0 ; i < n ; i++) {
		for(int j : edited_graph.neighbours(i)) {
			for(int k : edited_graph.neighbours(i)) {
				if(j < k and not edited_graph.has_edge(j, k)) {
					return NO_SOLUTION() ;
				}
			}
		}
	}
	
	std::vector<Cluster> clusters ; 
	std::set<int> visited_vertices ;
	for(int i = 0 ; i < n ; i++) {
		if(visited_vertices.find(i) == visited_vertices.end()) {
			Cluster c ; 
			c.push_back(i) ; 
			visited_vertices.insert(i) ;
			for(int j : edited_graph.neighbours(i)){
				c.push_back(j) ; 
				visited_vertices.insert(j) ;
			}
			clusters.push_back(c) ;
		}
	}
	return Solution(nb_edits, n, clusters) ;
}
