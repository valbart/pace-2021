#ifndef EXACTINSTANCE
#define EXACTINSTANCE

#include "../shared/graph.hpp"
#include "../shared/solution.hpp"
#include "cluster_graph.hpp"

using Cluster = std::vector<Node>;


class ExactInstance
{
private:
	//Graph *g;

	ClusterGraph cluster_graph; // Contains the cost of the instance

	int upper_bound;
	
	int nb_branch;

	int nr_graph_vertices;

public:
	//other missing arguments depending on how we want to implement
	//the cluster graph (if necessary?), the stars used for the lower bound.
	ExactInstance(Graph g, int upper_bound, int cost);

	// Copy constructor
	ExactInstance(const ExactInstance &);

	int get_upper_bound();

	int get_size();

	int get_cost();

	int degree(Node u);

	int twinness(Node u, Node v);

	void set_upper_bound(int value);

	void preprocess();
	void init_lower_bound();
	void recompute_lower_bound();
	void forced_moves();

	// Returns true iff there is an edge (u,v) in g
	bool has_edge(Node u, Node v);

	// Return the list of neighbors of u in g
	std::vector<Node> neighbours(Node u);

	// Return the number of vertices in g
	int nr_vertices();

	vector<pair<Node, Node>> merge(ClusterId c1, ClusterId c2, bool remark);

	vector<pair<Node, Node>> split(ClusterId c1, ClusterId c2, bool remark);
	
	void isolate_cluster(ClusterId c);
	bool can_merge(ClusterId c1, ClusterId c2) const;

	bool is_fixed(Node u, Node v) const;

	int lower_bound() const;

	bool is_feasible() const;

	bool is_finished() const;

	const Cluster get_cluster(Node u) const;

	ClusterId get_cluster_Id(Node u) const;

	void fix_pair(Node u, Node v);

	void display_cluster_adjacency(ClusterId cu) const;

	void display_cluster(ClusterId cu) const;

	void display_graph_adjacency() const;

	
	bool is_marked(Node u, Node v) const;

	const Star get_star(Node u, Node v) const;

	const std::vector<Star> all_stars() const;

	void mark_star(Star s);

	void unmark_star(Star s);

	Solution algo_brute();

	Solution heuristic();

	void clear_stars();

	int get_star_value() const;
	
	const std::vector<Star> get_stars_at(Node u) const;

	/*Computes a lower bound by
	* computing a packing of triples by computing an IS in a conflict graph
	*/
	void lower_bound_by_triple();

	void greedy_complete_lower_bound();

	const Cluster cluster(ClusterId i) const;

	bool check_packing();

	~ExactInstance();

	bool has_graph_edge(Node u, Node v) const;

	const vector<Node> graph_neighbors(Node u) const;

	void shuffle_graph();

	bool init_lower_bound_triple();

	bool init_lower_bound_greedy();

	void fast_lower_bound_triple();

	void fast_recompute_lower_bound();

	void fast_init_lower_bound();

	bool triple_to_claw();

};

#endif
