#ifndef CLUSTERGRAPH
#define CLUSTERGRAPH

#include <vector>
#include <list>
#include <unordered_map>
#include "../shared/graph.hpp"
#include "../shared/solution.hpp"

using namespace std;

using Node = int ; 

// Edges and non edges are represented as pairs to distinguish between
// the marked ones and the non-marked ones. 
// First element of a pair is marked, second element is non-marked.
using Edge = pair<int, int> ;
using NonEdge = pair<int, int> ;
using ClusterId = int;
using Node = int;
using Cluster = vector<Node> ;


class ClusterGraph {

    private :
        /*
        * Pointer to the corresponding graph
        */
        Graph g;
        /*
         Number of vertices of the instance graph (shortcut)
        */
        int n; 

        /* 
        Associate to each node the id of its cluster
        */       
        vector<ClusterId> node_to_cluster; 

        /*Associate to each cluster Id its (possibly empty) list of of nodes
        */
        vector<Cluster> clusters; 

        /*
        The actual cluster graph indicating which clusters can be merged or not.
        */
        vector<vector<bool>> cluster_adjacency;

        vector<vector<ClusterId>> cluster_adjacency_list;

        /*
        An array that associates to a pair of cluster id the numbers of markes/non-markes edges/non-edges
        between the nodes of the two clusters
        0: marked edges, 1: non-marked edges, 2: marked non-edges, 3: non-marked non-edges;
        */
        vector<vector<vector<int>>> cluster_edges_count;

        /*
        Associates each pair to the star containing it in the current lower bound (if exists)
        */

        vector<vector<Star>> marked_graph;

        /*
        * Cost associated to the current graph
        */
        int cost;

        /*
        Returns a list of affected pair after adding edge u,v
        */
        vector<pair<Node, Node>> add_graph_edge(Node u, Node v);
 

        /*
        Returns a list of affected pair after removing edge u,v
        */
        vector<pair<Node, Node>> remove_graph_edge(Node u, Node v);

        // Star related methods

        // Mark the edge u,v with start s
        void mark(Node u, Node v, const Star& s);

        // Unmark the edge u,v
        void unmark(Node u, Node v);

        // Value of the lower bound given by the marked stars
        int star_value;

        int nb_cluster_edges;

        int upper_bound;

    public :

        vector<int> degres;

        vector<vector<int>> twins; 

        ClusterGraph();
        
        ClusterGraph(Graph g, int cost, int upper_bound);

        ClusterGraph(const ClusterGraph&);

        ClusterGraph(const ClusterGraph &old_cluster_graph, Graph g);

        int get_cost() const;

        vector<Cluster> get_all_clusters() const;

        int degree(Node u);
        
        int twinness(Node u, Node v); 
        
        bool false_twins(Node u, Node v);

        int nr_edges(ClusterId cu, ClusterId cv, bool marked);
    
        int nr_non_edges(ClusterId cu, ClusterId cv, bool marked);
        
        int nr_edges(ClusterId cu, ClusterId cv);
        
        int nr_non_edges(ClusterId cu, ClusterId cv);
  
        pair<ClusterId,ClusterId> bestpair();
        
        /*
        Returns whether the pair can be modified in the current solution
        */
        bool is_fixed(Node u, Node v) const;

        /*
        Returns the id of the cluster containing the node u
        */
        ClusterId cluster_Id_of_node(Node u) const;

        vector<Node> neighboring_clusters(ClusterId i) const;

        /*
        Returns the cluster (list of node) with id i
        */
        const Cluster cluster(ClusterId i) const;

        void add_edge(ClusterId i, ClusterId j); 
	
	    void remove_edge(ClusterId i, ClusterId j); 
	
	    //inline bool has_edge(ClusterId i, ClusterId j) const ; 

        inline bool has_edge(ClusterId i, ClusterId j) const
        {
            return cluster_adjacency[i][j];
        }
        
        const Cluster get_cluster(ClusterId i) const;
        
        bool is_finished() const;
  	    
        Solution getSolution() ;

        /*
        Merges the clusters with id j to cluster with id i.
        The id j is freed (added to the empty_cluster list)
        Returns a list of pairs to remark for the LB.
        Does the remarking iff remark == true;
        */
        vector<pair<Node,Node>> merge(ClusterId i, ClusterId j, bool remark);

        /*
        Split clusters with id i and id j and prevent them to be 
        merged in the futur, by removing the (i,j) edge cluster_adjacency
        Returns a list of pairs to remark for the LB.
        Does the remarking iff remark == true;
        */
        vector<pair<Node, Node>> split(ClusterId i, ClusterId j, bool remark);

        pair<int,int> move_cost(Node u,ClusterId i);
        
        void move(Node u, ClusterId i);
  
        /* 
        For debugging purposes 
        */
        void display_cluster(ClusterId i) const;

        void display_cluster_adjacency(ClusterId cu) const;

        /*
        Star related methods
        */
        bool is_marked(Node u, Node v) const;

        Star get_star_at(Node u, Node v) const;

        vector<Star> all_stars() const;

        /*
        Adds the star s to the list of stars used in the lower bound. If one of the pairs in s is marked, this method will 
		throw an error. This will mark all the non-fixed edges of s.
        */
        void mark_star(const Star &s);

        /*
        Removes the star s from the list of stars used in the lower bound.
        */
        void unmark_star(const Star &s);

        int get_star_value() const;

        void clear_stars();

        Solution get_solution(int c) ;

        bool has_graph_edge(Node u, Node v) const;

        /*
        Unmark the start t = marked_graph[u][v] if exists and returns
        all the pairs of t - {u,v} to try to remark these
        */
        vector<pair<Node, Node>> unmark_to_modify(Node u, Node v);

        /*
        Try to remark the pair u,v
        */
        vector<pair<Node, Node>> try_mark(Node u, Node v, bool force_try, bool remark);     

        vector<pair<Node, Node>> try_mark_merge_star(Node u, Node v);   

        vector<pair<Node, Node>> try_mark_triple(Node u, Node v, bool fixed_only);

        vector<pair<Node, Node>> try_mark_split_star(Node u, Node v);

        bool is_feasible() const;

        void set_upper_bound(int ub);

        void init_lower_bound();

        ClusterGraph& operator=(const ClusterGraph& g);

        bool check_packing();

        const vector<Node> graph_neighbors(Node u) const;

        int nr_graph_vertices() const;

        vector<Edge> all_graph_edges() const;

        void shuffle_graph();

        void try_kernelize_lonely_vertex(Node u);
};


#endif
