#ifndef GRAPH
#define GRAPH

#include <vector>
#include <string>
#include <utility>
#include <list>
#include <algorithm>

using Node = int ; 
using Edge = std::pair<int, int> ;
using Cluster = std::vector<Node> ;

class Graph {
	private:
		
	std::vector<std::vector<bool>> adjacency_matrix;

	std::vector<std::vector<Node>> adjacency_list;

	int vertices;
  
	public : 
  
  	int edges;	
  	
	std::vector<int> labels;

  	// Constructor by nodes and edge list, typically used when 
	// splitting a graph according to its CC's
	Graph(std::vector<Node> nodes, std::vector<Edge> edges); 
	
	// Copy constructor
	Graph(const Graph&);

	Graph();	

	void add_edge(Node u, Node v) ; 
	
	void remove_edge(Node u, Node v) ; 
	
	//bool has_edge(Node u, Node v) const ; 
	inline bool has_edge(Node u, Node v) const 
	{
    	return this->adjacency_matrix[u][v];
	}
	
	const std::vector<Node> neighbours(Node u) const ;
	
	static Graph* load(std::string fileName) ;

	int nr_vertices() const;
  
  	int nb_edges() const;
	
	// Print adjacency list of every nodes. For debuging purposes.
	void display_adjacency() const; 

	std::vector<Edge> all_edges() const;

	void shuffle();

	static Graph* load_cin(std::istream &is);
	
} ;

class Star {
	private:

	Node center;

	std::vector<Node> leaves;

	public : 

	Star(const Star&);
	
	Star(Node center, std::vector<Node> leaves) ; 
	
	Node get_center() const ; 
	
	const std::vector<Node> get_leaves() const ;

	// Default constructor to initialize empty star.
	// A start is empty if and only if its center = -1
	Star() ; 

	bool is_empty() const;

	// Return all the pairs of nodes in the star
	std::vector<std::pair<Node,Node>> all_pairs() const;

	bool has_leaf(Node u) const;

	bool operator==(const Star &rhs) const;
	
	std::vector<Edge> edges() const;

	int nr_leaves() const;

	std::vector<std::pair<Node, Node>> all_leaves_pairs() const;
  
  	Star& operator=(const Star& s);
  
	void print() const;
} ;

#endif
