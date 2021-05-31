#include "graph.hpp"
#include <map>
#include <fstream>
#include <iostream>
#include <sstream>
#include <chrono>
#include <random>

Graph::Graph() {
    edges = 0;
    adjacency_matrix = {{}};
    adjacency_list = {{}};
    labels = {};
}

Graph::Graph(std::vector<Node> nodes, std::vector<Edge> edges) {
    this->vertices = nodes.size();
    this->edges = edges.size();
    this->labels = nodes;
    std::vector<std::vector<bool>> adjacency_matrix = std::vector<std::vector<bool>>(vertices);
    std::vector<std::vector<Node>> adjacency_list = std::vector<std::vector<Node>>(vertices);

    std::map<Node, Node> nodes_map;
    for (int i = 0; i < vertices; i++) {
        adjacency_matrix[i] = std::vector<bool>(vertices, false);
        nodes_map[nodes[i]] = i;
    }
    
    for (int i = 0; i < this->edges; i++) {
        Node u = nodes_map.find(edges[i].first)->second;
        Node v = nodes_map.find(edges[i].second)->second;
        adjacency_matrix[u][v] = true;
        adjacency_matrix[v][u] = true;
        adjacency_list[u].push_back(v);
        adjacency_list[v].push_back(u);
    }
    this->adjacency_matrix = adjacency_matrix;
    this->adjacency_list = adjacency_list;
}

Graph::Graph(const Graph& old_graph) {
    this->adjacency_matrix = old_graph.adjacency_matrix;
    this->adjacency_list = old_graph.adjacency_list;
    this->labels = old_graph.labels;
    this->edges = old_graph.edges;
    this->vertices = old_graph.vertices;
}

int Graph::nb_edges() const{
  return this->edges;
}


void Graph::add_edge(Node u, Node v) {
    if (has_edge(u,v))
        return;

    this->adjacency_matrix[u][v] = true;
    this->adjacency_matrix[v][u] = true;
    this->adjacency_list[u].push_back(v);
    this->adjacency_list[v].push_back(u);
    this->edges++;
}

void Graph::remove_edge(Node u, Node v) {
    if (!has_edge(u, v))
        return;

    this->adjacency_matrix[u][v] = false;
    this->adjacency_matrix[v][u] = false;

    auto pos_u = std::find(this->adjacency_list[v].begin(), this->adjacency_list[v].end(), u);
    if (pos_u != this->adjacency_list[v].end()) {
        this->adjacency_list[v].erase(pos_u);
    }
    
    auto pos_v = std::find(this->adjacency_list[u].begin(), this->adjacency_list[u].end(), v);
    if (pos_v != this->adjacency_list[u].end()) {
        this->adjacency_list[u].erase(pos_v);
    }

    this->edges--;
}

/*bool Graph::has_edge(Node u, Node v) const {
    return this->adjacency_matrix[u][v];
}*/

const std::vector<Node> Graph::neighbours(Node u) const {
    /*
    std::vector<Node> adjacency_list(0);
    for (Node v = 0; v < adjacency_matrix.size(); v++)
        if (this->adjacency_matrix[u][v]) adjacency_list.push_back(v);
    return adjacency_list;
    */
   return this->adjacency_list[u];
}

/**
 *  Returns a vector containing the substrsings of the input string s delimited by the charcter delim.
 *  Only used for reading instance files.
 */
std::vector<std::string> split_string (const std::string &s, char delim) {
    std::vector<std::string> result;
    std::stringstream ss (s);
    std::string item;

    while (getline (ss, item, delim)) {
        result.push_back (item);
    }

    return result;
}

Graph* Graph::load_cin(std::istream &is)
{
	int n, m, u, v;
	std::string s;
	is >> s >> s; // ignore header
	is >> n >> m;

	std::vector<Edge> edges;
    std::vector<Node> nodes(n);
	edges.reserve(m);
	for (int i = 0; i < m; ++i)
	{
		is >> u >> v;
		edges.emplace_back(u, v); // vertices are 0-indexed
	}
    for (int i = 0; i < n; i++) nodes[i] = i+1;

	return new Graph(nodes, edges);
}

Graph* Graph::load(std::string fileName) {
    std::string line;
    std::ifstream myfile (fileName);
    std::vector<Node> nodes;
    std::vector<Edge> edges;
    if (myfile.is_open()) {
        getline(myfile, line);
        std::vector<std::string> splitted_line = split_string(line, ' ');
        int nr_vertices = std::stoi(splitted_line[2]);
        int nr_edges = std::stoi(splitted_line[3]);
        nodes.resize(nr_vertices);
        edges.resize(nr_edges);
        // Do not handle comments in instance file. Should it?
        int edge_index = 0;
        while (getline (myfile,line)) {
                std::vector<std::string> splitted_line = split_string(line, ' ');
                Node u = stoi(splitted_line[0]);
                Node v = stoi(splitted_line[1]);
                edges[edge_index] = std::pair<Node,Node> (u,v);
                edge_index++;
            }
        for (int i = 0; i < nr_vertices; i++) nodes[i] = i+1; // Node labels start at 1. 
        myfile.close();
    }
    else std::cout << "Unable to open file" << std::endl;
    return new Graph(nodes, edges); // Applies an unnecessary mapping for nodes name. Implement specific constructor?
}

void Graph::display_adjacency() const {
    for (Node u = 0; u < this->adjacency_matrix.size(); u++) 
    {
        std::cout << " - ";
        for (Node v: this->neighbours(u)) std::cout << v << " ";
        std::cout << std::endl;
    }
    std::cout.flush();
}

int Graph::nr_vertices() const {
    return this->vertices;
}


void Graph::shuffle()
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine rng(seed);
    for (int i = 0; i < vertices; i++)
        std::shuffle(adjacency_list[i].begin(), adjacency_list[i].end(), rng);
}


std::vector<Edge> Graph::all_edges() const {
    std::vector<Edge> edges(0);
    for(Node u = 0; u < nr_vertices(); u++)
        for (Node v: neighbours(u))
            if(v > u)
                edges.push_back(Edge(u, v));
    return edges;
}

Star::Star(Node center, std::vector<Node> leaves):center(center), leaves(leaves) {}

Star::Star() {
    this->center = -1;
    //this->leaves = std::vector<Node> (0);
}

Node Star::get_center() const {
    return this->center;
}

const std::vector<Node> Star::get_leaves() const {
    return this->leaves;
}

bool Star::is_empty() const {
    return this->center == -1;
}

std::vector<std::pair<Node,Node>> Star::all_pairs() const {
    int nr_l = this->leaves.size();
    std::vector<std::pair<Node,Node>> pairs (0); // TODO pre-allocate since we know the size, but have to be carefull with index..
    for (int i = 0; i < nr_l; i++) {
        pairs.push_back(std::pair<Node,Node>(this->center, this->leaves[i]));
        for (int j = i+1; j < nr_l; j++) pairs.push_back(std::pair<Node,Node>(this->leaves[i], this->leaves[j]));
    }
    return pairs;
}

// == requires that leaves are sorted
bool Star::operator==(const Star &rhs) const {
    if (this->center != rhs.get_center() || this->leaves.size() != rhs.get_leaves().size())
        return false;
    std::vector<Node> rhs_leaves = rhs.get_leaves();
    for (unsigned int i = 0; i < this->leaves.size(); i++) {
        if (this->leaves[i] != rhs_leaves[i])
            return false;
    }
    return true;
}

bool Star::has_leaf(Node u) const {
    return std::find(leaves.begin(), leaves.end(), u) != leaves.end();
}

std::vector<Edge> Star::edges() const {
    std::vector<Edge> all_edges(0);
    for (Node u: leaves) all_edges.push_back(Edge(center, u));
    return all_edges;
} 

int Star::nr_leaves() const {
    return leaves.size();
}

std::vector<std::pair<Node, Node>> Star::all_leaves_pairs() const {
    std::vector<std::pair<Node, Node>> leaves_pairs;
    for (unsigned int i = 0; i < leaves.size(); i++)
        for (unsigned int j = i+1 ; j < leaves.size(); j++)
            leaves_pairs.push_back(std::pair<Node, Node>(leaves[i], leaves[j]));
    return leaves_pairs;
}

void Star::print() const
{
    std::cout << "C : " << center << " L : ";
    for (auto l: leaves) std::cout << l << ", ";
    std::cout << std::endl;
}

Star::Star(const Star& other)
{
    center = other.center;
    leaves = other.leaves;
}

Star& Star::operator=(const Star& s) {  
    center = s.center;
    leaves = s.leaves;
    return *this;
}
  
