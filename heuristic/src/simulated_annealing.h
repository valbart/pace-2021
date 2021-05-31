#pragma once
#include "instance.h"

#include <vector>
#include <random>

constexpr bool DEBUG_SA = false;

/******************** SIMULATED ANNEALING ***************/
// PARAMS, recommandés par un expert de la NASA
constexpr double w = 0.0001;
constexpr double decay_rate = 0.9999975;

double initial_T(int64_t cost);
bool accept(double T, int delta, double p);

// Simulated annealing that tries to move a random vertex to a random cluster
int64_t sa_single_vertex_move(Instance &g, int64_t it);

// Simulated annealing that tries to move a random vertex to a random cluster
// When a vertex is moved, try to move its neighbors to the same cluster next
int64_t sa_bfs(Instance &g, int64_t it);

// Simulated annealing that tries to move *each* vertex to a random cluster
// When a vertex is moved, try to move its neighbors to the same cluster next
int64_t sa_bfs2(Instance &g, int64_t it);

// Simulated annealing that destroys a random number of vertices (unif(1, ndestroy))
// and then greedily moves them back in the solution
int64_t random_destroy_repair_sa(Instance &g, int64_t it, int ndestroy, double Tinit);

void bfs_destroy_repair_sa_one_it(
	Instance &g, int n, int ndestroy, mt19937 &rng, 
	vector<bool> &seen, double T, int64_t &best_cost);
bool bfs_destroy_repair_sa_one_it_order(
	Instance &g, int n, int ndestroy, mt19937 &rng, 
	vector<int> &order, vector<bool> &seen, double T, vector<int> &best_sol, int64_t &best_cost);
int64_t bfs_destroy_repair_sa(Instance &g, int64_t it, int ndestroy, double Tinit);

// Simulated annealing that greedily moves a vertex and its neighbors
int64_t move_neighbors_sa(Instance &g, int64_t it, double Tinit);

// Simulated annealing that greedily moves a vertex and its neighbors
// that are in the same cluster
int64_t move_neighbors_same_cc_sa(Instance &g, int64_t it, double Tinit);
