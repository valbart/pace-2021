#pragma once
#include "instance.h"

#include <vector>
#include <queue>
#include <random>

constexpr bool DEBUG_LS = false;

int64_t low_mem_local_search(Instance &g);

/* BFS Local search: same as the low_mem_local_search,
 * except that when we move a vertex, we try to move its
 * neighbors next.
 * We do this in a BFS way: each vertex is moved at most once per iteration.
 * Maybe try without this requirement.
 *
 * This part is split in a few functions to avoid allocations and deallocations
 */
bool bfs_ls_process_queue(Instance &g, queue<int> &q, vector<bool> &seen);
int64_t bfs_ls_no_alloc(Instance &g, mt19937 &rng, vector<int> &order, vector<bool> &seen, queue<int> &q);
int64_t bfs_low_mem_local_search(Instance &g);
int64_t restart_low_mem_local_search(Instance &g, int it, vector<int> &best_sol);

/* When a solution is found,
 * move each vertex to a random cluster with proba p, 
 * then run local search again.
 */
int64_t random_break_low_mem_local_search(Instance &g, int64_t it, double p, vector<int> &best_sol);
int64_t destroy_repair_local_search(Instance &g, int64_t it, int ndestroy, vector<int> &best_sol);

int64_t move_neighbors_local_search(Instance &g);
int64_t move_neighbors_local_search_same_c(Instance &g);
