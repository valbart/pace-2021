#include "local_search.h"

#include <vector>
#include <queue>
#include <map>
#include <set>
#include <cassert>
#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>
#include <random>
#include <cmath>

int64_t low_mem_local_search(Instance &g)
{
	random_device rd;
	mt19937 rng(rd());
	vector<int> order(g.n);
	for (int i = 0; i < g.n; ++i)
		order[i] = i;
	bool improve;

	int i = 0;
	do
	{
		improve = false;
		shuffle(order.begin(), order.end(), rng);
		for (int j : order)
		{
			if (g.greedy_move(j))
				improve = true;
		}
		if constexpr (DEBUG_LS)
			if (i % 10 == 0)
				cerr << "\rLS | cc: " << g.cost() << " i: " << i << "               ";
		++i;
	} while (improve);

	if constexpr (DEBUG_LS)
	{
		cerr << "\rLS | cc: " << g.cost() << " i: " << i << "               ";
		cerr << endl;
	}
	return g.cost();
}

/* BFS Local search: same as the low_mem_local_search,
 * except that when we move a vertex, we try to move its
 * neighbors next.
 * We do this in a BFS way: each vertex is moved at most once per iteration.
 * Maybe try without this requirement.
 *
 * This part is split in a few functions to avoid allocations and deallocations
 */

bool bfs_ls_process_queue(Instance &g, queue<int> &q, vector<bool> &seen)
{
	bool improve = false;
	int v;
	while (!q.empty())
	{
		v = q.front(); q.pop();
		if (g.greedy_move(v))
		{
			improve = true;
			for (int u : g.neighbors(v))
			{
				if (seen[u])
					continue;
				q.push(u);
				seen[u] = true;
			}
		}
	}
	return improve;
}
	

int64_t bfs_ls_no_alloc(Instance &g, mt19937 &rng, vector<int> &order, vector<bool> &seen, queue<int> &q)
{
	bool improve;
	int64_t i = 0;
	do
	{
		improve = false;
		shuffle(order.begin(), order.end(), rng);
		fill(seen.begin(), seen.end(), false);
		for (int j : order)
		{
			q.push(j);
			seen[j] = true;
			if(bfs_ls_process_queue(g, q, seen))
				improve = true;
		}
		if constexpr (DEBUG_LS)
			if (i % 10 == 0)
				cerr << "\rBFS-LS | cc: " << g.cost() << " i: " << i << "               ";
		++i;
	} while (improve);


	if constexpr (DEBUG_LS)
	{
		cerr << "\rBFS-LS | cc: " << g.cost() << " i: " << i << "               ";
		cerr << endl;
	}
	return g.cost();

}

int64_t bfs_low_mem_local_search(Instance &g)
{
	int n = g.n;
	random_device rd;
	mt19937 rng(rd());
	vector<int> order(n);
	for (int i = 0; i < n; ++i)
		order[i] = i;

	// Bfs variables
	vector<bool> seen(n, 0);
	queue<int> q;

	return bfs_ls_no_alloc(g, rng, order, seen, q);	
}

int64_t restart_low_mem_local_search(Instance &g, int it, vector<int> &best_sol)
{
	int n = g.n;
	random_device rd;
	mt19937 rng(rd());

	// Bfs variables
	vector<int> order(n);
	for (int i = 0; i < n; ++i)
		order[i] = i;
	vector<bool> seen(n, 0);
	queue<int> q;

	int64_t mi = g.cost(), score;
	best_sol = g.sol();
	for (int i = 0; i < it; ++i)
	{
		g.reinit_all_zero();
		score = bfs_ls_no_alloc(g, rng, order, seen, q);
		if (score < mi)
		{
			mi = score;
			best_sol = g.sol();
		}
	}

	return mi;
}

/* When a solution is found,
 * move each vertex to a random cluster with proba p, 
 * then run local search again.
 */
int64_t random_break_low_mem_local_search(Instance &g, int64_t it, double p, vector<int> &best_sol)
{
	int n = g.n;
	random_device rd;
	mt19937 rng(rd());
	uniform_int_distribution unif(0, n - 1);
	uniform_real_distribution zero_one_unif(0.0, 1.0);

	// Bfs variables
	vector<int> order(n);
	for (int i = 0; i < n; ++i)
		order[i] = i;
	vector<bool> seen(n, 0);
	queue<int> q;

	int64_t best_cost = bfs_ls_no_alloc(g, rng, order, seen, q), score, delta;
	best_sol = g.sol();
	int nc;

	for (int64_t i = 0; i < it; ++i)
	{
		for (int v = 0; v < n; ++v)
		{
			if (zero_one_unif(rng) < p)
			{
				nc = unif(rng);
				delta = g.delta_cost(v, nc);
				g.move_with_delta(v, nc, delta);
			}
		}


		score = bfs_ls_no_alloc(g, rng, order, seen, q);
		if (score < best_cost)
		{
			best_cost = score;
			best_sol = g.sol();
		}
		if constexpr (DEBUG_LS)
			if (i%10 == 0)
				cerr << "\rRBLS | cc: " << best_cost << " i: " << i << "            ";
	}
	if constexpr (DEBUG_LS)
		cerr << "\rRBLS | cc: " << best_cost << " i: " << it << "            " << endl;

	return best_cost;
}

int64_t destroy_repair_local_search(Instance &g, int64_t it, int ndestroy, vector<int> &best_sol)
{
	if constexpr (DEBUG_LS)
		cerr << "Starting DRLS with nd: " << ndestroy << endl;

	int n = g.n;
	random_device rd;
	mt19937 rng(rd());
	uniform_int_distribution unif(0, n - 1);

	// Bfs variables
	vector<int> order(n);
	for (int i = 0; i < n; ++i)
		order[i] = i;
	vector<bool> seen(n, 0);
	queue<int> q;

	int64_t best_cost = g.cost(), score;
	best_sol = g.sol();
	vector<int> vs(ndestroy);

	for (int64_t i = 0; i < it; ++i)
	{
		for (int j = 0; j < ndestroy; ++j)
			vs[j] = unif(rng);
		
		g.move_to_zero_size(vs);

		// Process these vertices first
		for (int v: vs)
			q.push(v);

		score = bfs_ls_no_alloc(g, rng, order, seen, q);
		if (score < best_cost)
		{
			best_cost = score;
			best_sol = g.sol();
		}
		g.reinit_state(best_sol, best_cost);
		if constexpr (DEBUG_LS)
			cerr << "\rLNS | cc: " << best_cost << " i: " << i << "            ";
	}
	if constexpr (DEBUG_LS)
	{
		cerr << "\rLNS | cc: " << best_cost << " i: " << it << "            ";
		cerr << endl;
	}
	return g.cost();
}


int64_t move_neighbors_local_search(Instance &g)
{
	int n = g.n;
	random_device rd;
	mt19937 rng(rd());
	vector<int> order(n);
	for (int i = 0; i < n; ++i)
		order[i] = i;

	vector<bool> seen(n);

	bool improve;
	int64_t i = 0;
	vector<int> to_move;
	do
	{
		improve = false;
		shuffle(order.begin(), order.end(), rng);
		fill(seen.begin(), seen.end(), false);
		for (int j : order)
		{
			if (seen[j])
				continue;

			to_move.clear();
			to_move.push_back(j);
			for (int v: g.neighbors(j))
				if (!seen[v])
					to_move.push_back(v);
			if (g.greedy_move_many(to_move))
				improve = true;
			for (int v: to_move)
				seen[v] = true;
		}
		if constexpr (DEBUG_LS)
			if (i % 10 == 0)
				cerr << "\rMN-LS | cc: " << g.cost() << " i: " << i << "               ";
		++i;
	} while (improve);


	if constexpr (DEBUG_LS)
	{
		cerr << "\rMN-LS | cc: " << g.cost() << " i: " << i << "               ";
		cerr << endl;
	}
	return g.cost();
}


int64_t move_neighbors_local_search_same_c(Instance &g)
{
	int n = g.n;
	random_device rd;
	mt19937 rng(rd());
	vector<int> order(n);
	for (int i = 0; i < n; ++i)
		order[i] = i;

	vector<bool> seen(n);

	bool improve;
	int64_t i = 0;
	vector<int> to_move;
	do
	{
		improve = false;
		shuffle(order.begin(), order.end(), rng);
		fill(seen.begin(), seen.end(), false);
		for (int j : order)
		{
			if (seen[j])
				continue;

			to_move.clear();
			to_move.push_back(j);
			for (int v: g.neighbors(j))
				if (!seen[v] && g.sol()[j] == g.sol()[v])
					to_move.push_back(v);
			if (g.greedy_move_many(to_move))
				improve = true;
			for (int v: to_move)
				seen[v] = true;
		}
		if constexpr (DEBUG_LS)
			if (i % 10 == 0)
				cerr << "\rMNSC-LS | cc: " << g.cost() << " i: " << i << "               ";
		++i;
	} while (improve);


	if constexpr (DEBUG_LS)
	{
		cerr << "\rMNSC-LS | cc: " << g.cost() << " i: " << i << "               ";
		cerr << endl;
	}
	return g.cost();
}