#include "simulated_annealing.h"
#include "local_search.h"
#include "profiler.h"

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

/******************** SIMULATED ANNEALING ***************/
double initial_T(int64_t cost)
{
	return w * cost / log(2);
}

bool accept(double T, int delta, double p)
{
	return p < exp(-delta / T);
}

// Simulated annealing that tries to move a random vertex to a random cluster
int64_t sa_single_vertex_move(Instance &g, int64_t it)
{
	int n = g.n;
	int64_t best_cost = g.cost();
	vector<int> best_sol = g.sol();
	double T = initial_T(best_cost);

	random_device rd;
	mt19937 rng(42);
	uniform_int_distribution unif(0, n - 1);
	uniform_real_distribution zero_one_unif(0.0, 1.0);

	int v, c;
	int64_t delta;
	for (int64_t i = 0; i < it; ++i)
	{
		v = unif(rng);
		c = unif(rng);
		delta = g.delta_cost(v, c);
		if (accept(T, delta, zero_one_unif(rng)))
		{
			g.move_with_delta(v, c, delta);
			if (g.cost() < best_cost)
			{
				best_cost = g.cost();
				best_sol = g.sol();
			}
		}

		if constexpr (DEBUG_SA)
			if (i % 10'000 == 0)
				cerr << "\rSA | bc: " << best_cost << " cc: " << g.cost() << " t: " << T << " i: " << i << "/" << it << "                 ";

		T *= decay_rate;
	}
	if constexpr (DEBUG_SA)
		cerr << endl;

	g.reinit_state(best_sol, best_cost);
	int64_t tmp = low_mem_local_search(g);

	cerr << "sa cost: " << best_cost << " ls cost: " << tmp << endl;

	if (tmp < best_cost)
	{
		best_cost = tmp;
		best_sol = g.sol();
	}
	return best_cost;
}

// Simulated annealing that tries to move a random vertex to a random cluster
// When a vertex is moved, try to move its neighbors to the same cluster next
int64_t sa_bfs(Instance &g, int64_t it)
{
	int n = g.n;
	int64_t best_cost = g.cost();
	vector<int> best_sol = g.sol();
	double T = initial_T(best_cost);

	// Randomness sources
	random_device rd;
	mt19937 rng(rd());
	uniform_int_distribution unif(0, n - 1);
	uniform_real_distribution zero_one_unif(0.0, 1.0);

	// BFS variables
	vector<bool> seen(n, false);
	vector<int> visited;
	queue<int> q;

	int v, c;
	int64_t delta;
	for (int64_t i = 0; i < it; ++i)
	{
		v = unif(rng);
		c = unif(rng);
		q.push(v);
		while (!q.empty())
		{
			v = q.front();
			q.pop();
			if (seen[v])
				continue;
			seen[v] = true;
			visited.push_back(v);
			delta = g.delta_cost(v, c);
			if (accept(T, delta, zero_one_unif(rng)))
			{
				g.move_with_delta(v, c, delta);
				if (g.cost() < best_cost)
				{
					best_cost = g.cost();
					best_sol = g.sol();
				}
				// try to move neighbors next
				for (int u : g.neighbors(v))
				{
					q.push(u);
				}
			}

			if constexpr (DEBUG_SA)
				if (i % 10'000 == 0)
					cerr << "\rSA-BFS | bc: " << best_cost << " cc: " << g.cost() << " t: " << T << " i: " << i << "/" << it << "                 ";
			++i;
			T *= decay_rate;
		}
		// Cheap reset for seen
		for (int u : visited)
			seen[u] = false;
		visited.clear();
	}
	if constexpr (DEBUG_SA)
		cerr << "\rSA-BFS | bc: " << best_cost << " cc: " << g.cost() << " t: " << T << " i: " << it << "/" << it << "                 " << endl;

	g.reinit_state(best_sol, best_cost);
	int64_t tmp = bfs_low_mem_local_search(g);

	cerr << "sa-bfs cost: " << best_cost << " ls cost: " << tmp << endl;

	if (tmp < best_cost)
	{
		best_cost = tmp;
		best_sol = g.sol();
	}
	return best_cost;
}

// Simulated annealing that tries to move *each* vertex to a random cluster
// When a vertex is moved, try to move its neighbors to the same cluster next
int64_t sa_bfs2(Instance &g, int64_t it)
{
	int n = g.n;
	int64_t best_cost = g.cost();
	vector<int> best_sol = g.sol();
	double T = initial_T(best_cost);

	// Randomness sources
	random_device rd;
	mt19937 rng(rd());
	// uniform_int_distribution unif(0, n - 1);
	uniform_real_distribution zero_one_unif(0.0, 1.0);

	vector<int> order(n);
	for (int i = 0; i < n; ++i)
		order[i] = i;

	// BFS variables
	vector<bool> seen(n, false);
	queue<int> q;

	int v, c;
	int64_t delta;
	for (int64_t i = 0; i < it; ++i)
	{
		shuffle(order.begin(), order.end(), rng);
		fill(seen.begin(), seen.end(), false);

		for (int j : order)
		{
			c = g.random_neighboring_cluster(j, rng);
			q.push(j);
			while (!q.empty())
			{
				v = q.front();
				q.pop();
				if (seen[v])
					continue;
				seen[v] = true;
				delta = g.delta_cost(v, c);
				if (accept(T, delta, zero_one_unif(rng)))
				{
					g.move_with_delta(v, c, delta);
					if (g.cost() < best_cost)
					{
						best_cost = g.cost();
						best_sol = g.sol();
					}
					// try to move neighbors next
					for (int u : g.neighbors(v))
					{
						q.push(u);
					}
				}

				if constexpr (DEBUG_SA)
					if (i % 10'000 == 0)
						cerr << "\rSA-BFS2 | bc: " << best_cost << " cc: " << g.cost() << " t: " << T << " i: " << i << "/" << it << "                 ";
				++i;
				T *= decay_rate;
			}
		}
	}
	if constexpr (DEBUG_SA)
		cerr << "\rSA-BFS2 | bc: " << best_cost << " cc: " << g.cost() << " t: " << T << " i: " << it << "/" << it << "                 " << endl;

	int64_t tmp = bfs_low_mem_local_search(g);

	cerr << "sa-bfs2 cost: " << best_cost << " ls cost: " << tmp << endl;

	if (tmp < best_cost)
	{
		best_cost = tmp;
		best_sol = g.sol();
	}
	else
	{
		g.reinit_state(best_sol, best_cost);
		tmp = bfs_low_mem_local_search(g);

		if (tmp < best_cost)
		{
			best_cost = tmp;
			best_sol = g.sol();
		}
	}

	cerr << "sa-bfs2 cost: " << best_cost << " ls cost: " << tmp << endl;

	return best_cost;
}

int64_t destroy_repair_sa(Instance &g, int64_t it, int ndestroy, double Tinit)
{
	if constexpr (DEBUG_SA)
		cerr << "Starting rand DRSA with nd: " << ndestroy << endl;

	int n = g.n;
	int64_t best_cost = g.cost();
	double T = Tinit;

	random_device rd;
	mt19937 rng(rd());
	uniform_int_distribution unif(0, n - 1);
	uniform_int_distribution unif_destroy(1, ndestroy);
	uniform_real_distribution zero_one_unif(0.0, 1.0);

	vector<int> vs(ndestroy);
	vector<int> old_cluster_of_vs(ndestroy);

	for (int64_t i = 0; i < it; ++i)
	{
		vs.clear(), old_cluster_of_vs.clear();
		for (int j = unif_destroy(rng); j > 0; --j)
		{
			int v = unif(rng);
			vs.push_back(v);
			old_cluster_of_vs.push_back(g.sol()[v]);
		}

		int64_t old_cost = g.cost();
		g.destroy_greedy_repair(vs);

		int64_t delta = g.cost() - old_cost;
		if (accept(T, delta, zero_one_unif(rng)))
		{
			if (g.cost() < best_cost)
				best_cost = g.cost();
		}
		else
			g.revert_cluster_of_with_cost(vs, old_cluster_of_vs, old_cost);

		if constexpr (DEBUG_SA)
			if (i % 100 == 0)
				cerr << "\rrand DRSA | bc: " << best_cost << " cc: " << g.cost() << " t: " << T << " i: " << i << "                 ";

		T *= decay_rate;
	}
	if constexpr (DEBUG_SA)
		cerr << endl;

	return best_cost;
}


void bfs_destroy_repair_sa_one_it(Instance &g, int n, int ndestroy, mt19937 &rng, vector<bool> &seen, double T, int64_t &best_cost)
{

	uniform_int_distribution unif(0, n - 1);
	uniform_int_distribution unif_destroy(1, ndestroy);
	uniform_real_distribution zero_one_unif(0.0, 1.0);
	vector<int> vs;
	vector<int> old_cluster_of_vs;

	if (true || zero_one_unif(rng) < 0.9)
	{
		fill(seen.begin(), seen.end(), false);
		// int nv = unif_destroy(rng) / 20;
		int nv = 5 + unif_destroy(rng);
		// nv = nv / 4 + nv * 3 * i / (4 * it);
		g.bfs_fill_vs(unif(rng), nv, vs, old_cluster_of_vs, seen);
	}
	else
	{
		for (int j = unif_destroy(rng); j > 0; --j)
		{
			int v = unif(rng);
			vs.push_back(v);
			old_cluster_of_vs.push_back(g.sol()[v]);
		}
	}

	int64_t old_cost = g.cost();
	g.destroy_greedy_repair(vs);

	int64_t delta = g.cost() - old_cost;
	if (accept(T, delta, zero_one_unif(rng)))
	{
		if (g.cost() < best_cost)
			best_cost = g.cost();
	}
	else
		g.revert_cluster_of_with_cost(vs, old_cluster_of_vs, old_cost);
}


bool bfs_destroy_repair_sa_one_it_order(
	Instance &g, int n, int ndestroy, mt19937 &rng, 
	vector<int> &order, vector<bool> &seen, double T, vector<int> &best_sol, int64_t &best_cost)
{
	uniform_int_distribution unif_destroy(1, ndestroy);
	uniform_real_distribution zero_one_unif(0.0, 1.0);
	vector<int> vs;
	vector<int> old_cluster_of_vs;

	shuffle(order.begin(), order.end(), rng);
	fill(seen.begin(), seen.end(), false);

	bool res = false;

	for (int i = 0; i < n; ++i)
	{
		if (seen[order[i]])
			continue;

		vs.clear();
		old_cluster_of_vs.clear();
		int nv = unif_destroy(rng) + 5;
		g.bfs_fill_vs(order[i], nv, vs, old_cluster_of_vs, seen);

		if (vs.size() == 0)
			continue;

		int64_t old_cost = g.cost();
		g.destroy_greedy_repair(vs);

		int64_t delta = g.cost() - old_cost;
		if (accept(T, delta, zero_one_unif(rng)))
		{
			if (g.cost() < best_cost)
			{
				res = true;
				best_cost = g.cost();
				best_sol = g.sol();
			}
		}
		else
		{
			g.revert_cluster_of_with_cost(vs, old_cluster_of_vs, old_cost);
		}
	}
	return res;
}


int64_t bfs_destroy_repair_sa(Instance &g, int64_t it, int ndestroy, double Tinit)
{
	if constexpr (DEBUG_SA)
		cerr << "Starting bfs DRSA with nd: " << ndestroy << endl;

	int n = g.n;
	int64_t best_cost = g.cost();
	double T = Tinit;

	random_device rd;
	mt19937 rng(rd());

	vector<bool> seen(n, false);

	for (int64_t i = 0; i < it; ++i)
	{
		bfs_destroy_repair_sa_one_it(g, n, ndestroy, rng, seen, T, best_cost);

		if constexpr (DEBUG_SA)
			if (i % 100 == 0)
				cerr << "\rbfs DRSA | bc: " << best_cost << " cc: " << g.cost() << " t: " << T << " i: " << i << "                 ";

		T *= decay_rate;
	}
	if constexpr (DEBUG_SA)
		cerr << endl;

	return best_cost;
}


int64_t move_neighbors_sa(Instance &g, int64_t it, double Tinit)
{
	int n = g.n;
	int64_t best_cost = g.cost();
	double T = Tinit;

	random_device rd;
	mt19937 rng(rd());

	uniform_int_distribution unif(0, n - 1);
	uniform_real_distribution zero_one_unif(0.0, 1.0);

	vector<int> vs, old_cluster_of_vs;

	for (int64_t i = 0; i < it; ++i)
	{
		old_cluster_of_vs.clear();

		int iv = unif(rng);

		vs = g.neighbors(iv);
		vs.push_back(iv);
		for (int v : vs)
			old_cluster_of_vs.push_back(g.sol()[v]);

		int64_t old_cost = g.cost();
		g.greedy_move_many(vs, false);

		int64_t delta = g.cost() - old_cost;
		if (accept(T, delta, zero_one_unif(rng)))
		{
			if (g.cost() < best_cost)
				best_cost = g.cost();
		}
		else
			g.revert_cluster_of_with_cost(vs, old_cluster_of_vs, old_cost);

		if constexpr (DEBUG_SA)
			if (i % 100 == 0)
				cerr << "\rmove_neigh_SA | bc: " << best_cost << " cc: " << g.cost() << " t: " << T << " i: " << i << "                 ";

		T *= decay_rate;
	}

	if constexpr (DEBUG_SA)
		cerr << endl;

	return best_cost;
}

int64_t move_neighbors_same_cc_sa(Instance &g, int64_t it, double Tinit)
{
	int n = g.n;
	int64_t best_cost = g.cost();
	double T = Tinit;

	random_device rd;
	mt19937 rng(rd());

	uniform_int_distribution unif(0, n - 1);
	uniform_real_distribution zero_one_unif(0.0, 1.0);

	vector<int> vs;

	for (int64_t i = 0; i < it; ++i)
	{
		int iv = unif(rng);

		int civ = g.sol()[iv];
		vs.clear();
		vs.push_back(iv);
		for (int v : g.neighbors(iv))
			if (g.sol()[v] == civ)
				vs.push_back(v);

		int64_t old_cost = g.cost();
		g.greedy_move_many(vs, false);

		int64_t delta = g.cost() - old_cost;
		if (accept(T, delta, zero_one_unif(rng)))
		{
			if (g.cost() < best_cost)
				best_cost = g.cost();
		}
		else
			g.revert_single_cluster_of_with_cost(vs, civ, old_cost);

		if constexpr (DEBUG_SA)
			if (i % 100 == 0)
				cerr << "\rmove_neigh_same_cc_SA | bc: " << best_cost << " cc: " << g.cost() << " t: " << T << " i: " << i << "                 ";

		T *= decay_rate;
	}

	if constexpr (DEBUG_SA)
		cerr << endl;

	return best_cost;
}
