#include "instance.h"
#include "profiler.h"

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
#include <functional>
#include <cmath>

Instance::Instance(int64_t _n, const vector<Edge> &edges):
	g(_n, edges), 
	cluster_of(g.n(), 0), cluster_size(g.n(), 0), candidate_clusters_adj(g.n(), 0), n(g.n())
{
	cluster_size[0] = n;
	for (int i = 1; i < n; ++i)
		zero_size_cluster.push_back(i);
	_cost = ((int64_t)n * (n - 1)) / 2 - g.m();
}


// void Instance::print_sol(const vector<Edge> &initial_edges, function<int(int)> f_relabel)
// {

// 	vector<vector<int>> clusters(n);
// 	for (int u = 0; u < n; ++u)
// 		clusters[cluster_of[u]].push_back(u);

// 	for (const auto &c: clusters)
// 		for (int u : c)
// 			for (int v : c)
// 				if (u < v)
// 					// if (!g.is_initial(u, v))
// 					if (!binary_search(initial_edges.begin(), 
// 										initial_edges.end(), 
// 										make_pair(f_relabel(u), f_relabel(v))))
// 						cout << f_relabel(u) + 1 << " " << f_relabel(v) + 1 << "\n";

// 	for (auto &[u, v]: g.initial_edges_vec())
// 			if (cluster_of[u] != cluster_of[v])
// 				cout << f_relabel(u) + 1 << " " << f_relabel(v) + 1 << "\n";

// }

void Instance::reinit_all_zero()
{
	fill(cluster_size.begin(), cluster_size.end(), 0);
	fill(cluster_of.begin(), cluster_of.end(), 0);
	cluster_size[0] = n;

	zero_size_cluster.clear();
	for (int i = 1; i < n; ++i)
		zero_size_cluster.push_back(i);
	_cost = ((n * (n - 1)) / 2) - g.m();
}

// void Instance::reinit_all_alone()
// {
// 	for (int i = 0; i < n; ++i)
// 	{
// 		cluster_of[i] = i;
// 		cluster_size[i] = 1;
// 	}
// 	_cost = g.m();
// 	zero_size_cluster.clear();
// }

void Instance::reinit_state(const vector<int> &v, int64_t cost)
{
	for (int i = 0; i < n; ++i)
		cluster_size[i] = 0;
	for (int i = 0; i < n; ++i)
	{
		cluster_of[i] = v[i];
		cluster_size[v[i]]++;
	}
	_cost = cost;
	zero_size_cluster.clear();
	for (int i = 0; i < n; ++i)
	{
		if (cluster_size[i] == 0)
			zero_size_cluster.push_back(i);
	}
}

void Instance::debug() const
{
	for (int i = 0; i < n; ++i)
	{
		cout << i << " --> " << cluster_of[i] << " \tadjs: ";
		for (int v : g.neighbors(i))
			cout << v << " ";

		cout << endl;
	}
}


// amortized management of zero_size clusters
int Instance::get_zero_size()
{
	assert(!zero_size_cluster.empty());

	int res = zero_size_cluster.back();
	
	while (cluster_size[res] != 0 && !zero_size_cluster.empty())
	{
		zero_size_cluster.pop_back();
		res = zero_size_cluster.back();
	}
	return res;
}



// Move v and handle zero size clusters changes
void Instance::_move(int v, int c)
{
	int cv = cluster_of[v];
	if (cv == c)
		return;

	cluster_of[v] = c;
	cluster_size[cv]--;
	cluster_size[c]++;

	if (cluster_size[cv] == 0)
		zero_size_cluster.push_back(cv);

}


// Move v to the cluster that reduces the most the cost,
// if any.
// Returns True iff there was a move.
// O(d(v))
bool Instance::greedy_move(int v)
{
	int cu, cv = cluster_of[v];

	for (int u : g.neighbors(v))
	{
		cu = cluster_of[u];
		++candidate_clusters_adj[cu];
	}

	int64_t self_edges = candidate_clusters_adj[cv];
	int64_t best_cost = 0, cost;
	int best_cluster = -1;
	int64_t self_cost = -(cluster_size[cv] - 1 - 2 * self_edges);
	for (int u : g.neighbors(v))
	{
		cu = cluster_of[u];
		if (cu == cv)
			continue;
		cost = (cluster_size[cu] - 2 * candidate_clusters_adj[cu]) + self_cost;
		if (cost < best_cost)
		{
			best_cost = cost;
			best_cluster = cu;
		}
	}

	// Reset the content of candidate_clusters_adj
	for (int u : g.neighbors(v))
	{
		cu = cluster_of[u];
		candidate_clusters_adj[cu] = 0;
	}

	// Check whether it is better to put it in an empty cluster
	// Todo: if there is no empty cluster, create a new one ! (that cannot happen as long as we have n clusters)
	if (self_cost < best_cost && !zero_size_cluster.empty())
	{
		best_cost = self_cost;
		best_cluster = get_zero_size();
	}

	if (best_cluster == -1)
		return false;

	// Move v
	_move(v, best_cluster);
	_cost += best_cost;

	return true;
}


// Move vs, all at once to the same cluster, greedily
// Returns True iff there was a move.
// O(\sum_{v\in vs} d(v))
bool Instance::greedy_move_many(const vector<int> &vs, bool revert)
{
	if (vs.empty())
		return false;

	int64_t initial_cost = _cost;
	vector<int> old_clusters;
	if (revert)
	{
		old_clusters.resize(vs.size());
		for (size_t i = 0; i < vs.size(); ++i)
			old_clusters[i] = cluster_of[vs[i]];
	}

	int cu;
	int zero_c = cluster_of[vs[0]];
	if (!zero_size_cluster.empty())
		zero_c = get_zero_size();

	bool will_move = false;
	for (int v: vs)
		if (cluster_of[v] != zero_c)
			will_move = true;

	if (!will_move)
		return false;

	// strategy: move all to a zero_size cluster, then try to merge with a neighboring one.
	for (int v : vs)
		move_with_delta(v, zero_c, delta_cost(v, zero_c));

	for (int v : vs)
	{
		for (int u: g.neighbors(v))
		{
			cu = cluster_of[u];
			if (cu != zero_c)
				++candidate_clusters_adj[cu];
		}
	}

	int64_t best_delta = 0, tmp_delta;
	int best_cluster = -1;
	int64_t zero_c_size = cluster_size[zero_c];
	// assert(zero_c_size == vs.size());
	for (int v : vs)
	{
		for (int u: g.neighbors(v))
		{
			cu = cluster_of[u];
			if (cu == zero_c)
				continue;
			tmp_delta = (((int64_t)cluster_size[cu])* zero_c_size) - 2*candidate_clusters_adj[cu];
			if (tmp_delta < best_delta)
			{
				best_delta = tmp_delta;
				best_cluster = cu;
			}
		}
	}

	// reset
	for (int v : vs)
		for (int u: g.neighbors(v))
			candidate_clusters_adj[cluster_of[u]] = 0;

	// was it better initially ?
	if (initial_cost <= _cost + best_delta)
	{
		if (revert)
			revert_cluster_of_with_cost(vs, old_clusters, initial_cost);
		return false;
	}

	// was it better to put vs in a new cluster >
	if (best_delta >= 0) 
	{
		return true;
	}


	for (int v : vs)
		_move(v, best_cluster);

	_cost += best_delta;
	return true;
}


bool Instance::greedy_move_neighbors(int v, bool revert)
{
	vector<int> vs = g.neighbors(v);
	vs.push_back(v);
	return greedy_move_many(vs, revert);
}

bool Instance::greedy_move_neighbors_same_c(int v, bool revert)
{
	int cv = cluster_of[v];
	vector<int> vs = {v};
	
	for (int u: g.neighbors(v))
		if (cluster_of[u] == cv)
			vs.push_back(u);

	return greedy_move_many(vs, revert);
}

// Cost of moving v to c
// O(d(v)), unless c is the cluster of v, in which case it is O(1)
int64_t Instance::delta_cost(int v, int c)
{
	int cv = cluster_of[v];
	if (cv == c)
		return 0;
	
	int64_t self_edges = 0, to_edges = 0;
	for (int u : g.neighbors(v))
	{
		if (cluster_of[u] == cv)
			++self_edges;
		if (cluster_of[u] == c)
			++to_edges;
	}

	return (cluster_size[c] - 2 * to_edges) - (cluster_size[cv] - 1 - 2 * self_edges);
}

void Instance::move_with_delta(int v, int c, int64_t delta)
{
	_move(v, c);
	_cost += delta;
}

void Instance::move_to_zero_size(const vector<int> &vs)
{
	for (int v : vs)
	{
		if (cluster_size[cluster_of[v]] == 1)
			continue;
		int c = get_zero_size();
		move_with_delta(v, c, delta_cost(v, c));
	}
}


void Instance::revert_cluster_of_with_cost(const vector<int> &vs, const vector<int> &old_cluster_of_vs, int64_t old_cost)
{
	for (size_t i = 0; i < vs.size(); ++i)
	{
		int v = vs[i], c = old_cluster_of_vs[i];
		_move(v, c);
	}
	_cost = old_cost;
}

void Instance::revert_single_cluster_of_with_cost(const vector<int> &vs, int old_cluster_of_vs, int64_t old_cost)
{
	for (size_t i = 0; i < vs.size(); ++i)
	{
		int v = vs[i];
		_move(v, old_cluster_of_vs);
	}
	_cost = old_cost;
}

void Instance::destroy_greedy_repair(const vector<int> &vs)
{	
	move_to_zero_size(vs);
	for (int v : vs)
	{
		greedy_move(v);
	}
}

// Returns a random neighboring cluster, or a zero size cluster,
// uniformly at random among NEIGHBORS.
// Todo: uniform over neighboring clusters ?
int Instance::random_neighboring_cluster(int v, mt19937 &rng)
{
	// We can change -1 to -x to increase the proba of new cluster
	uniform_int_distribution<int> unif(-1, g.deg(v) - 1);
	int r = unif(rng);
	if (r < 0)
		return get_zero_size();
	else
		return cluster_of[g.neighbors(v)[r]];
}


void Instance::bfs_fill_vs(int v, unsigned int nv, vector<int> &vs, vector<int> &cluster_of_vs, vector<bool> &seen) const
{
	queue<int> q;
	q.push(v);
	while (!q.empty() && vs.size() < nv)
	{
		v = q.front(); q.pop();
		if (seen[v])
			continue;
		seen[v] = true;

		vs.push_back(v);
		cluster_of_vs.push_back(cluster_of[v]);
		for (int u : g.neighbors(v))
		{
			// if (seen[u])
				// continue;
			q.push(u);
			// seen[u] = true;
		}
	}
}


// DEBUG FUNCTIONS //

int64_t Instance::actual_cost()
{
	int64_t r = 0;
	for (int i = 0; i < n; ++i)
	{
		r += cluster_size[i] * (cluster_size[i] - 1);
	}
	for (int i = 0; i < n; ++i)
		for (int v : g.neighbors(i))
			if (cluster_of[v] != cluster_of[i])
				++r;
			else
				--r;
	return r / 2;
}

bool Instance::sanity_check()
{
	auto ac = actual_cost();
	if (ac != cost())
	{
		cerr << "cost pb" << cost() << " vs " << ac << endl;
		return false;
	}
	for (int i = 0; i < n; ++i)
	{
		int count = 0;
		for (int j = 0; j < n; ++j)
			if (cluster_of[j] == i)
				++count;
		if (count != cluster_size[i])
		{
			cerr << "pb at " << i << endl;
			return false;
		}
	}
	return true;
}


Instance Instance::from_istream(istream &is)
{
	int n, m, u, v;
	string s;
	is >> s >> s; // ignore header
	is >> n >> m;

	vector<Edge> edges;
	edges.reserve(m);
	for (int i = 0; i < m; ++i)
	{
		is >> u >> v;
		edges.emplace_back(u - 1, v - 1); // vertices are 0-indexed
	}

	return Instance(n, edges);
}

Instance Instance::from_cin()
{
	return from_istream(cin);
}

Instance Instance::from_file(const string &fname)
{

	ifstream graph_file(fname);
	return from_istream(graph_file);
}
