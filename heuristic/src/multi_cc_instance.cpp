#include "multi_cc_instance.h"
#include "union_find.h"
#include "local_search.h"
#include "simulated_annealing.h"

#include <iostream>
#include <fstream>
#include <functional>
#include <cassert>

/* Generates a LowMemGraph for each CC of the graph described by `edges`.
 * Generates additional info:
 * ccs[i]: list of the vertices in the i-th connected component.
 * vertex_to_cc[i] = {j,k} means that i is the k-th vertex of the j-th cc.
 * Hence the edge u,v should be added to the instance vertex_to_cc[u].first
 * between vertex_to_cc[u].second and vertex_to_cc[u].second. 
 */
KernelizedMultiCCInstance::KernelizedMultiCCInstance(int n, const vector<Edge> &edges): vertex_to_cc(n)
{
	// sort edges and format them as min/max for easy lookup
	initial_edges.reserve(edges.size());
	for (auto &[u,v]: edges)
		initial_edges.emplace_back(min(u, v), max(u, v));
	sort(initial_edges.begin(), initial_edges.end());

	// build initial graph and kernelize it
	LowMemGraph g(n, initial_edges);
	g.kernelize();

	// Compute cc info
	ccs = g.connected_components();
	// Sort them by decreasing size
	sort(ccs.begin(), ccs.end(), 
		[](const ConnectedComponent &cc1, const ConnectedComponent &cc2)
		{
			return cc1.size() > cc2.size();
		});
	int n_cc = ccs.size();
	for (int i = 0; i < n_cc; ++i)
		for (uint j = 0; j < ccs[i].size(); ++j)
			vertex_to_cc[ccs[i][j]] = {i, j};

	// build a list of edges and an instance graph for each cc
	vector<vector<Edge>> cc_edges(n_cc);
	for (auto &[u, v]: initial_edges)
	{
		auto &[cc_u, id_u] = vertex_to_cc[u];
		auto &[cc_v, id_v] = vertex_to_cc[v];
		if (cc_u == cc_v)
			cc_edges[cc_u].emplace_back(id_u, id_v);
	}

	for (int i = 0; i < n_cc; ++i)
		cc_instances.emplace_back(ccs[i].size(), cc_edges[i]);


	solutions.reserve(n_cc);
	costs.reserve(n_cc);
	total_cost = 0;
	for (auto &is: cc_instances)
	{
		solutions.push_back(is.sol());
		costs.push_back(is.cost());
		total_cost += is.cost();
	}

}

// Exits when exit_condition returns false, gets as input a number of iterations
void KernelizedMultiCCInstance::solve(function<bool(int)> exit_condition)
{
	int n_cc = cc_instances.size();
	int64_t res;
	int64_t it = 0;
	while (exit_condition(it))
	{
		for (int i = 0; i < n_cc; ++i)
		{
			int n = cc_instances[i].n;
			if (n <= 3)
				continue;
			cc_instances[i].reinit_all_zero();
			res = bfs_low_mem_local_search(cc_instances[i]);
			// res = move_neighbors_local_search_same_c(cc_graphs[i]);
			if (res < costs[i])
			{
				total_cost += res - costs[i];
				costs[i] = res;
				solutions[i] = cc_instances[i].sol();
			}
		}
		++it;
	}
}

void KernelizedMultiCCInstance::sa_multi_cc(
	function<bool(int)> exit_condition,
	function<double(double, int64_t, int64_t)> temp_iter,
	int initial_it,
	int ndestroy, double Tinit,
	int small_cc_threshold)
{
	int n_cc = cc_instances.size();
	int64_t res;

	random_device rd;
	mt19937 rng(rd());
	vector<vector<bool>> seen;
	for (auto &cc: cc_instances)
		seen.emplace_back(cc.n);

	vector<vector<int>> orders;

	for (int i = 0; i < n_cc; ++i)
	{
		auto &cc = cc_instances[i];
		orders.emplace_back(cc.n);
		for (int j = 0; j < cc.n; ++j)
			orders[i][j] = j;
	}


	vector<bool> interesting(n_cc, true);
	// first, do one iteration of local search on each cc
	// and mark thoses that we should not look at.
	for (int i = 0; i < n_cc; ++i)
	{
		auto &instance = cc_instances[i];
		int64_t n = instance.n, m = instance.m();
		if (n <= 3) break;
		if (m == ((n * (n - 1)) / 2))
		{
			interesting[i] = false;
			continue;
		}
		
		for (int j = 0; j < initial_it; ++j)
		{
			instance.reinit_all_zero();
			res = bfs_low_mem_local_search(instance);
			if (res < costs[i])
			{
				total_cost += res - costs[i];
				costs[i] = res;
				solutions[i] = instance.sol();
			}
		}
	}

	vector<double> T(n_cc, Tinit);
	vector<uint> last_improv(n_cc, 0);
	int64_t it = 0;
	while (exit_condition(it))
	{
		for (int i = 0; i < n_cc; ++i)
		{
			if (!interesting[i])
				continue;
			
			auto &cc = cc_instances[i];
			int n = cc.n;
			// Stop iteration when we reach small ccs
			if (n <= 3)
				break;

			total_cost -= costs[i];
			if (n > small_cc_threshold)
			{
				bool improv = bfs_destroy_repair_sa_one_it_order(
					cc, n, ndestroy, rng, 
					orders[i], seen[i], T[i], 
					solutions[i], costs[i]);
				if (!improv)
				{
					if (++last_improv[i] > 5*cc.n)
					{
						last_improv[i] = 0;
						cc.reinit_all_zero();
					}
				}
				else
					last_improv[i] = 0;
			}
			else
			{
				cc.reinit_all_zero();
				res = bfs_low_mem_local_search(cc);
				if (res < costs[i])
				{
					costs[i] = res;
					solutions[i] = cc.sol();
				}
			}
			total_cost += costs[i];
			if (costs[i] <= 1)
				interesting[i] = false;

			T[i] = temp_iter(T[i], it, n);
		}
		
		++it;
		//if constexpr (DEBUG_MULTI)
			//if (it % 1 == 0)
				//cerr << "\rmulti-ccbfs DRSA | cc: " << cost() << " t: " << T[0] << " i: " << it << "                 ";

	}
	//if constexpr (DEBUG_MULTI)
		//cerr << "\rmulti-ccbfs DRSA | cc: " << cost() << " t: " << T[0] << " i: " << it << "                 \n";

	for (int i = 0; i < n_cc; ++i)
		cc_instances[i].reinit_state(solutions[i], costs[i]);

}

void KernelizedMultiCCInstance::print_sol()
{
	int n_cc = cc_instances.size();
	// Compute edge to add
	for (int i = 0; i < n_cc; ++i)
	{	
		int n = cc_instances[i].n;
		Labeling &cluster_of = solutions[i];
		vector<vector<int>> clusters(n);
		for (int u = 0; u < n; ++u)
			clusters[cluster_of[u]].push_back(u);

		for (const auto &c: clusters)
			for (int x : c)
				for (int y : c)
				{
					int u = ccs[i][x];
					int v = ccs[i][y];
					if (u < v)
						if (!binary_search(initial_edges.begin(), 
											initial_edges.end(), 
											make_pair(u, v)))
							cout << u + 1 << " " << v + 1 << "\n";
				}
	}

	// Compute edges to delete
	for (auto &e: initial_edges)
	{
		auto &[u,v] = e;
		auto &[cc_u, id_u] = vertex_to_cc[u];
		auto &[cc_v, id_v] = vertex_to_cc[v];
		if ((cc_u != cc_v || (cc_u == cc_v && solutions[cc_u][id_u] != solutions[cc_v][id_v])))
			cout << u + 1 << " " << v + 1 << "\n";
	}

	cout.flush();
}


int64_t KernelizedMultiCCInstance::count_sol()
{
	int64_t res = 0;
	int n_cc = cc_instances.size();
	// Compute edge to add
	for (int i = 0; i < n_cc; ++i)
	{	
		int n = cc_instances[i].n;
		Labeling &cluster_of = solutions[i];
		vector<vector<int>> clusters(n);
		for (int u = 0; u < n; ++u)
			clusters[cluster_of[u]].push_back(u);

		for (const auto &c: clusters)
			for (int x : c)
				for (int y : c)
				{
					int u = ccs[i][x];
					int v = ccs[i][y];
					if (u < v)
						if (!binary_search(initial_edges.begin(), 
											initial_edges.end(), 
											make_pair(u, v)))
							++res;
				}
	}

	// Compute edges to delete
	for (auto &e: initial_edges)
	{
		auto &[u,v] = e;
		auto &[cc_u, id_u] = vertex_to_cc[u];
		auto &[cc_v, id_v] = vertex_to_cc[v];
		if ((cc_u != cc_v || (cc_u == cc_v && solutions[cc_u][id_u] != solutions[cc_v][id_v])))
			++res;
	}

	return res;
}

vector<Cluster> KernelizedMultiCCInstance::get_sol()
{
	vector<Cluster> res;
	int n_cc = cc_instances.size();
	
	// Concatenate clusters of each instance
	for (int i = 0; i < n_cc; ++i)
	{	
		int n = cc_instances[i].n;
		Labeling &cluster_of = solutions[i];
		vector<vector<int>> clusters(n);
		for (int u = 0; u < n; ++u)
			clusters[cluster_of[u]].push_back(ccs[i][u]);

		for (auto &c : clusters)
			res.emplace_back(move(c));
	}

	return res;
}

KernelizedMultiCCInstance KernelizedMultiCCInstance::from_istream(istream &is)
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

	return KernelizedMultiCCInstance(n, edges);
}

KernelizedMultiCCInstance KernelizedMultiCCInstance::from_cin()
{
	return from_istream(cin);
}

KernelizedMultiCCInstance KernelizedMultiCCInstance::from_file(const string &fname)
{

	ifstream graph_file(fname);
	return from_istream(graph_file);
}
