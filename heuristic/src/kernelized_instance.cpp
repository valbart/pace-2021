#include "kernelized_instance.h"
#include "simulated_annealing.h"

#include <fstream>

KernelizedInstance::KernelizedInstance(int n, const vector<Edge> &edges) : vertex_to_instance(n)
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
	auto ccs = g.connected_components();

	vector<int> vertex_to_cc(n);
	int n_cc = ccs.size();
	for (int i = 0; i < n_cc; ++i)
		for (int v : ccs[i])
			vertex_to_cc[v] = i;

	vector<int64_t> ccs_m(n_cc, 0);
	for (auto &[u, v]: initial_edges)
	{
		int cc_u = vertex_to_cc[u], cc_v = vertex_to_cc[v];
		if (cc_u == cc_v)
			++ccs_m[cc_u];
	}

	vector<bool> cc_interesting(n_cc);
	int instance_n = 0;
	for (int i = 0; i < n_cc; ++i)
	{
		int64_t m = ccs_m[i], n = (int64_t)ccs[i].size();
		cc_interesting[i] = ((n > 2) && (m != (n * (n-1) / 2)));
		if (cc_interesting[i])
			instance_n += n;
	}

	instance_to_vertex.resize(instance_n);
	int instance_v = 0;
	for (int v = 0; v < n; ++v)
	{
		if (!cc_interesting[vertex_to_cc[v]])
			vertex_to_instance[v] = -1 - vertex_to_cc[v];
		else
		{
			instance_to_vertex[instance_v] = v;
			vertex_to_instance[v] = instance_v;
			++instance_v;
		}
	}

	vector<Edge> instance_edges;
	for (auto &[u, v] : initial_edges)
	{
		int iu = vertex_to_instance[u], iv = vertex_to_instance[v];
		if ((iu >= 0) && (iv >= 0) && (vertex_to_cc[u] == vertex_to_cc[v]))
			instance_edges.emplace_back(iu, iv);
	}

	instance = Instance(instance_n, instance_edges);
}


int64_t KernelizedInstance::bfs_destroy_repair_sa(
	function<bool(int)> exit_condition,
	int ndestroy,
	double Tinit,
	double _decay_rate)
{
	random_device rd;
	mt19937 rng(rd());

	int64_t best_cost = instance.cost();
	vector<int> best_sol = instance.sol();
	int64_t res = best_cost;

	vector<bool> seen(instance.n);

	double T = Tinit;
	int it = 0;
	while (exit_condition(it))
	{
		bfs_destroy_repair_sa_one_it(instance, instance.n, ndestroy, rng, seen, T, res);

		if (it > 100'000 && res < best_cost)
			best_cost = res, best_sol = instance.sol();

		if constexpr (DEBUG_SA)
			if (it % 100 == 0)
				cerr << "\rker bfs DRSA | bc: " << best_cost << " cc: " << instance.cost() << " t: " << T << " i: " << it << "                 ";

		T *= _decay_rate;
		++it;
	}
	if constexpr (DEBUG_SA)
		cerr << endl;
		
	if (instance.cost() < best_cost)
		best_cost = instance.cost(), best_sol = instance.sol();

	instance.reinit_state(best_sol, best_cost);

	return best_cost;
}


void KernelizedInstance::print_sol()
{
	//compute edges to add
	int n = instance.n;
	const Labeling &cluster_of = instance.sol();
	vector<Cluster> clusters(n);
	for (int u = 0; u < n; ++u)
		clusters[cluster_of[u]].push_back(u);

	for (const auto &c : clusters)
	{
		for (int x : c)
			for (int y : c)
			{
				int u = instance_to_vertex[x], v = instance_to_vertex[y];
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
		auto &[u, v] = e;
		int iu = vertex_to_instance[u], iv = vertex_to_instance[v];
		if ((iu < 0 && iv < 0 && iu != iv)
			|| ((iu ^ iv) < 0)					// iu < 0 and iv >= 0 or iu >= 0 and iv < 0
			|| (iu >= 0 && iv >= 0 && cluster_of[iu] != cluster_of[iv]))
			cout << u + 1 << " " << v + 1 << "\n";
	}

	cout.flush();
}


int64_t KernelizedInstance::count_sol()
{
	int64_t res = 0;
	//compute edges to add
	int n = instance.n;
	const Labeling &cluster_of = instance.sol();
	vector<Cluster> clusters(n);
	for (int u = 0; u < n; ++u)
		clusters[cluster_of[u]].push_back(u);

	for (const auto &c : clusters)
	{
		for (int x : c)
			for (int y : c)
			{
				int u = instance_to_vertex[x], v = instance_to_vertex[y];
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
		auto &[u, v] = e;
		int iu = vertex_to_instance[u], iv = vertex_to_instance[v];
		if ((iu < 0 && iv < 0 && iu != iv)
			|| ((iu ^ iv) < 0)					// iu < 0 and iv >= 0 or iu >= 0 and iv < 0
			|| (iu >= 0 && iv >= 0 && cluster_of[iu] != cluster_of[iv]))
			++res;
	}

	return res;
}


KernelizedInstance KernelizedInstance::from_istream(istream &is)
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

	return KernelizedInstance(n, edges);
}


KernelizedInstance KernelizedInstance::from_cin()
{
	return from_istream(cin);
}


KernelizedInstance KernelizedInstance::from_file(const string &fname)
{
	ifstream graph_file(fname);
	return from_istream(graph_file);
}