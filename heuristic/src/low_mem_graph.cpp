#include "low_mem_graph.h"
#include "union_find.h"

#include <vector>
#include <cassert>
#include <iostream>
#include <string>
#include <fstream>
#include <algorithm>

LowMemGraph::LowMemGraph(int n, const vector<Edge> &edges): adjs(n), _n(n), _m(0)
{
	for (auto &[u,v]: edges)
		add_edge(u, v);

	sort_adjs();
}

void LowMemGraph::add_edge(int u, int v)
{
	assert(u != v);
	adjs[u].push_back(v);
	adjs[v].push_back(u);
	++_m;
}

void LowMemGraph::sort_adjs()
{
	for (int i = 0; i < _n; ++i)
		sort(adjs[i].begin(), adjs[i].end());
}

bool LowMemGraph::remove_excess_degree_one()
{
	bool deleted = false;
	bool delete_next;

	vector<bool> del(_n, false);
	for (int u = 0; u < _n; ++u)
	{
		if (deg(u) == 1)
			continue;

		// delete all but one neigbor of degree 1 of v
		delete_next = false;
		for (int v : adjs[u])
		{
			if (del[v])
				continue;
			if (deg(v) == 1)
			{
				if (delete_next)
				{
					del[v] = true;
					deleted = true;
				}
				else
					delete_next = true;
			}
		}
	}

	_m = 0;
	for (int u = 0; u < _n; ++u)
	{
		if (del[u])
			adjs[u].clear();
		else
			adjs[u].erase(remove_if(adjs[u].begin(), adjs[u].end(), [&](int x) { return del[x]; }), adjs[u].end());
		_m += adjs[u].size();
	}
	_m /= 2;

	return deleted;
}


// ajdacency lists are assumed to be sorted in increasing order 
int LowMemGraph::neighborhoods_intersection_size(int u, int v)
{
	int r = 0;
	uint i = 0, j = 0;
	const vector<int> &nu = adjs[u];
	const vector<int> &nv = adjs[v];
	while (i < nu.size() && j < nv.size())
	{
		if (nu[i] == nv[j]) ++r;
		if (nu[i] < nv[j])
			++i;
		else
			++j;
	}
	return r;
}

void LowMemGraph::remove_edge(int u, int v)
{

	adjs[u].erase(remove(adjs[u].begin(), adjs[u].end(), v), adjs[u].end());
	adjs[v].erase(remove(adjs[v].begin(), adjs[v].end(), u), adjs[v].end());
	--_m;
}
// ajdacency lists are assumed to be sorted in increasing order 
bool LowMemGraph::disjoint_neighborhoods(int u, int v)
{
	return neighborhoods_intersection_size(u, v) == 0;
}

bool LowMemGraph::remove_edge_disjoint_neighbors()
{
	bool deleted = false;
	vector<int> candidate(_n, -1);

	// Mark vertices that have a degree 1 adjacent or two degree 2 that are adjacent, adjacents
	for (int u = 0; u < _n; ++u)
	{
		if (deg(u) == 1)
		{
			int v = adjs[u][0];
			candidate[v] = u;
			continue;
		}
		if (deg(u) == 2)
		{
			int v = adjs[u][0], w = adjs[u][1];
			if (deg(v) != 2)
			{
				if (deg(w) == 2)
					swap(v, w);
				else
					continue;
			}
			// v has degree 2
			int x = adjs[v][0], y = adjs[v][1];
			if (y == u)
				swap(x, y);
			// x == u
			if (w == y)
				candidate[w] = u;
		}
	}

	vector<int> tmp;
	for (int u = 0; u < _n; ++u)
	{
		if (candidate[u] == -1)
			continue;
		tmp.clear();
		for (int v : adjs[u])
		{
			if (v != candidate[u] && disjoint_neighborhoods(u, v))
			{
				auto it = lower_bound(adjs[v].begin(), adjs[v].end(), u);
				adjs[v].erase(it);
			}
			else
			{
				tmp.push_back(v);
			}
		}
		adjs[u] = tmp;
	}


	_m = 0;
	for (int u = 0; u < _n; ++u)
	{
		_m += deg(u);
	}
	_m /= 2;

	return deleted;
}

// If 2 degree 2 vertices v,w are
// adjacent to u,x that are not adjacent,
// remove two non adjacent edges in this c4. 
bool LowMemGraph::remove_C4()
{
	// vector<Edge> to_delete;
	bool deleted = false;
	int x, y;
	for (int u = 0; u < _n; ++u)
	{
		for (int v: adjs[u])
		{
			if (deg(v) != 2)
				continue;
			// x is the other neighbor of v
			if ((x = adjs[v][0]) == u)
				x = adjs[v][1];

			// u and x must be non-adjacent
			if (adjacent(u, x))
				continue;

			for (int w: adjs[u])
			{
				if (deg(w) != 2 || w == v)
					continue;

				if ((y = adjs[w][0]) == u)
					y = adjs[w][1];
				
				if (y != x)
					continue;

				remove_edge(u, v);
				remove_edge(w, x);
				deleted = true;
				break;
			}
		}
	}

	// sort(to_delete.begin(), to_delete.end());

	// for (int u = 0; u < _n; ++u)
	// {
	// 	adjs[u].erase(remove_if(adjs[u].begin(), adjs[u].end(),
	// 					[&](int z)
	// 					{
	// 						return binary_search(to_delete.begin(), to_delete.end(), make_pair(min(u, z), max(u, z)));
	// 					}
	// 					), adjs[u].end());
	// }

	int64_t m = 0;
	for (int u = 0; u < _n; ++u)
		m += deg(u);
	m /= 2;
	assert(m == _m);

	return deleted;
}


// If 3 degree <= 3 vertices u,v,w form a triangle
// which is not in any diamond,
// isolate them.
bool LowMemGraph::remove_deg3_triangles()
{
	bool deleted = false;
	vector<Edge> to_delete;
	int a, b, c;
	int nb, min_d;
	for (int u = 0; u < _n; ++u)
	{
		if (deg(u) != 3)
			continue;

		to_delete.clear();

		a = adjs[u][0];
		b = adjs[u][1];
		c = adjs[u][2];
		
		min_d = min(deg(a), min(deg(b), deg(c)));
		if (min_d > 3) // we need at least a degree 2 or 3 in neighbors
			continue;

		nb = adjacent(a, b) + adjacent(b, c) + adjacent(c, a);
		if (nb == 3) // K_4
		{
			if (deg(b) == 3)
				swap(a, b);
			if (deg(c) == 3)
				swap(a, c);
			// a has degree 3

			if (deg(b) <= 5 && deg(c) <= 5)
			{
				for (int x: neighbors(b))
					if (x != a && x != c && x != u)
						to_delete.emplace_back(b, x);

				for (int x: neighbors(c))
					if (x != a && x != b && x != u)
						to_delete.emplace_back(c, x);
			}

		}
		else if (nb == 2) // Diamond
		{
			if (deg(a) <= 3 && deg(b) <= 3 && deg(c) <= 3)
			{
				if (adjacent(a, b))
					swap(a, c);
				if (adjacent(a, b))
					swap(b, c);
				// a and b are not adjacent

				if (neighborhoods_intersection_size(a, b) == 2) // not a diamond
				{
					for (int x: neighbors(a))
						if (x != c && x != u)
							to_delete.emplace_back(a, x);

					for (int x: neighbors(b))
						if (x != c && x != u)
							to_delete.emplace_back(b, x);
				}
			}
		}
		else if (nb == 1) // Triangle
		{
			if (adjacent(a, c))
				swap(b, c);
			if (adjacent(b, c))
				swap(a, c);
			// a and b are adjacent, the others are not

			if (deg(a) <= 3 && deg(b) <= 3 && neighborhoods_intersection_size(a, b) == 1) // not a diamond
			{
				to_delete.emplace_back(u, c);

				for (int x: neighbors(a))
					if (x != b && x != u)
						to_delete.emplace_back(a, x);

				for (int x: neighbors(b))
					if (x != a && x != u)
						to_delete.emplace_back(b, x);
			}
		}

		for (auto &[x, y] : to_delete)
		{
			remove_edge(x, y);
			deleted = true;
		}
	}

	return deleted;
}


int64_t LowMemGraph::kernelize()
{
	int64_t m_removed = _m;
	int i = 0;
	bool cont;
	do
	{
		cont = false;
		while (remove_excess_degree_one())
			cont = true;

		while (remove_edge_disjoint_neighbors())
			cont = true;

		while (remove_C4())
			cont = true;

		while (remove_deg3_triangles())
			cont = true;

		++i;
	} while(cont);

	m_removed -= _m;
	return m_removed;
}

void LowMemGraph::debug() const
{
	cerr << "n: " << _n << ", m: " << _m << endl;
	for (int i = 0; i < _n; ++i)
	{
		cerr << i << ": ";
		for (int v : adjs[i])
			cerr << v << " ";

		cerr << endl;
	}
}

bool LowMemGraph::check_adj_lists() const
{
	for (int u = 0; u < _n; ++u)
	{
		if (!is_sorted(adjs[u].begin(), adjs[u].end()))
			return false;
		for (int v: adjs[u])
		{
			if (!binary_search(adjs[v].begin(), adjs[v].end(), u))
				return false;
		}
	}
	return true;
}


vector<vector<int>> LowMemGraph::connected_components()
{
	// Use an union find to build all the connected components
	// from edges without building a graph
	union_find uf(_n);
	for (int u = 0; u < _n; ++u)
		for (int v: adjs[u])
			uf.merge(u, v);
	
	int cc_count = 0;
	vector<int> cc_id_aux(_n);
	for (int i = 0; i < _n; ++i)
		if (uf.parent[i] == i)
			cc_id_aux[i] = cc_count++;

	vector<vector<int>> ccs(cc_count);
	for (int i = 0; i < _n; ++i)
		ccs[cc_id_aux[uf.find(i)]].push_back(i);

	return ccs;
}


/* INPUT FUNCTIONS */
// Does not handle comments
LowMemGraph LowMemGraph::from_istream(istream &is)
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
		edges.emplace_back(u - 1, v - 1);
	}
	
	return LowMemGraph(n, edges);
}

LowMemGraph LowMemGraph::from_cin()
{
	return from_istream(cin);
}

LowMemGraph LowMemGraph::from_file(const string &fname)
{
	ifstream graph_file(fname);
	return from_istream(graph_file);
}