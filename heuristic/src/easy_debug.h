#pragma once

#include <vector>

using namespace std;
using Edge = pair<int,int>;

template <typename T>
inline void pv(const vector<T> &v)
{
	for (uint i = 0; i < v.size(); ++i)
	{
		cerr << v[i] << ", ";
	}
	cerr << endl;
}

inline void pv(const vector<Edge> &v)
{
	for (uint i = 0; i < v.size(); ++i)
	{
		cerr << "(" << v[i].first << ", " << v[i].second << ") " ;
	}
	cerr << endl;
}
