#include "union_find.h"

union_find::union_find(int n) : parent(n), count(n, 1) 
{
	for (int i = 0; i < n; ++i)
		parent[i] = i;
}

int union_find::find(int a) 
{ 
	while (parent[parent[a]] != parent[a])
	{
		parent[a] = parent[parent[a]];
	}
	return parent[a];
}

bool union_find::merge(int a, int b) 
{
	a = find(a); b = find(b);
	if (a == b) return false;
	if (count[a] < count[b]) std::swap(a,b);
	parent[b] = a; 
	count[a] += count[b];
	return true; 
}	