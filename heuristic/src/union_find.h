#pragma once

#include <vector>

class union_find 
{
public:
	std::vector<int> parent, count;
	
	union_find(int n);

	int find(int a); 
	bool merge(int a, int b);	
};
