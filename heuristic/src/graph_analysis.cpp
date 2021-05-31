#include "low_mem_graph.h"
#include "local_search.h"
#include "multi_cc_instance.h"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

using Edge = pair<int,int>;

// Outputs: name, n, m, min/max deg, nb CC and sizes
void analyze(string fname)
{
	auto i = KernelizedMultiCCInstance::from_file(fname);
	cout << fname << "\t" << 0 << "\t" << i.m();
	cout << "\t" << i._n_cc() << "\t";
	cout << "[";
	vector<int> tmp;
	for (auto &cc: i._ccs())
		if (cc.size() > 5)
			cout << cc.size() << ", ";
	cout << "]";
	cout << endl;
}


int main()
{
	cout << "name\tn\tm\tnb CC\tCC sizes" << endl;
	int instance_i = 1;
	for (int i = instance_i; i < 200; i += 2)
	{
		string s;
		if (i < 100)
			s += "0";
		if (i < 10)
			s += "0";
		s += to_string(i);
		// string fname = "./exact/exact" + s + ".gr";
		string fname = "./heur/heur" + s + ".gr";
		// cerr << "starting " << fname << endl;
		analyze(fname);
	}
	return 0;
}

