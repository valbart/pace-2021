#include "low_mem_graph.h"
#include "local_search.h"
#include "simulated_annealing.h"
#include "instance.h"
#include "multi_cc_instance.h"
#include "profiler.h"

#include <chrono>
#include <functional>
#include <iostream>
#include <cassert>
#include <map>

using namespace std;
using namespace std::chrono;

void test_kernelizing(string fname)
{
	auto gs = KernelizedMultiCCInstance::from_file(fname);

}

void test_multi(string fname)
{
	auto gs = KernelizedMultiCCInstance::from_file(fname);
	auto timeout = [](int i){ return i < 1000; };
	double initial_t = 20;
	double decay_r = 0.99975;
	auto temp 	 = [&](double cur_T, int64_t it, int64_t n) -> double
	{
		// if (cur_T < 0.01)
			// return initial_t;
		// return cur_T * decay_r;
		it = it % 1500;
		double res = initial_t / (1 + it + log(n));
		return res;
	};
	gs.sa_multi_cc(timeout, temp, 10, 50, initial_t, 0);
	cout << fname << " :" << gs.count_sol() << endl;
	// gs.print_sol();
}

int main()
{
	int instance_i = 199;
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

		cerr << "starting " << fname << ":\n";
		test_multi(fname);
		// test_kernelizing(fname);
		// test_greedy_on_small_cc(fname);
	}

	Profiler::print_timings();
	return 0;
}
