#include "low_mem_graph.h"
#include "local_search.h"
#include "simulated_annealing.h"
#include "instance.h"
#include "multi_cc_instance.h"

#include <chrono>
#include <csignal>
#include <cstdlib>
#include <cstring>
#include <iostream>

using namespace std;
using namespace std::chrono;

constexpr bool DEBUG_STRAT = false;
constexpr bool DESCRIBE_OPTION = true;

// Signal handling
volatile sig_atomic_t tle = 0;

void term(int signum)
{
	(void) signum;
	tle = 1;
}

int main(int argc, char* argv[])
{
	int initial_it = 10;
	int ndestroy = 50;
	double Tinit = 20;
	int small_cc_threshold = 0;
	if constexpr (DESCRIBE_OPTION)
	{
		if (argc == 2 && argv[1] == string("--describe"))
		{
			cout << "KernelizedMultiCCInstance\n"
				 << "gs.sa_multi_cc([&](int i){ (void)i; return tle == 0; }, "
				 << initial_it << ", "
				 << ndestroy << ", "
				 << Tinit << ", "
				 << small_cc_threshold << ");"
				 << endl;
			return 0;
		}
	}

	// Signal handling
	struct sigaction action;
	memset(&action, 0, sizeof(struct sigaction));
	action.sa_handler = term;
	sigaction(SIGTERM, &action, NULL);

	// Fast IO
	cin.tie(0);
	ios::sync_with_stdio(false);

	auto gs = KernelizedMultiCCInstance::from_cin();
	auto timeout = [&](int i){ (void) i; return tle == 0; };
	auto temp 	 = [&](double cur_T, int64_t it, int64_t n) -> double
	{
		(void) cur_T;
		it = it % 2000;
		double res = Tinit / (1 + it + log(n));
		return res;
	};
	gs.sa_multi_cc(timeout, temp, initial_it, ndestroy, Tinit, small_cc_threshold);
	gs.print_sol();

	return 0;
}
