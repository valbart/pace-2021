#include "profiler.h"

#include <iostream>

Profiler::Profiler(const string &c_name): name(c_name) 
{
	start = high_resolution_clock::now();
}

Profiler::~Profiler()
{
	const auto end = high_resolution_clock::now();
	duration<double, micro> ms = end - start;
	double mu_s = ms.count(); //duration_cast<double, milliseconds>(end - start).count();
	// cerr << name << " " << mu_s << endl;
	auto [it, b] = timings.try_emplace(name, 0, 0);
	auto &[s, p] = (*it);
	p.first  += mu_s;
	p.second += 1;
}

void Profiler::print_timings()
{
	for (auto &[k, p]: timings)
	{
		double avg = ((double)p.first)/((double)p.second);
		cerr << k << ": avg " << avg << "µs over " << p.second << " calls ("
			 << p.first/1'000'000 << "s)" << endl;
	}
}
