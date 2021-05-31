#pragma once

#include <string>
#include <unordered_map>
#include <chrono>

using namespace std;
using namespace std::chrono;

#define PROFILE_FUNC Profiler tmp(__FUNCTION__);
#define PROFILE_NAME(s) Profiler tmp(s);

class Profiler
{
private:
	inline static unordered_map<string, pair<double, int64_t>> timings;

	string name;
	time_point<high_resolution_clock> start;
public:
	Profiler(const string &name);
	~Profiler();
	
	static void print_timings();
};