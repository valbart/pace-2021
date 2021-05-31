#include <vector>
#include <assert.h>
#include "time.h"

#include "shared/graph.hpp"
#include "shared/solution.hpp"
#include "exact/instance.hpp"
#include "exact/lp_solve.hpp"
#include "heuristic/src/multi_cc_instance.h"
#include "heuristic/src/union_find.h"
#include "heuristic/src/local_search.h"
#include "heuristic/src/simulated_annealing.h"


using namespace std;

void dfs_aux(Graph *g, vector<int> &marks, Node v)
{
  for (Node w : g->neighbours(v))
  {
    if (marks[w] == -1)
    {
      marks[w] = marks[v];
      dfs_aux(g, marks, w);
    }
  }
}

vector<Graph> dfs(Graph *g)
{
  int n = g->nr_vertices();
  vector<int> marks(n, -1);
  int i = 0;
  for (int u = 0; u < n; u++)
  {
    if (marks[u] == -1)
    {
      marks[u] = i;
      dfs_aux(g, marks, u);
      i++;
    }
  }
  vector<Graph> cc;
  for (int j = 0; j < i; j++)
  {
    vector<Node> nodes;
    vector<Edge> edges;
    for (int u = 0; u < n; u++)
    {
      if (marks[u] == j)
      {
        nodes.push_back(u + 1);
        for (Node v : g->neighbours(u))
        {
          assert(marks[v] == j);
          if (u < v)
            edges.push_back(make_pair(u + 1, v + 1));
        }
      }
    }
    cc.push_back(Graph(nodes, edges));
  }
  return cc;
}

Solution run_algo(Graph* g)
{
  // Create instance from graph
  ExactInstance instance(*g, g->nb_edges(), 0);
  
  // Preprocess
  instance.preprocess();

  // Create instance for heuristic
  KernelizedMultiCCInstance heuristic(g->nr_vertices(), g->all_edges());
  auto timeout = [](int i){ return i < 5000; };
	double initial_t = 20;
	double decay_r = 0.99975;
	auto temp 	 = [&](double cur_T, int64_t it, int64_t n) -> double
	{
		it = it % 150;
		double res = initial_t / (1 + it + log(n));
		return res;
	};
  
  // Compute heuristic
	heuristic.sa_multi_cc(timeout, temp, 40, 50, initial_t, 0);

  // Get Solution from heuristic
	int ub_heuristic = heuristic.count_sol();
  vector<Cluster> cluster_heuristic = heuristic.get_sol();
  //cout << "UB" << " :" << heuristic.count_sol() << endl;
  Solution heur = Solution(ub_heuristic, g->nr_vertices(), cluster_heuristic);


  //ExactInstance copy(instance);


  //Solution heur = copy.heuristic();

  //return heur;

  //cout << "Heuristic done" << endl;
  
  // Set upper bound for the instance
  instance.set_upper_bound(heur.get_cost());
  

  // Initialize lower bound
  clock_t start = clock();
  instance.init_lower_bound();
  clock_t end = clock();
  cout << "First LB :" << instance.lower_bound() << endl;

  // Recompute lower bound if enough time left
  double time_for_lb = double(end-start)/CLOCKS_PER_SEC;
  int nr_recompute_lb = 0;
  if (time_for_lb <= 5)
    nr_recompute_lb = 20;
  else if (time_for_lb > 5 && time_for_lb <= 20)
    nr_recompute_lb = 10;
  else if (time_for_lb > 20 && time_for_lb <= 60)
    nr_recompute_lb = 4;

  cout << "Computing LB " << nr_recompute_lb << " times" << endl;

  for (int i = 0; i < nr_recompute_lb; i++)
  {
    instance.shuffle_graph();
    instance.recompute_lower_bound();
    //cout << "LB: " << instance.lower_bound() << " - UB: " << instance.get_upper_bound() << endl;
  }

  // Apply forced moves
  instance.forced_moves();

  //cout << "Starting algo" << endl;

  cout << "LB: " << instance.lower_bound() << " - UB: " << instance.get_upper_bound() << endl;

  // Compute solution
  Solution res = instance.algo_brute();

  // Get optimal solution
  if (res.get_cost() < heur.get_cost())
  {
    heur = res;
  }
  return heur;
}

void show_usage() {
	std::cout << "Usage:" << std::endl ;
	std::cout << "exact [file_name]" << std::endl ;
	std::cout << "\t- file_name : path to a file with a .gr extension." << std::endl ;
} ;

int main(int argc, char *argv[])
{
	if(argc != 2) {
		std::cerr << "Error: Incorrect number of arguments." << std::endl ;
		show_usage() ;
		return 1 ; 
	}

	string file = argv[1];
	Graph *G = Graph::load(file);
  
  
	int cost = 0;
  //Solution opt = lp_solve(G);
  //if (opt.get_cost() != std::numeric_limits<int>::max()) 
  if (false)
  {
    cout << "Try LP" << endl;
		//cost = opt.get_cost() ;
    cout << "LP found opt "  << cost << endl;
		//opt.print(G) ;
	} else 
  {
		//cout << "Branch search" << endl ;
		vector<Graph> cc = dfs(G);
		vector<int> cluster;
		for (Graph g : cc)
    {
			Graph save(g);
			Solution sol = run_algo(&g);
			cost += sol.get_cost();
			//sol.print(&save);
	  }
  }

	cout << "Final cost " << cost << endl;

	return 0;
}
