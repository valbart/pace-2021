#include <vector>
#include <assert.h>
#include "time.h"

#include "../shared/graph.hpp"
#include "../shared/solution.hpp"
#include "instance.hpp"
#include "lp_solve.hpp"
#include "../heuristic/src/multi_cc_instance.h"
#include "../heuristic/src/union_find.h"
#include "../heuristic/src/local_search.h"
#include "../heuristic/src/simulated_annealing.h"

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

Solution run_algo(Graph *g)
{
    // Create instance from graph
    ExactInstance instance(*g, g->nb_edges(), 0);

    // Preprocess
    instance.preprocess();

    // Create instance for heuristic
    KernelizedMultiCCInstance heuristic(g->nr_vertices(), g->all_edges());
    auto timeout = [](int i) { return i < 5000; };
    double initial_t = 20;
    double decay_r = 0.99975;
    auto temp = [&](double cur_T, int64_t it, int64_t n) -> double {
        it = it % 150;
        double res = initial_t / (1 + it + log(n));
        return res;
    };

    // Compute heuristic
    heuristic.sa_multi_cc(timeout, temp, 40, 50, initial_t, 0);

    // Get Solution from heuristic
    int ub_heuristic = heuristic.count_sol();
    vector<Cluster> cluster_heuristic = heuristic.get_sol();
    Solution heur = Solution(ub_heuristic, g->nr_vertices(), cluster_heuristic);

    // Set upper bound for the instance
    instance.set_upper_bound(heur.get_cost());

    // Initialize lower bound
    clock_t start = clock();
    instance.init_lower_bound();
    clock_t end = clock();

    // Recompute lower bound if enough time left
    double time_for_lb = double(end - start) / CLOCKS_PER_SEC;
    int nr_recompute_lb = 0;
    if (time_for_lb <= 3)
        nr_recompute_lb = 50;
    else if (time_for_lb > 3 && time_for_lb <= 5)
        nr_recompute_lb = 30;
    else if (time_for_lb > 5 && time_for_lb <= 20)
        nr_recompute_lb = 10;
    else if (time_for_lb > 20 && time_for_lb <= 60)
        nr_recompute_lb = 4;

    if (instance.get_upper_bound() == instance.lower_bound())
        nr_recompute_lb = 0;

    for (int i = 0; i < nr_recompute_lb; i++)
    {
        instance.shuffle_graph();
        instance.recompute_lower_bound();
    }

    //cout << "LB: " << instance.lower_bound() << " UB: " << instance.get_upper_bound() << endl;

    // Apply forced moves
    instance.forced_moves();

    // Compute solution
    Solution res = instance.algo_brute();

    // Get optimal solution
    if (res.get_cost() < heur.get_cost())
    {
        heur = res;
    }
    return heur;
}

int main(int argc, char *argv[])
{
    cin.tie(0);
    ios::sync_with_stdio(false);
    Graph *G = Graph::load_cin(cin);

    vector<Graph> cc = dfs(G);
    for (Graph g : cc)
    {
        Graph save(g);
        /*Trying LP*/


        g.shuffle();


        if (save.nr_vertices() > 1 && save.nr_vertices() < 120) // < 120
        {
            Solution lp_sol = lp_solve(&g);
            if (lp_sol.get_cost() != std::numeric_limits<int>::max())
            {
                lp_sol.print(&save);
                continue;
            }
        }
        /*Else go normal*/
        Solution sol = run_algo(&g);
        sol.print(&save);
    }
    return 0;
}
