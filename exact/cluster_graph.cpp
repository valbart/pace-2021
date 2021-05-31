#include "cluster_graph.hpp"
#include <iostream>
#include <unordered_set>
#include <algorithm>
#include <assert.h>

// Auxiliary methods declaration. Implemented at end of file.
// Need to move that in another hpp file for clean code.
vector<int> intersection(vector<int> &nums1, vector<int> &nums2);
vector<int> set_minus(vector<int> a, vector<int> b);

ClusterGraph::ClusterGraph() {}

ClusterGraph& ClusterGraph::operator=(const ClusterGraph& rhs){
  this->g = rhs.g;
  this->n = rhs.n;
  this->node_to_cluster = rhs.node_to_cluster;
  this->clusters = rhs.clusters;
  this->cluster_edges_count = rhs.cluster_edges_count;
  this->marked_graph = rhs.marked_graph;
  this->cluster_adjacency = rhs.cluster_adjacency;
  this->cluster_adjacency_list = rhs.cluster_adjacency_list;
  this->degres = rhs.degres;
  this->twins = rhs.twins;
  this->cost = rhs.cost;
  this->star_value = rhs.star_value;
  this->nb_cluster_edges = rhs.nb_cluster_edges;
  this->upper_bound = rhs.upper_bound;
  return *this;
}

  
// Carefull: when looking for edges between clusters, first element should be < than second element
ClusterGraph::ClusterGraph(Graph g, int _cost, int upper_bound)
{
  this->g = g;
  this->n = g.nr_vertices();
  this->node_to_cluster = vector<ClusterId>(n);
  this->clusters = vector<Cluster>(n);
  this->cluster_edges_count = vector<vector<vector<int>>>(n);
  this->cluster_adjacency = vector<vector<bool>>(n);
  this->cluster_adjacency_list = vector<vector<ClusterId>>(n);
  this->marked_graph = vector<vector<Star>>(n);
  this->star_value = 0;
  this->degres = vector<int>(n, 0);
  this->cost = _cost;
  this->nb_cluster_edges = (n * (n - 1)) / 2;
  for (int u = 0; u < this->n; u++)
  {
    this->node_to_cluster[u] = u;
    this->clusters[u].push_back(u);
    this->cluster_edges_count[u] = vector<vector<int>>(n);
    this->cluster_adjacency[u] = vector<bool>(n, true);
    this->cluster_adjacency[u][u] = false;
    this->marked_graph[u] = vector<Star>(this->n);
    for (int v = 0; v < n; v++)
    {
      if (v != u)
        cluster_adjacency_list[u].push_back(v);
      this->marked_graph[u][v] = Star();
      if (v > u)
      {
        if (g.has_edge(u, v))
        {
          this->cluster_edges_count[u][v] = {0, 1, 0, 0};
          this->degres[u]++;
          this->degres[v]++;
        }
        else
          this->cluster_edges_count[u][v] = {0, 0, 0, 1};
      }
    }
  }
  this->twins = vector<vector<int>>(this->n, vector<int>(this->n));

  for (int u = 0; u < this->n; u++)
  {
    for (int v = u + 1; v < this->n; v++)
    {
      int a = 0;
      if (this->g.has_edge(u, v))
        a -= 2;
      for (Node w : this->g.neighbours(u))
      {
        if (!this->g.has_edge(w, v))
          a++;
      }
      for (Node w : this->g.neighbours(v))
      {
        if (!this->g.has_edge(w, u))
          a++;
      }
      this->twins[u][v] = a;
      this->twins[v][u] = a;
    }
  }
}

ClusterGraph::ClusterGraph(const ClusterGraph &old_cluster_graph)
{
  this->g = old_cluster_graph.g;
  this->n = old_cluster_graph.n;
  this->node_to_cluster = old_cluster_graph.node_to_cluster;
  this->clusters = old_cluster_graph.clusters;
  this->cluster_edges_count = old_cluster_graph.cluster_edges_count;
  this->marked_graph = old_cluster_graph.marked_graph;
  this->cluster_adjacency = old_cluster_graph.cluster_adjacency;
  this->cluster_adjacency_list = old_cluster_graph.cluster_adjacency_list;
  this->degres = old_cluster_graph.degres;
  this->twins = old_cluster_graph.twins;
  this->cost = old_cluster_graph.cost;
  this->star_value = old_cluster_graph.star_value;
  this->nb_cluster_edges = old_cluster_graph.nb_cluster_edges;
  this->upper_bound = old_cluster_graph.upper_bound;
}

ClusterGraph::ClusterGraph(const ClusterGraph &old_cluster_graph, Graph new_graph)
{
  this->g = new_graph;
  this->n = old_cluster_graph.n;
  this->node_to_cluster = old_cluster_graph.node_to_cluster;
  this->clusters = old_cluster_graph.clusters;
  this->cluster_edges_count = old_cluster_graph.cluster_edges_count;
  this->marked_graph = old_cluster_graph.marked_graph;
  this->star_value = old_cluster_graph.star_value;
  this->cluster_adjacency = old_cluster_graph.cluster_adjacency;
  this->cluster_adjacency_list = old_cluster_graph.cluster_adjacency_list;
  this->degres = old_cluster_graph.degres;
  this->twins = old_cluster_graph.twins;
  this->cost = old_cluster_graph.cost;
  this->nb_cluster_edges = old_cluster_graph.nb_cluster_edges;
  this->upper_bound = old_cluster_graph.upper_bound;
}

int ClusterGraph::get_cost() const
{
  return this->cost;
}

int ClusterGraph::degree(Node u)
{
  return this->degres[u];
}

bool ClusterGraph::false_twins(Node u, Node v)
{
  return (!g.has_edge(u, v) && (unsigned int) twins[u][v] == g.neighbours(u).size()+g.neighbours(v).size());
}

int ClusterGraph::twinness(Node u, Node v)
{
  if (this->g.has_edge(u, v))
  {
    return this->twins[u][v];
  }
  else
  {
    return numeric_limits<int>::max();
  }
}

pair<ClusterId, ClusterId> ClusterGraph::bestpair()
{
  pair<ClusterId, ClusterId> best(0, 1);
  pair<ClusterId, ClusterId> value(-1, -1);
  for (ClusterId cu = 0; cu < this->n; cu++)
  {
    for (ClusterId cv = 0; cv < this->n; cv++)
    {
      //if (cu < cv && has_edge(cu, cv))
      if (cu < cv && cluster_adjacency[cu][cv])
      {
        pair<ClusterId, ClusterId> new_value(0, 0);
        new_value.first = cluster_edges_count[cu][cv][0] + cluster_edges_count[cu][cv][1];
        if (new_value.first < value.first)
          continue;
        for (Node u : clusters[cu])
        {
          new_value.second += this->degres[u];
        }
        for (Node u : clusters[cv])
        {
          new_value.second += this->degres[u];
        }
        if (new_value.first > value.first || new_value.second >= value.second)
        {
          best = make_pair(cu, cv);
          value = new_value;
        }
      }
    }
  }
  return best;
}

int ClusterGraph::nr_edges(ClusterId cu, ClusterId cv)
{
  //if (!has_edge(cu, cv)) //TODO: immplem can_merge here
  if (!cluster_adjacency[cu][cv])
    return 0;
  ClusterId c1 = min(cu, cv);
  ClusterId c2 = max(cu, cv);
  return cluster_edges_count[c1][c2][0] + cluster_edges_count[c1][c2][1];
}

int ClusterGraph::nr_edges(ClusterId cu, ClusterId cv, bool marked)
{
  //if (!has_edge(cu, cv))
  if (!cluster_adjacency[cu][cv])
    return 0;
  ClusterId c1 = min(cu, cv);
  ClusterId c2 = max(cu, cv);
  if (marked)
    return cluster_edges_count[c1][c2][0];
  else
    return cluster_edges_count[c1][c2][1];
}

int ClusterGraph::nr_non_edges(ClusterId cu, ClusterId cv)
{
  //if (!has_edge(cu, cv))
  if (!cluster_adjacency[cu][cv])
    return clusters[cu].size() * clusters[cv].size();
  ClusterId c1 = min(cu, cv);
  ClusterId c2 = max(cu, cv);
  return cluster_edges_count[c1][c2][2] + cluster_edges_count[c1][c2][3];
}

int ClusterGraph::nr_non_edges(ClusterId cu, ClusterId cv, bool marked)
{
  //if (!has_edge(cu, cv))
  if (!cluster_adjacency[cu][cv])
    return clusters[cu].size() * clusters[cv].size();
  ClusterId c1 = min(cu, cv);
  ClusterId c2 = max(cu, cv);
  if (marked)
    return cluster_edges_count[c1][c2][2];
  else
    return cluster_edges_count[c1][c2][3];
}

ClusterId ClusterGraph::cluster_Id_of_node(Node u) const
{
  return this->node_to_cluster[u];
}

const vector<Node> ClusterGraph::cluster(ClusterId i) const
{
  return this->clusters[i];
}

vector<ClusterId> ClusterGraph::neighboring_clusters(ClusterId i) const
{
  /*vector<ClusterId> c_neighbors(0);
  for (ClusterId cu = 0; cu < n; cu++)
    if (cluster_adjacency[i][cu])
      c_neighbors.push_back(cu);
  return c_neighbors;*/
  return this->cluster_adjacency_list[i];
}

void ClusterGraph::add_edge(ClusterId i, ClusterId j)
{
  cluster_adjacency[i][j] = true;
  cluster_adjacency[j][i] = true;
  this->cluster_adjacency_list[i].push_back(j);
  this->cluster_adjacency_list[j].push_back(i);
  this->nb_cluster_edges++;
}

void ClusterGraph::remove_edge(ClusterId i, ClusterId j)
{
  cluster_adjacency[i][j] = false;
  cluster_adjacency[j][i] = false;

  auto pos_i = find(cluster_adjacency_list[j].begin(), cluster_adjacency_list[j].end(), i);
  if (pos_i != cluster_adjacency_list[j].end()) {
    cluster_adjacency_list[j].erase(pos_i);
  }
    
  auto pos_j = find(cluster_adjacency_list[i].begin(), cluster_adjacency_list[i].end(), j);
  if (pos_j != cluster_adjacency_list[i].end()) {
    cluster_adjacency_list[i].erase(pos_j);
  }

  this->nb_cluster_edges--;
}

/*inline bool ClusterGraph::has_edge(ClusterId i, ClusterId j) const
{
  return cluster_adjacency[i][j];
}*/

const Cluster ClusterGraph::get_cluster(ClusterId i) const
{
  return this->clusters[i];
}

vector<pair<Node, Node>> ClusterGraph::merge(ClusterId i, ClusterId j, bool remark)
{
  vector<ClusterId> neighbors_i = neighboring_clusters(i);
  vector<ClusterId> neighbors_j = neighboring_clusters(j);
  vector<ClusterId> common_n = intersection(neighbors_i, neighbors_j);
  vector<ClusterId> private_i = set_minus(neighbors_i, neighbors_j);
  vector<ClusterId> private_j = set_minus(neighbors_j, neighbors_i);
  vector<pair<Node, Node>> to_remark;
  for (auto &it : private_i)
  {
    if (it == j) continue;
    vector<pair<Node,Node>> extension = this->split(i, it, false);
    to_remark.insert(to_remark.end(), extension.begin(), extension.end());
  }

  for (auto &it : private_j)
  {
    if (it == i) continue;
    vector<pair<Node,Node>> extension = this->split(j, it, false);
    to_remark.insert(to_remark.end(), extension.begin(), extension.end());
  }

  for (auto &it : common_n)
  {
    this->remove_edge(it, j);
  }

  for (auto &it : common_n)
  {
    for (int s = 0; s < 4; s++)
    {
      this->cluster_edges_count[min(i, it)][max(i, it)][s] +=
          this->cluster_edges_count[min(j, it)][max(j, it)][s];
    }
  }

  for (auto &it_i : this->clusters[i])
  {
    for (auto &it_j : this->clusters[j])
    {
      vector<pair<Node, Node>> extension = this->add_graph_edge(it_i, it_j);
      to_remark.insert(to_remark.end(), extension.begin(), extension.end());
      unmark(it_i, it_j);
    }
  }
  this->remove_edge(i, j);

  for (auto &it_j : this->clusters[j])
    this->node_to_cluster[it_j] = i;
  this->clusters[i].insert(this->clusters[i].end(), this->clusters[j].begin(), this->clusters[j].end());

  this->clusters[j].clear();

  if (remark) 
  {
    for (pair<Node, Node> p: to_remark) try_mark(p.first, p.second, false, true);
    to_remark.clear();
  }

  return to_remark;
}

vector<pair<Node, Node>> ClusterGraph::split(ClusterId i, ClusterId j, bool remark)
{
  vector<pair<Node,Node>> to_remark (0);
  //if (has_edge(i, j))
  if (cluster_adjacency[i][j])
  {
    for (auto &vertex_in_i : this->clusters[i])
    {
      for (auto &vertex_in_j : this->clusters[j])
      {
        vector<pair<Node, Node>> extension = this->remove_graph_edge(vertex_in_i, vertex_in_j);
        to_remark.insert(to_remark.end(), extension.begin(), extension.end());
        unmark(vertex_in_i, vertex_in_j);
      }
    }
    this->remove_edge(i, j);
  }

  if (remark) 
  {
    for (pair<Node, Node> p: to_remark) try_mark(p.first, p.second, false, true);
    to_remark.clear();
  }

  return to_remark;
}

vector<pair<Node, Node>> ClusterGraph::add_graph_edge(Node u, Node v)
{
  assert(!is_fixed(u, v));
  if (this->g.has_edge(u, v))
    return {};

  vector<pair<Node,Node>> to_remark = unmark_to_modify(u, v); 
  assert(!is_marked(u,v));

  
  this->cost += 1;
  
  this->degres[u]++;
  this->degres[v]++;
  for (int w = 0; w < this->n; w++)
  {
    if (g.has_edge(w, v))
    {
      this->twins[u][w]--;
      this->twins[w][u]--;
    }
    else
    {
      this->twins[u][w]++;
      this->twins[w][u]++;
    }
    if (g.has_edge(w, u))
    {
      this->twins[v][w]--;
      this->twins[w][v]--;
    }
    else
    {
      this->twins[v][w]++;
      this->twins[w][v]++;
    }
  }
  this->g.add_edge(u, v);
  //cout << "Add " << u << "," << v << endl;
  
  ClusterId cu = this->node_to_cluster[u];
  ClusterId cv = this->node_to_cluster[v];
  this->cluster_edges_count[min(cu, cv)][max(cu, cv)][1]++; // Add a non-marked edge
  this->cluster_edges_count[min(cu, cv)][max(cu, cv)][3]--; // Delete a non-marked non-edge

  return to_remark;
  //return{};
}

vector<pair<Node,Node>> ClusterGraph::remove_graph_edge(Node u, Node v)
{
  assert(!is_fixed(u, v));
  if (!this->g.has_edge(u, v))
    return {};

  vector<pair<Node, Node>> to_remark = unmark_to_modify(u, v);
  assert(!is_marked(u,v));

  this->cost += 1;

  this->degres[u]--;
  this->degres[v]--;
  for (int w = 0; w < this->n; w++)
  {
    if (g.has_edge(w, v))
    {
      this->twins[u][w]++;
      this->twins[w][u]++;
    }
    else
    {
      this->twins[u][w]--;
      this->twins[w][u]--;
    }
    if (g.has_edge(w, u))
    {
      this->twins[v][w]++;
      this->twins[w][v]++;
    }
    else
    {
      this->twins[v][w]--;
      this->twins[w][v]--;
    }
  }
  
  this->g.remove_edge(u, v);
  //cout << "Delete " << u << "," << v << endl;

  ClusterId cu = this->node_to_cluster[u];
  ClusterId cv = this->node_to_cluster[v];
  this->cluster_edges_count[min(cu, cv)][max(cu, cv)][1]--; // Delete a non-marked edge
  this->cluster_edges_count[min(cu, cv)][max(cu, cv)][3]++; // Add a non-marked non-edge

  return to_remark;
  //return {};
}

bool ClusterGraph::is_fixed(Node u, Node v) const
{
  assert(u != v);
  ClusterId cu = this->cluster_Id_of_node(u);
  ClusterId cv = this->cluster_Id_of_node(v);
  //return (cu == cv || !this->has_edge(cu, cv));
  return (cu == cv || !cluster_adjacency[cu][cv]);
}

void ClusterGraph::display_cluster(ClusterId i) const
{
  for (auto &it : this->clusters[i])
    cout << it << " ";
  cout << endl;
}

void ClusterGraph::display_cluster_adjacency(ClusterId cu) const
{
  cout << "- ";
  vector<ClusterId> c_neighbors = neighboring_clusters(cu);
  for (auto &it : c_neighbors)
    cout << it << " ";
  cout << endl;
}

bool ClusterGraph::is_finished() const
{
  return (this->nb_cluster_edges == 0);
}

Solution ClusterGraph::getSolution()
{
  return Solution(this->cost, this->node_to_cluster);
}

// Stars related methods

int ClusterGraph::get_star_value() const
{
  return this->star_value;
}

bool ClusterGraph::is_marked(Node u, Node v) const
{
  return !this->marked_graph[u][v].is_empty() && !this->marked_graph[v][u].is_empty();
}

Star ClusterGraph::get_star_at(Node u, Node v) const
{
  assert(this->is_marked(u, v));
  return this->marked_graph[u][v];
}

vector<Star> ClusterGraph::all_stars() const
{
  vector<Star> stars(0);
  for (Node u = 0; u < this->n; u++)
  {
    for (Node v = u + 1; v < this->n; v++)
    {
      if (this->is_marked(u, v))
      {
        Star s = marked_graph[u][v];
        vector<Node> sorted_leaves = s.get_leaves();
        sort(sorted_leaves.begin(), sorted_leaves.end());
        Star sorted_s = Star(s.get_center(), sorted_leaves);
        if (find(stars.begin(), stars.end(), sorted_s) == stars.end()) 
        {
          stars.push_back(sorted_s);
        }
      }
    }
  }
  return stars;
}

vector<pair<Node, Node>> ClusterGraph::unmark_to_modify(Node u, Node v) {

  if (!is_marked(u, v))
    return {};

  Star t = marked_graph[u][v];
  vector<Node> new_leaves(0);

  for (Node l: t.get_leaves())
  {
    if (l != u  && l != v)
      new_leaves.push_back(l);
  }

  unmark_star(t);

  if (new_leaves.size() >= 2) 
  {
    Star new_star = Star(t.get_center(), new_leaves);
    mark_star(new_star);
  }

  return t.all_pairs();
}

void ClusterGraph::mark(Node u, Node v, const Star &s)
{
  assert(!this->is_marked(u, v));
  if (this->is_fixed(u, v))
    return;
  this->marked_graph[u][v] = s;
  this->marked_graph[v][u] = s; // TODO: use only half of the matrix has for cluster_edges_count?
  ClusterId cu = this->cluster_Id_of_node(u);
  ClusterId cv = this->cluster_Id_of_node(v);
  ClusterId min_c = min(cu, cv);
  ClusterId max_c = max(cu, cv);

  //if (cu == cv) return;

  if (this->g.has_edge(u, v))
  {
    this->cluster_edges_count[min_c][max_c][0] += 1;
    this->cluster_edges_count[min_c][max_c][1] -= 1;
  }
  else
  {
    this->cluster_edges_count[min_c][max_c][2] += 1;
    this->cluster_edges_count[min_c][max_c][3] -= 1;
  }
}

void ClusterGraph::unmark(Node u, Node v)
{
  if (!this->is_marked(u, v))
    return;
  this->marked_graph[u][v] = Star(); // TODO: Maybe create a unique empty star..
  this->marked_graph[v][u] = Star();
  ClusterId cu = this->cluster_Id_of_node(u);
  ClusterId cv = this->cluster_Id_of_node(v);
  ClusterId min_c = min(cu, cv);
  ClusterId max_c = max(cu, cv);

  //if (cu == cv) return;

  if (this->g.has_edge(u, v))
  {
    this->cluster_edges_count[min_c][max_c][0] -= 1;
    this->cluster_edges_count[min_c][max_c][1] += 1;
  }
  else
  {
    this->cluster_edges_count[min_c][max_c][2] -= 1;
    this->cluster_edges_count[min_c][max_c][3] += 1;
  }
}

void ClusterGraph::mark_star(const Star &s)
{
  vector<Node> leaves = s.get_leaves();
  Node center = s.get_center();
  this->star_value += leaves.size() - 1;
  for (auto pairs: s.all_pairs())
    mark(pairs.first, pairs.second, s);
}

void ClusterGraph::unmark_star(const Star &s)
{
  vector<Node> leaves = s.get_leaves();
  Node center = s.get_center();
  this->star_value -= (leaves.size() - 1);
  for (unsigned int i = 0; i < leaves.size(); i++)
  {
    this->unmark(center, leaves[i]);
    for (unsigned int j = i + 1; j < leaves.size(); j++)
      this->unmark(leaves[i], leaves[j]);
  }
}

// Auxiliary methods
vector<int> intersection(vector<int> &nums1, vector<int> &nums2)
{
  unordered_set<int> m(nums1.begin(), nums1.end());
  vector<int> res;
  for (auto a : nums2)
  {
    if (m.count(a))
    {
      res.push_back(a);
      m.erase(a);
    }
  }
  return res;
}

// Removes from a element of a that are also in b
vector<int> set_minus(vector<int> a, vector<int> b)
{
  vector<int> a_copy = vector<int>(a);
  sort(b.begin(), b.end()); //TODO: sort before
  a_copy.erase(remove_if(a_copy.begin(),
                         a_copy.end(),
                         [&](int x) { return binary_search(begin(b), end(b), x); }),
               a_copy.end());
  return a_copy;
}

pair<int, int> ClusterGraph::move_cost(Node u, ClusterId i)
{
  if (u == -1 && i == -1)
    return make_pair(-1, -1);
  ClusterId j = cluster_Id_of_node(u);
  assert(j != i);
  int before = 1 - cluster(j).size();
  for (Node v : cluster(j))
  {
    if (v != u && g.has_edge(u, v))
      before += 2;
  }
  int after = -cluster(i).size();
  for (Node v : cluster(i))
  {
    if (g.has_edge(u, v))
      after += 2;
  }
  return make_pair(after - before, cluster(i).size() - cluster(j).size());
}

void ClusterGraph::move(Node u, ClusterId i)
{
  // il faudrait gérer les degrés et twins, mais pas besoin : cette méthode n'est appelée que pour l'heuristique
  ClusterId old = cluster_Id_of_node(u);
  if (old != i)
  {
    node_to_cluster[u] = i;
    clusters[i].push_back(u);

    // clusters[old].remove(clusters[old].begin(),clusters[old].end(),u);
    // Replaced remove by find + erase
    auto pos_u = find(clusters[old].begin(), clusters[old].end(), u);
    clusters[old].erase(pos_u);

    for (Node v = 0; v < g.nr_vertices(); v++)
    {
      if (v == u)
        continue;
      ClusterId cv = cluster_Id_of_node(v);
      if (cv != i && cv != old)
      {
        int mo = min(old, cv);
        int Mo = max(old, cv);
        int mi = min(i, cv);
        int Mi = max(i, cv);

        if (g.has_edge(u, v))
        {
          if (is_marked(u, v))
          {
            cluster_edges_count[mo][Mo][0]--;
            cluster_edges_count[mi][Mi][0]++;
          }
          else
          {
            cluster_edges_count[mo][Mo][1]--;
            cluster_edges_count[mi][Mi][1]++;
          }
        }
        else
        {
          if (is_marked(u, v))
          {
            cluster_edges_count[mo][Mo][2]--;
            cluster_edges_count[mi][Mi][2]++;
          }
          else
          {
            cluster_edges_count[mo][Mo][3]--;
            cluster_edges_count[mi][Mi][3]++;
          }
        }
      }
      if (cv == old)
      {
        int m = min(old, i);
        int M = max(old, i);
        if (g.has_edge(u, v))
        {
          if (is_marked(u, v))
          {
            cluster_edges_count[m][M][0]++;
          }
          else
          {
            cluster_edges_count[m][M][1]++;
          }
        }
        else
        {
          if (is_marked(u, v))
          {
            cluster_edges_count[m][M][2]++;
          }
          else
          {
            cluster_edges_count[m][M][3]++;
          }
        }
      }
      if (cv == i)
      {
        int m = min(old, i);
        int M = max(old, i);
        if (g.has_edge(u, v))
        {
          if (is_marked(u, v))
          {
            cluster_edges_count[m][M][0]--;
          }
          else
          {
            cluster_edges_count[m][M][1]--;
          }
        }
        else
        {
          if (is_marked(u, v))
          {
            cluster_edges_count[m][M][2]--;
          }
          else
          {
            cluster_edges_count[m][M][3]--;
          }
        }
      }
    }
  }
}

vector<Cluster> ClusterGraph::get_all_clusters() const
{
  return this->clusters;
}

Solution ClusterGraph::get_solution(int c)
{
  return Solution(c+cost, node_to_cluster);
}

void ClusterGraph::clear_stars() {
  for (Node u = 0; u < n; u++) {
    for (Node v = u+1; v < n; v++) 
      if (is_marked(u, v)) 
        unmark_star(marked_graph[u][v]);
  }
}

bool ClusterGraph::has_graph_edge(Node u, Node v) const {
  return this->g.has_edge(u,v);
}

vector<pair<Node, Node>> ClusterGraph::try_mark(Node u, Node v, bool force_try, bool remark) 
{
  if (!is_feasible())
    return {};
  if (!force_try && is_marked(u, v))
    return {};
  
  vector<pair<Node, Node>> to_remark(0);

  if (is_fixed(u, v))
  {
    if (g.has_edge(u, v)) 
    {
      vector<pair<Node, Node>> extension = try_mark_split_star(u, v);
      to_remark.insert(to_remark.end(), extension.begin(), extension.end());
    }
  }
  else 
  {
    vector<pair<Node, Node>> extension_1 = try_mark_triple(u, v, true);
    to_remark.insert(to_remark.end(), extension_1.begin(), extension_1.end());

    if (g.has_edge(u, v))
    {
      vector<pair<Node, Node>> extension_2 = try_mark_merge_star(u, v);
      to_remark.insert(to_remark.end(), extension_2.begin(), extension_2.end());
    }
  }
  vector<pair<Node, Node>> extension_3 = try_mark_triple(u, v, false);
  to_remark.insert(to_remark.end(), extension_3.begin(), extension_3.end());


  
  if (!remark)
    return to_remark;
  
  for (pair<Node, Node> p: to_remark)
    try_mark(p.first, p.second, false, true);

  return {};
  
}

vector<pair<Node, Node>> ClusterGraph::try_mark_merge_star(Node a, Node b) {
  assert(!is_fixed(a, b));
  assert(g.has_edge(a, b));

  vector<pair<Node, Node>> to_remark (0);
  vector<Node> edge {a, b};
  for (int i = 0; i < 2; i++) 
  {
    Node u = edge[i];
    Node v = edge[1-i];
    for (Node x = 0; x < n; x++) 
    {
      if (x == u || marked_graph[u][x].is_empty()) 
        continue;

      Star t = marked_graph[u][x];
      if (t.get_center() != u)
        continue;

      if (t.has_leaf(v))
        continue;
      
      vector<Star> conflicting_stars (0);
      if (is_marked(a, b))
      {
        Star temp_t = marked_graph[a][b];
        bool has_marked_edge = false;
        for (Edge temp_t_edge: temp_t.edges()) 
        {
          if (is_marked(temp_t_edge.first, temp_t_edge.second))
          {
            has_marked_edge = true;
            break;
          }    
        }

        if(has_marked_edge) {
          conflicting_stars.push_back(temp_t);
        }
        else
          continue;
      }

      bool has_edge = false;

      for (Node w: t.get_leaves())
      {
        if (g.has_edge(v, w)) 
        {
          has_edge = true;
          break;
        }
        if (is_marked(v, w))
          conflicting_stars.push_back(marked_graph[v][w]);
      }

      if (has_edge)
        continue;

      if (conflicting_stars.size() > 1 || 
          (conflicting_stars.size() == 1 && conflicting_stars[0].nr_leaves() > 2))
          continue;
      
      if (conflicting_stars.size() == 1) 
      {
        Star t_del = conflicting_stars[0];
        unmark_star(t_del);
        vector<pair<Node,Node>> extension = t_del.all_pairs();
        
        int pos_uv = -1;
        for (int j = 0; j < extension.size(); j++)
          if (min(extension[j].first, extension[j].second) == min(u, v) && max(extension[j].first, extension[j].second) == max(u, v) )
            pos_uv = j;
        if (pos_uv != -1)
          extension.erase(extension.begin()+pos_uv, extension.begin() + pos_uv + 1);

        to_remark.insert(to_remark.end(), extension.begin(), extension.end());
      }

      unmark_star(t);

      vector<Node> new_leaves = t.get_leaves();
      new_leaves.push_back(v);
      Star new_t (t.get_center(), new_leaves);
      mark_star(new_t);

      return to_remark;
    }
  }

  return {};

}

vector<pair<Node, Node>> ClusterGraph::try_mark_triple(Node a, Node b, bool fixed_only) 
{
  assert(!fixed_only || !is_fixed(a, b));
  bool uv_fixed = is_fixed(a, b) && g.has_edge(a, b);
  vector<pair<Node, Node>> to_remark(0);
  vector<Node> e {a, b};

  for (int i = 0; i < 2; i++)
  {
    Node u = e[i];
    Node v = e[1-i];
    for (Node w: g.neighbours(u))
    {
      if (w  == v)
        continue;
      if ((int) g.has_edge(u, v) + (int) g.has_edge(v,w) != 1)
        continue;
      
      vector<pair<Node, Node>> marked_pair(0);
      if (is_marked(u, v)) marked_pair.push_back(pair<Node,Node>(u, v));
      if (is_marked(u, w)) marked_pair.push_back(pair<Node,Node>(u, w));
      if (is_marked(v, w)) marked_pair.push_back(pair<Node,Node>(v, w));

      if (!uv_fixed && marked_pair.size() > 0)
        continue;
      if (uv_fixed && marked_pair.size() > 1)
        continue;
      if (fixed_only && !is_fixed(u,w))
        continue;
      
      if (marked_pair.size() == 1)
      {
        pair<Node, Node> e_2 = marked_pair[0];
        Star t_del = marked_graph[e_2.first][e_2.second];

        if (t_del.nr_leaves() > 2)
          continue;
        
        vector<Node> t_del_fixed (0);
        for (Node x: t_del.get_leaves())
          if (is_fixed(t_del.get_center(), x))
            t_del_fixed.push_back(x);
        
        if (t_del_fixed.size() > 0)
          continue;
        
        unmark_star(t_del);
        vector<pair<Node, Node>> extension = t_del.all_pairs();
        to_remark.insert(to_remark.end(), extension.begin(), extension.end());
      }
    
    Node center = g.has_edge(u,v) ? u : w;
    vector<Node> leaves(0);
    vector<Node> candidate_leaves {u, v, w};
    for (Node l: candidate_leaves)
      if (l != center)
        leaves.push_back(l);
    
    Star new_t(center, leaves);
    mark_star(new_t);
    
    if (!is_fixed(a, b))
      return to_remark;
    }

  }
  return to_remark;
}

vector<pair<Node, Node>> ClusterGraph::try_mark_split_star(Node a, Node b)
{ 
  assert(g.has_edge(a, b));
  assert(is_fixed(a, b));
  vector<Node> e {a, b};
  vector<pair<Node, Node>> to_remark(0);

  for (int i = 0; i < 2; i++)
  {
    Node u = e[i];
    Node v = e[1-i];

    for (Node x = 0; x < n; x++)
    {
      if (x == u || !is_marked(u, x))
        continue;
      
      Star t = marked_graph[u][x];
      if (t.get_center() != u)
        continue;

      if (t.has_leaf(v))
        continue;
      
      vector<Node> fixed_t_edges (0);
      for (Node l: t.get_leaves())
        if (is_fixed(u, l))
          fixed_t_edges.push_back(l);
      if (t.nr_leaves() == 2 && fixed_t_edges.size() > 0)
        continue;
      
      vector<Node> can_split (0);
      for (Node l: t.get_leaves())
        if (!is_marked(v, l) && !g.has_edge(v, l) && !is_fixed(u, l))
          can_split.push_back(l);
      
      if (can_split.size() == 0)
        continue;

      unmark_star(t);
      vector<pair<Node, Node>> extension = t.all_leaves_pairs();
      to_remark.insert(to_remark.end(), extension.begin(), extension.end());

      if (can_split.size() < /*(unsigned long)*/ t.nr_leaves() - 1)
      {
        vector<Node> leaves = set_minus(t.get_leaves(), can_split);
        Star new_t (t.get_center(), leaves);
        for (auto pa: new_t.all_pairs())
          assert (!is_marked(pa.first, pa.second));
        mark_star(new_t);
      }
      else 
      {
        vector<Node> remaining = set_minus(t.get_leaves(), can_split);
        if (remaining.size() == 1)
          to_remark.push_back(pair<Node, Node>(u, remaining[0]));
      }

      for (Node l: can_split)
      {
        Star new_s (u, {v, l}); 
        mark_star(new_s);
        break;
      }

    }

  }

  return to_remark;
}

bool ClusterGraph::is_feasible() const
{
  return cost + star_value < upper_bound;
}

void ClusterGraph::set_upper_bound(int ub)
{
  upper_bound = ub;
}

void ClusterGraph::init_lower_bound() 
{
  vector<pair<Node,Node>> to_remark = g.all_edges();
  while(to_remark.size() > 0)
  {
    vector<pair<Node, Node>> to_remark_again(0);
    for (pair<Node,Node> e: to_remark)
    {
      vector<pair<Node, Node>> extension = try_mark(e.first, e.second, true, false);
      to_remark_again.insert(to_remark_again.end(), extension.begin(), extension.end());
    }
    to_remark = to_remark_again;
  }
}

bool ClusterGraph::check_packing() 
{
  int n = this->n;
  vector<vector<Star>> marked_copy(n);
  for (Node u = 0; u < n; u++)
  {
    marked_copy[u] = vector<Star> (n);
    for (Node v = u+1; v < n; v++)
    {
      marked_copy[u][v] = Star();
      if (is_marked(u,v))
      {
        Node center = marked_graph[u][v].get_center();
        vector<Node> leaves = marked_graph[u][v].get_leaves();
        Star new_s = Star(center, leaves);
        marked_copy[u][v] = new_s;
      }
    }
  }

  for (Node u = 0; u < n; u++) 
  {
    for (Node v = u+1; v < n; v++) 
    {
      Star s_uv = marked_copy[u][v];
      vector<pair<Node,Node>> allpairs = s_uv.all_pairs();
      for (unsigned int i = 1; i < allpairs.size(); i++) 
      {
        Node min_n = min(allpairs[i].first, allpairs[i].second);
        Node max_n = max(allpairs[i].first, allpairs[i].second);
        marked_copy[min_n][max_n] = Star();
      }
    }
  }
  vector<pair<Node,Node>> seen (0);

  for (Node u = 0; u < n; u++)
  {
    for (Node v = u+1; v < n; v++)
    {
      if (!marked_copy[u][v].is_empty())
      { 
        for (auto p: marked_copy[u][v].all_pairs())
        {
          Node q = min(p.first, p.second);
          Node r = max(p.first, p.second);
          if (is_fixed(q,r))
            continue;
          for (unsigned int i = 0; i < seen.size(); i++)
          {
            if (seen[i].first == q && seen[i].second == r)
            {
              cout << "ERROR : " << q << "," << r << endl;
              return false;
            }
          }   
          seen.push_back(pair<Node,Node>(q,r));         
        }
      }
    }
  }

  return true;
    
}

const vector<Node> ClusterGraph::graph_neighbors(Node u) const
{
  return g.neighbours(u);
}

int ClusterGraph::nr_graph_vertices() const
{
  return n;
}

vector<Edge> ClusterGraph::all_graph_edges() const
{
  return g.all_edges();
}

void ClusterGraph::shuffle_graph() 
{
  this->g.shuffle();
}

void ClusterGraph::try_kernelize_lonely_vertex(Node u)
{
  bool is_reducible = false;
  ClusterId cu = cluster_Id_of_node(u);

  for (Node v: graph_neighbors(u))
    if (cluster_Id_of_node(v) != cu)
      return;

  int nr_edges_to_cu = 0;
  for (ClusterId cv  = 0; cv < n; cv++)
  {
    if (cv == cu || clusters[cv].size() == 0)
      continue;
    else
      nr_edges_to_cu += nr_edges(cu, cv);
  }

  if (nr_edges_to_cu < clusters[cu].size()) // Vertex is reducible
  {
    for (Node v = 0;  v < n; v++)
    {
      if (v == u || cluster_Id_of_node(v) == cu)
        continue;
      int nr_neighbor_in_cu = 0;

      for (Node w: graph_neighbors(v))
        if (cluster_Id_of_node(w) == cu)
          nr_neighbor_in_cu++;

      if (nr_neighbor_in_cu <= clusters[cu].size()/2 && cluster_adjacency[cu][cluster_Id_of_node(v)])
      {
        split(cluster_Id_of_node(v), cu, true);
        //cout << "kernerlize" << endl;
      }
    }
  }
  
}