#include "instance.hpp"
#include <assert.h>
#include <unordered_set>

struct star_hasher {
  
  std::size_t operator()(const Star& s) const {
      std::size_t seed = s.nr_leaves()+1;
      seed ^= s.get_center() + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      //std::size_t seed = s.get_center();
      for(auto& i : s.get_leaves()) 
      {
        seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      }
      return seed;
  }

};

bool is_packing(vector<Star> &p);

extern vector<int> intersection(vector<int> &nums1, vector<int> &nums2);
extern vector<int> set_minus(vector<int> a, vector<int> b);

ExactInstance::ExactInstance(Graph g, int upper_bound, int cost)
{
  //this->g = new Graph(g);
  this->cluster_graph = ClusterGraph(g, cost, upper_bound);
  this->upper_bound = upper_bound;
  this->nb_branch = 0;
  this->nr_graph_vertices = g.nr_vertices();
}

ExactInstance::ExactInstance(const ExactInstance &old_instance)
{
  //this->g = new Graph(*(old_instance.g));
  this->cluster_graph = ClusterGraph(old_instance.cluster_graph); //old_instance.g);
  this->upper_bound = old_instance.upper_bound;
  this->nb_branch = old_instance.nb_branch;
  this->nr_graph_vertices = old_instance.nr_graph_vertices;
}


int ExactInstance::get_upper_bound()
{
  return this->upper_bound;
}

int ExactInstance::get_size()
{
  return this->nr_graph_vertices;
}

int ExactInstance::get_cost()
{
  return this->cluster_graph.get_cost();
}

int ExactInstance::degree(Node u)
{
  return this->cluster_graph.degree(u);
}

int ExactInstance::twinness(Node u, Node v)
{
  return this->cluster_graph.twinness(u, v);
}

void ExactInstance::set_upper_bound(int value)
{
  this->upper_bound = value;
  this->cluster_graph.set_upper_bound(value);
}

void ExactInstance::preprocess()
{
  int n = get_size();
  for (int u = 0; u < n; u++)
  {
    if (degree(u) == 1)
    {
      Node v = graph_neighbors(u)[0];
      if (degree(v) <= 2)
      {
        if (get_cluster(v).size() > 1)
          continue;
        fix_pair(u, v);
        
        ClusterId cv = get_cluster_Id(v);
        for (Node w : graph_neighbors(v))
        {
          if (w != u)
          {
            split(get_cluster_Id(w), cv, true);
            
          }
        }
      }
    }

    for (int v = 0; v < n; v++)
    {
      if (u >= v || is_fixed(u, v))
      {
        continue;
      }
      if (cluster_graph.false_twins(u, v))
      {
	      fix_pair(u,v);
      }
      if (twinness(u, v) <= 1)
      {
        fix_pair(u, v);
      }
    }
  }


  std::vector<std::vector<int>> tw2;
  for (int u = 0; u < n; u++)
  {
    std::vector<int> tw2u;
    for (int v = 0; v < n; v++)
    {
      if (twinness(u, v) <= 2){
        tw2u.push_back(v);
      }
    }
    tw2.push_back(tw2u);
  }
  for (int u = 0; u < n; u++)
  {
    if (tw2[u].size() >= 7)
    {
      for (int v : tw2[u])
      {
        fix_pair(u, v);
      }
    }
    else
    {
      for (int v : tw2[u])
      {
        for (int w : tw2[u])
        {
          if (w > v && twinness(v, w) <= 2)
          {
            for (int x : tw2[u])
            {
              if (x > w && twinness(v, x) <= 2 && twinness(w, x) <= 2)
              {
                fix_pair(u, v);
                fix_pair(u, w);
                fix_pair(u, x);
              }
            }
          }
        }
      }
    }
  }
  for (int i = 3; i <= 8; i++)
  {
    for (int u = 0; u < n; u++)
    {
      std::vector<int> tw;
      for (int v = 0; v < n; v++)
      {
        if (u != v && twinness(u, v) <= i)
          tw.push_back(v);
      }
      if ((int)tw.size() >= 4 * i - 1)
      {
        for (int v : tw)
        {
          fix_pair(u, v);
        }
      }
    }
  }
  

}

void ExactInstance::init_lower_bound()
{
  if (nr_graph_vertices >= 600) //>=600
    greedy_complete_lower_bound();
  else if (nr_graph_vertices < 600 && nr_graph_vertices > 250) // > 250
    fast_lower_bound_triple();
  else
    lower_bound_by_triple();
  cluster_graph.init_lower_bound();
}

bool ExactInstance::init_lower_bound_triple()
{
  return true;
}

bool ExactInstance:: init_lower_bound_greedy()
{
  return true;
}

void ExactInstance::recompute_lower_bound()
{
  ExactInstance temp(*this);
  temp.clear_stars();
  temp.init_lower_bound();
  if (temp.lower_bound() > this->lower_bound())
  {
    clear_stars();
    for (Star t : temp.all_stars()) {
      mark_star(t);
    }
  }
}

void ExactInstance::fast_init_lower_bound()
{
  fast_lower_bound_triple();
  //greedy_complete_lower_bound();
  cluster_graph.init_lower_bound(); 
}

void ExactInstance::fast_recompute_lower_bound() 
{
  ExactInstance temp(*this);
  temp.clear_stars();
  temp.fast_init_lower_bound();
  if (temp.lower_bound() > this->lower_bound())
  {
    clear_stars();
    for (Star t : temp.all_stars()) {
      mark_star(t);
    }
  }
}

void ExactInstance::forced_moves()
{
  
  int diff = this->upper_bound - this->lower_bound();
  
  
  /*if (nb_branch % 30 == 0)
  {
    //cout << "Recompute " << lower_bound() << endl;
    fast_recompute_lower_bound();
    //cout << "After " << lower_bound() << endl;
    // TODO: include that in the algo?
  }*/
  

  bool modified = true;
  while (modified && this->is_feasible())
  {
    modified = false;

    for (Node u = 0; u < nr_graph_vertices; u++)
    {
      cluster_graph.try_kernelize_lonely_vertex(u);
    }

    for (ClusterId cu = 0; cu < nr_graph_vertices; cu++)
    {
      if (cluster_graph.cluster(cu).size() == 0)
        continue;
      for (ClusterId cv = cu + 1; cv < nr_graph_vertices; cv++)
      {
        if (cluster_graph.cluster(cv).size() == 0 || !can_merge(cu, cv))
          continue;
        if (!is_feasible())
          break;
        
        if (cluster_graph.nr_edges(cu, cv, false) >= upper_bound - lower_bound())
        {
          merge(cu, cv, true);
          modified = true;
          continue;
        }
        
        if (cluster_graph.nr_non_edges(cu, cv, false) >= upper_bound - lower_bound())
        {
          split(cu, cv, true);
          modified = true;
          continue;
        }
        
        
        for (Node u : cluster_graph.cluster(cu))
        {
          int dbaru = 0;
          for (Node v : cluster_graph.cluster(cv))
          {
            if (!has_graph_edge(u, v))
              dbaru++;
          }
          if (dbaru > degree(u))
          {
            split(cu, cv, true);
            modified = true;
            break;
          }
        }
        

        if (!can_merge(cu, cv))
          continue;
        
        for (Node v : cluster_graph.cluster(cv))
        {
          int dbarv = 0;
          for (Node u : cluster_graph.cluster(cu))
          {
            if (!has_graph_edge(u, v))
              dbarv++;
          }
          if (dbarv > degree(v))
          {
            split(cu, cv, true);
            modified = true;
            break;
          }
        }

        if (!is_feasible() || !can_merge(cu, cv))
          continue;
        
        int augment_lb = 0;
        
        for (Node u2 = 0; u2 < nr_vertices(); u2++)
        {
          ClusterId cu2 = get_cluster_Id(cu);
          if (cu2 == cu || cu2 == cv)
            continue;
          int nr_non_marked_edges_cu = 0;
          int nr_non_marked_non_edges_cu = 0;
          int nr_non_marked_edges_cv = 0;
          int nr_non_marked_non_edges_cv = 0;
          for (Node u: cluster_graph.cluster(cu))
          {
            if (is_marked(u, u2))
              continue;
            if (has_graph_edge(u, u2))
              nr_non_marked_edges_cu++;
            else
              nr_non_marked_non_edges_cu++;
          }
          for (Node v: cluster_graph.cluster(cv))
          {
            if (is_marked(v, u2))
              continue;
            if (has_graph_edge(v, u2))
              nr_non_marked_edges_cv++;
            else
              nr_non_marked_non_edges_cv++;
          }
          augment_lb += min(nr_non_marked_edges_cu, nr_non_marked_non_edges_cv)
            + min(nr_non_marked_edges_cv, nr_non_marked_non_edges_cu);
        }

        int merge_lb = lower_bound() + augment_lb + cluster_graph.nr_non_edges(cu, cv, false);
        
        for (ClusterId cu2: cluster_graph.neighboring_clusters(cu))
        {
          if (cu2 == cv || cluster_graph.has_edge(cu2, cv))
            continue;
          merge_lb += cluster_graph.nr_edges(cu2, cu, false);
        }
        
        for (ClusterId cv2: cluster_graph.neighboring_clusters(cv))
        
        {
          if (cv2 == cu || cluster_graph.has_edge(cv2, cu))
            continue;
          merge_lb += cluster_graph.nr_edges(cv2, cv, false);
        }
        
        if (merge_lb >= upper_bound)
        {
          split(cu, cv, true);
          modified = true;
          continue;
        }
        
        
        /*
        std::vector<int> ncu;
        for (Node u : cluster_graph.cluster(cu))
        {


          for (Node u2 : graph_neighbors(u))
          {

            if (!is_marked(u,u2) && !cluster_graph.cluster_Id_of_node(u2) == cv)
            //if (!is_marked(u, u2) && cluster_graph.cluster_Id_of_node(u2) != cv)
            if (is_marked(u, u2) || cluster_graph.cluster_Id_of_node(u2) == cv)
              continue;
            bool keep = false;
            
            for (Node v : cluster_graph.cluster(cv))
            {
              if (!is_marked(u2, v) && !has_graph_edge(u2, v))
              {
                keep = true;
                break;
              }
            }
            if (keep)
              ncu.push_back(u2);
          }
        }
        std::vector<int> ncv;
        for (Node v : cluster_graph.cluster(cv))
        {
          for (Node v2 : graph_neighbors(v))
          {
            //if (!is_marked(v, v2) && cluster_graph.cluster_Id_of_node(v2) != cu)
            if (is_marked(v, v2) || cluster_graph.cluster_Id_of_node(v2) == cu)
              continue;
            bool keep = false;
            for (Node u : cluster_graph.cluster(cu))
            {
              if (!is_marked(v2, u) && !has_graph_edge(v2, u))
              {
                keep = true;
                break;
              }
            }
            if (keep)
              ncv.push_back(v2);
          }
        }
        
        
        int merge_lb = lower_bound() + ncu.size() + ncv.size() + cluster_graph.nr_non_edges(cu, cv, false);
        
        //for (ClusterId cu2 = 0; cu2 < nr_graph_vertices(); cu2++)
        for (ClusterId cu2: cluster_graph.neighboring_clusters(cu))
        {
          if (cu2 == cv || cluster_graph.has_edge(cu2, cv))
            continue;
          merge_lb += cluster_graph.nr_edges(cu2, cu, false);
        }
        
        for (ClusterId cv2: cluster_graph.neighboring_clusters(cv))
        //for (ClusterId cv2 = 0; cv2 < nr_graph_vertices(); cv2++)
        {
          if (cv2 == cu || cluster_graph.has_edge(cv2, cu))
            continue;
          merge_lb += cluster_graph.nr_edges(cv2, cv, false);
        }
        
        if (merge_lb >= upper_bound)
        {
          split(cu, cv, true);
          modified = true;
          continue;
        }
        */
      }
    }
  }
  
  
  
  
  for (ClusterId cu = 0; cu < nr_graph_vertices; cu++)
  {
    if (cluster_graph.cluster(cu).size() == 0 || cluster_graph.neighboring_clusters(cu).size() == 0)
      continue;
    int max_diff = 0;
    
    //for (ClusterId cv = 0; cv < nr_graph_vertices(); cv++)
    for (ClusterId cv: cluster_graph.neighboring_clusters(cu))
    {
      //if (cu == cv || !cluster_graph.has_edge(cu, cv))
        //continue;
      max_diff = max(max_diff, cluster_graph.nr_edges(cu, cv) - cluster_graph.nr_non_edges(cu, cv));
    }
    if (max_diff <= 0) // TODO: check with < 
    {
      isolate_cluster(cu);
    }
  }
  

  
  
}

const Cluster ExactInstance::get_cluster(Node u) const
{
  return this->cluster_graph.cluster(this->get_cluster_Id(u));
}

ClusterId ExactInstance::get_cluster_Id(Node u) const
{
  return this->cluster_graph.cluster_Id_of_node(u);
}

vector<pair<Node,Node>> ExactInstance::merge(ClusterId c1, ClusterId c2, bool remark)
{
  return this->cluster_graph.merge(c1, c2, remark);
}

vector<pair<Node,Node>> ExactInstance::split(ClusterId c1, ClusterId c2, bool remark)
{
  return this->cluster_graph.split(c1, c2, remark);
}

void ExactInstance::isolate_cluster(ClusterId c)
{
  for (ClusterId c2 = 0; c2 < nr_graph_vertices; c2++)
  {
    if (c2 != c)
      split(c, c2, true);
  }
}

bool ExactInstance::can_merge(ClusterId c1, ClusterId c2) const
{
  // ASSERT c1 != c2
  return this->cluster_graph.has_edge(c1, c2);
}

bool ExactInstance::is_fixed(Node u, Node v) const
{
  // ASSERT u != v
  ClusterId cu = this->get_cluster_Id(u);
  ClusterId cv = this->get_cluster_Id(v);
  return (cu == cv || !this->cluster_graph.has_edge(cu, cv));
}

void ExactInstance::fix_pair(Node u, Node v)
{
  if (this->is_fixed(u, v))
    return;
  ClusterId cu = this->get_cluster_Id(u);
  ClusterId cv = this->get_cluster_Id(v);
  if (this->has_graph_edge(u, v))
    this->merge(cu, cv, true);
  else
    this->split(cu, cv, true);
}

void ExactInstance::display_cluster_adjacency(ClusterId cu) const
{
  this->cluster_graph.display_cluster_adjacency(cu);
}

void ExactInstance::display_cluster(ClusterId cu) const
{
  this->cluster_graph.display_cluster(cu);
}

void ExactInstance::display_graph_adjacency() const
{
  //this->g->display_adjacency();
}

int ExactInstance::lower_bound() const
{
  return this->cluster_graph.get_star_value() + this->cluster_graph.get_cost();
}

bool ExactInstance::is_feasible() const
{
  return this->lower_bound() < this->upper_bound;
}


Solution ExactInstance::algo_brute()
{
  if (nb_branch % 10 == 0)
    preprocess();
  nb_branch++;

  if (!is_feasible())
  {
    return NO_SOLUTION();
  }


  if (cluster_graph.is_finished())
  {
    return cluster_graph.getSolution();
  }

  Solution best = NO_SOLUTION();

  std::pair<ClusterId, ClusterId> c = cluster_graph.bestpair();

  ExactInstance c_merge(*this);

  c_merge.merge(c.first, c.second, true);

  if (nb_branch % 2 == 0)
    c_merge.forced_moves();

  if (c_merge.is_feasible())
  {
    c_merge.set_upper_bound(std::min(c_merge.upper_bound, this->upper_bound));
    best = c_merge.algo_brute();
  }

  ExactInstance c_split(*this);
  
  c_split.split(c.first, c.second, true);

  if (nb_branch % 2  == 0)
    c_split.forced_moves();

  if (c_split.is_feasible())
  {
    c_split.set_upper_bound(std::min(c_split.upper_bound, this->upper_bound));
    Solution res = c_split.algo_brute();
    if (res.get_cost() < this->upper_bound)
    {
      upper_bound = res.get_cost();
      if (res.get_cost() < best.get_cost()) 
      { 
        best = res;
      }
    }
  }

  return best;
}

bool ExactInstance::has_edge(Node u, Node v)
{
  return this->has_graph_edge(u, v);
}

std::vector<Node> ExactInstance::neighbours(Node u)
{
  return this->graph_neighbors(u);
}

int ExactInstance::nr_vertices()
{
  return this->nr_graph_vertices;
}

Solution ExactInstance::heuristic()
{
  int n = nr_graph_vertices;
  if (n == 1)
  {
    std::vector<int> v(1, 0);
    return Solution((int64_t)0, v);
  }

  int current = 0;
  for (Node i = 0; i < n; i++) {
    for (Node j = i+1; j < n; j++) {
      if (get_cluster_Id(i) != get_cluster_Id(j) && has_edge(i,j)) current++;
    }
  }
  while (true)
  {
    int k = -1;
    int clk = -1;
    for (Node i = 0; i < n; i++)
    {
      for (ClusterId cli = 0; cli < n; cli++)
      {
        if (cluster_graph.cluster_Id_of_node(i) != cli && cluster_graph.move_cost(i, cli) > cluster_graph.move_cost(k, clk))
        {
          k = i;
          clk = cli;
        }
      }
    }

    std::pair<int, int> m = cluster_graph.move_cost(k, clk);

    if (m.first < 0 || (m.first == 0 && m.second < 0))
      break;
    current -= m.first;
    cluster_graph.move(k, clk);
  }
  return cluster_graph.get_solution(current);
}

bool ExactInstance::is_marked(Node u, Node v) const
{
  return this->cluster_graph.is_marked(u, v);
}

void ExactInstance::mark_star(Star s) {
  this->cluster_graph.mark_star(s);
}

void ExactInstance::unmark_star(Star s) {
  this->cluster_graph.unmark_star(s);
}

const Star ExactInstance::get_star(Node u, Node v) const {
  return this->cluster_graph.get_star_at(u, v);
}

const vector<Star> ExactInstance::all_stars() const {
  return this->cluster_graph.all_stars();
}


void ExactInstance::lower_bound_by_triple() {
  
  assert (this->all_stars().size() == 0);
  int n = this->nr_vertices();
  vector<vector<unordered_set<Star, star_hasher>>> triples_for (n);
  //list<Star> triples (0);

  std::unordered_set<Star, star_hasher> triples;

  for (Node u = 0; u < n; u++) {
    triples_for[u] = vector<unordered_set<Star, star_hasher>> (n);
    for (Node v = u+1; v < n; v++)
    triples_for[u][v] = unordered_set<Star, star_hasher> (0);
  }
  
  for (Node u = 0; u < n ; u++) {
    for (Node v: this->graph_neighbors(u)) {
      for (Node w: this->graph_neighbors(u)) {
        if (v >= w || this->has_graph_edge(v, w)) 
          continue;
        Star t = Star (u, {v,w});
        triples.insert(t);
        //vector<Edge> pairs = t.all_pairs();
        for (Edge e: t.all_pairs())
          triples_for[min(e.first, e.second)][max(e.first, e.second)].insert(Star(u,{v, w}));
      }
    }

  }

  
  for (Node u = 0; u < n; u++) {
    for (Node v = u+1; v < n; v++) 
      if (this->is_fixed(u, v)) {
        for (Star to_delete: triples_for[u][v]) 
        {
          triples.erase(to_delete);
        }
        triples_for[u][v].clear();
      }
  }
  
  while (triples.size() > 0) {
    Star t_min = *triples.begin();
    int deg_t_min = 0;
    for (Edge e: t_min.all_pairs()) 
      deg_t_min += triples_for[min(e.first, e.second)][max(e.first,e.second)].size();
    for (Star t: triples) {
      int deg = 0;
      for (Edge e: t.all_pairs()) 
        deg += triples_for[min(e.first, e.second)][max(e.first, e.second)].size();
      if (deg < deg_t_min) {
        t_min = t;
        deg_t_min = deg;
      }
    }

    for (Edge e: t_min.all_pairs()) {
      int u = min(e.first, e.second);
      int v = max(e.first, e.second);

      for (Star s: triples_for[u][v]) {
        triples.erase(s);
        for (Edge e_s: s.all_pairs()) {
          int u_s = min(e_s.first, e_s.second);
          int v_s = max(e_s.first, e_s.second);
          if (u_s == u  && v_s == v)
            continue;
          triples_for[u_s][v_s].erase(s);
        }
      }
      triples_for[u][v].clear();
    }

    this->mark_star(t_min);

  }
  
}

bool sortbysec(const pair<Star,int> &a, const pair<Star,int> &b)
{
    return (a.second < b.second);
}

bool ExactInstance::triple_to_claw()
{
    bool modified = false;
    int n = nr_vertices();
    for (Edge e: cluster_graph.all_graph_edges())
    {
      vector<Node> sorted_e {min(e.first, e.second), max(e.first, e.second)};
      for (int i = 0; i < 2; i++)
      {
        Node u = sorted_e[i];
        Node v = sorted_e[1-i];
        for (Node x = 0; x < n; x++)
        {
          if (x == u || !is_marked(u,x))
            continue;

          Star t = get_star(u,x);

          if (t.has_leaf(v))
            continue;
          
          if (is_marked(u,v) && get_star(u,v).nr_leaves() > 2)
            break;

          vector<Star> marked_star(0);
          if (is_marked(e.first, e.second))
            marked_star.push_back(get_star(e.first, e.second));
          
          int nb_edges = 0;

          for (Node w: t.get_leaves())
          {
            if (is_marked(v, w))
              marked_star.push_back(get_star(v, w));
            if (has_graph_edge(v,w))
              nb_edges++;
          }
          
          if (nb_edges == 0 && (marked_star.size() == 0 || (marked_star.size() == 1 && marked_star[0].nr_leaves() == 2)))
          {
            modified = true;
            for (Star t_del: marked_star)
              unmark_star(t_del);
            vector<Node> new_leaves = t.get_leaves();
            new_leaves.push_back(v);
            Star t_new (u, new_leaves);
            unmark_star(t);
            mark_star(t_new);
            break;
          }

        }

      }
    }
    return modified;
}


void ExactInstance::fast_lower_bound_triple()
{
  assert (this->all_stars().size() == 0);
  bool modified = true;

  while(modified)
  {
    modified = false;
    int n = this->nr_vertices();

    vector<pair<Star, int>> triples (0);
    vector<vector<int>> nr_pair_used(n);
    for (int i = 0; i < n; i++)
      nr_pair_used[i] = vector<int>(n, 0);
    
    for (Node u = 0; u < n ; u++) {
      for (Node v: this->graph_neighbors(u)) {
        for (Node w: this->graph_neighbors(u)) {
          if (v >= w || this->has_graph_edge(v, w)) 
            continue;
          if (is_marked(u,v) || is_marked(u,w) || is_marked(v,w))
            continue;
          Star t = Star (u, {v,w});
          triples.push_back(pair<Star,int>(t, 0));
          vector<Edge> pairs = t.all_pairs();
          for (Edge e: t.all_pairs())
            nr_pair_used[e.first][e.second]++;
        }
      }

    }

    for (pair<Star, int> p: triples)
    {
      for (Edge e: p.first.all_pairs())
        p.second += nr_pair_used[e.first][e.second];
    }

    sort(triples.begin(), triples.end(), sortbysec);
    
    for (pair<Star, int> p: triples)
    {
      Star t = p.first;
      bool to_mark = true;
      for (Edge e: t.all_pairs())
      {
        if (is_marked(e.first, e.second))
          to_mark = false;
      }
      if (to_mark)
        mark_star(t);
    }
  
    //modified = triple_to_claw();
    modified = false;
  }  
}


void ExactInstance::clear_stars() {
  cluster_graph.clear_stars();
}

int ExactInstance::get_star_value() const {
  return cluster_graph.get_star_value();
}

const Cluster ExactInstance::cluster(ClusterId i) const {
  return this->cluster_graph.cluster(i);
}

bool ExactInstance::check_packing() 
{
  int n = nr_vertices();
  vector<vector<Star>> marked_copy(n);
  for (Node u = 0; u < n; u++)
  {
    marked_copy[u] = vector<Star> (n);
    for (Node v = u+1; v < n; v++)
    {
      marked_copy[u][v] = Star();
      if (is_marked(u,v))
      {
        Node center = get_star(u,v).get_center();
        vector<Node> leaves = get_star(u,v).get_leaves();
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

void ExactInstance::greedy_complete_lower_bound() 
{
  bool modified = true;
  int n = nr_graph_vertices;
  while (modified)
  {
    modified = false;
    for (Node u = 0; u < n; u++)
    {
      for (Node v: graph_neighbors(u))
      {
        for (Node w: graph_neighbors(u))
        {
          if (v >= w || has_graph_edge(v,w))
            continue;
          if (is_marked(u,v) || is_marked(u,w) || is_marked(v,w))
            continue;
          Star new_s (u, {v,w});
          mark_star(new_s);
        }
      }
    }

    for (Edge e: cluster_graph.all_graph_edges())
    {
      vector<Node> sorted_e {min(e.first, e.second), max(e.first, e.second)};
      for (int i = 0; i < 2; i++)
      {
        Node u = sorted_e[i];
        Node v = sorted_e[1-i];
        for (Node x = 0; x < n; x++)
        {
          if (x == u || !is_marked(u,x))
            continue;

          Star t = get_star(u,x);

          if (t.has_leaf(v))
            continue;
          
          if (is_marked(u,v) && get_star(u,v).nr_leaves() > 2)
            break;

          vector<Star> marked_star(0);
          if (is_marked(e.first, e.second))
            marked_star.push_back(get_star(e.first, e.second));
          
          int nb_edges = 0;

          for (Node w: t.get_leaves())
          {
            if (is_marked(v, w))
              marked_star.push_back(get_star(v, w));
            if (has_graph_edge(v,w))
              nb_edges++;
          }
          
          if (nb_edges == 0 && (marked_star.size() == 0 || (marked_star.size() == 1 && marked_star[0].nr_leaves() == 2)))
          {
            modified = true;
            for (Star t_del: marked_star)
              unmark_star(t_del);
            vector<Node> new_leaves = t.get_leaves();
            new_leaves.push_back(v);
            Star t_new (u, new_leaves);
            unmark_star(t);
            mark_star(t_new);
            break;
          }

        }

      }
    }
  }
}

ExactInstance::~ExactInstance() {
  //delete g;
}


bool ExactInstance::has_graph_edge(Node u, Node v) const 
{
  return cluster_graph.has_graph_edge(u, v);
}

const vector<Node> ExactInstance::graph_neighbors(Node u) const
{
  return cluster_graph.graph_neighbors(u);
}

void ExactInstance::shuffle_graph()
{
  this->cluster_graph.shuffle_graph();
}