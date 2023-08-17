#include "ClusterVariants.h"

ClusterVariants::ClusterVariants(string title,
				 double p_threshold,
				 vector< ClusterVariantItem >& vars,
				 int nlibs,
				 int thread_count)
  : m_title(title), m_p_threshold(p_threshold), m_vars(vars)
{
  m_nvars = m_vars.size();
  init_matrix(thread_count);
};

void dump_vec(string title, vector<double> v)
{  
  cout << title << "=(" << v[0]; for (unsigned int i=1; i<v.size(); ++i) cout << "," << v[i]; cout  << ") ";
}  

void ClusterVariants::init_matrix(int index, int from_ind, int to_ind)
{
  // cout << "init matrix: " << m_title << endl;
  for (int i1=from_ind; i1<to_ind; ++i1) {
    vector<bool> rr(m_nvars);
    for (int i2=0; i2<m_nvars; ++i2) {
      massert(i1 < (int)m_selected.size() && i2 < (int)m_selected[i1].size(), "matrix out of range");
      double p_value = chi_square_test_p(m_vars[i1].counts, m_vars[i2].counts);
      bool selected = p_value > m_p_threshold;
      rr[i2] = selected;
    }
    // work with mtx to protect shared data structure
    {
      std::lock_guard<std::mutex> lk(m_mtx);
      m_selected[i1] = rr;
    }
  }
  // cout << "init matrix done: " << m_title << endl;
}

void cluster_variants_init_matrix(ClusterVariants* cv, int index, int from_ind, int to_ind)
{
  cv->init_matrix(index, from_ind, to_ind);
}

void ClusterVariants::init_matrix(int thread_count)
{
  m_selected.resize(m_nvars);
  for (int i=0; i<m_nvars; ++i)
    m_selected[i].resize(m_nvars);
  
  if (thread_count > m_nvars)
    thread_count = m_nvars;

  int step = floor(m_nvars / thread_count);
  // cout << "number of threads used to initialize distance matrix: " << thread_count << endl;

  vector<thread> threads;
  for (int i=0; i<thread_count; ++i) {
    int from_ind = i * step;
    int to_ind = (i < (thread_count-1)) ? (i+1) * step : m_nvars;
    // cout << "init matrix, index=" << i << ", from_index=" << from_ind << ", to_index=" << to_ind << endl;
    threads.push_back(thread(cluster_variants_init_matrix, this, i, from_ind, to_ind));
  }

  // cout << "waiting for all threads to finish\n";
  for (auto& th : threads) th.join();
}

int ClusterVariants::cluster_variants(vector<int>& vclusters)
{
  vclusters.resize(m_nvars);

  using namespace boost;
  {
    typedef adjacency_list <vecS, vecS, undirectedS> Graph;
    Graph G;

    // vertices
    for (int i=0; i<m_nvars; ++i)
      add_vertex(G);

    int edge_count = 0;
    for (int i1=0; i1<m_nvars; ++i1) {
      for (int i2=i1; i2<m_nvars; ++i2) {
	if (m_selected[i1][i2]) {
	  edge_count++;
	  add_edge(i1, i2, G);
	}
      }
    }
    // massert(edge_count > 0, "no pairs of variants were associated, consider reducing p threshold");

    // cout << "connected_components: " << m_title << endl;
    int num = connected_components(G, &vclusters[0]);
    // cout << "connected_components done: " << m_title << endl;

    return num;
  }
}
