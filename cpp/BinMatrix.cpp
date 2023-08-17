#include "BinMatrix.h"

BinMatrix::BinMatrix(vector< BinSegmentSide >& segs, int nlibs,
		     double p_threshold, bool debug)
  : m_segs(segs), m_chi_squared(nlibs-1), m_debug(debug)
{
  m_nsegs = m_segs.size();
  m_selected.resize(m_nsegs);
  for (int i=0; i<m_nsegs; ++i)
    m_selected[i].resize(m_nsegs);
  m_chi_threshold = quantile(m_chi_squared, 1-p_threshold);
  cout << "Threshold: P<" << p_threshold << " CS<=" << m_chi_threshold << endl;
};

double BinMatrix::get_chi_value(vector < double >& c1, vector < double >& c2)
{
  massert(c1.size() == c2.size(), "vector length not equal");
  int n = c1.size();

  // get totals
  double tot1 = 0;
  double tot2 = 0;
  for (int i=0; i<n; ++i) {
    tot1 += c1[i];
    tot2 += c2[i];
  }
  double tot = tot1 + tot2;
  double f1 = tot1/tot;
  double f2 = tot2/tot;

  // get lib freq
  vector <double> freq_lib(n);
  for (int i=0; i<n; ++i)
    freq_lib[i] = (c1[i] + c2[i]) / tot;

  double result = 0;
  for (int i=0; i<n; ++i) {
    double exp1 = tot * f1 * freq_lib[i];
    double exp2 = tot * f2 * freq_lib[i];
    result += ((c1[i] - exp1) * (c1[i] - exp1)) / exp1;
    result += ((c2[i] - exp2) * (c2[i] - exp2)) / exp2;
  }
  return result;
}

double BinMatrix::get_p_value(vector < double >& c1, vector < double >& c2)
{
  double stat = get_chi_value(c1, c2);
  // cout << "chi=" << stat << endl;
  return (1 - cdf(m_chi_squared, stat));
}

double BinMatrix::compute_seg_chi_distance(BinSegmentSide& seg1, BinSegmentSide& seg2)
{
  return get_chi_value(seg1.counts, seg2.counts);
  // return get_p_value(seg1.counts, seg2.counts);
}

double BinMatrix::compute_coverage_factor(BinSegmentSide& seg1, BinSegmentSide& seg2)
{
  double slope = slope_origin(seg1.counts, seg2.counts);
  double length_ratio = seg1.side_length / seg2.side_length;
  return slope * length_ratio;
}

void BinMatrix::init_matrix(int index, int from_ind, int to_ind)
{
  for (int i1=from_ind; i1<to_ind; ++i1) {
    for (int i2=0; i2<m_nsegs; ++i2) {
      massert(i1 < (int)m_selected.size() && i2 < (int)m_selected[i1].size(), "matrix out of range");
      bool selected = false;
      if (m_segs[i1].seed_bin == m_segs[i2].seed_bin) {
	double chi_value = compute_seg_chi_distance(m_segs[i1], m_segs[i2]);
	selected = chi_value <= m_chi_threshold;
	if (m_debug) {
	  double P = 1 - cdf(m_chi_squared, chi_value);
	  cout << "(" << m_segs[i1].id << "," << m_segs[i2].id << ") CS=" << chi_value << " P=" << P;
	  cout << " S=" << (selected ? "T" : "F") << endl;
	}
      }
      m_selected[i1][i2] = selected;
    }
  }
}

void global_init_matrix(BinMatrix* bin_matrix, int index, int from_ind, int to_ind)
{
  bin_matrix->init_matrix(index, from_ind, to_ind);
}

void BinMatrix::init_matrix(int thread_count)
{
  if (m_debug)
    thread_count = 1;
  
  if (thread_count > m_nsegs)
    thread_count = m_nsegs;

  int step = floor(m_nsegs / thread_count);
  // cout << "number of segments: " << m_nsegs << endl;
  cout << "number of threads used to initialize distance matrix: " << thread_count << endl;

  vector<thread> threads;
  for (int i=0; i<thread_count; ++i) {
    int from_ind = i * step;
    int to_ind = (i < (thread_count-1)) ? (i+1) * step : m_nsegs;
    // cout << "init matrix, index=" << i << ", from_index=" << from_ind << ", to_index=" << to_ind << endl;
    threads.push_back(thread(global_init_matrix, this, i, from_ind, to_ind));
  }

  cout << "waiting for all threads to finish\n";
  for (auto& th : threads) th.join();
}

int BinMatrix::cluster_segments(vector<int>& bins, vector< pair<int,int> >& seg_pairs)
{
  bins.resize(m_nsegs);

  using namespace boost;
  {
    typedef adjacency_list <vecS, vecS, undirectedS> Graph;
    Graph G;

    // vertices
    for (int i=0; i<m_nsegs; ++i)
      add_vertex(G);

    // internal edges
    for (unsigned int i=0; i<seg_pairs.size(); ++i) {
      int index1 = seg_pairs[i].first;
      int index2 = seg_pairs[i].second;
      massert(m_segs[index1].id == m_segs[index2].id,
	      "expecting segment sides to be adjacent, %s != %s",
	      m_segs[index1].id.c_str(), m_segs[index2].id.c_str());
      massert(index1<m_nsegs && index2<m_nsegs, "one of the indices is out of range: i1=%d,i2=%d", index1, index2);
      add_edge(index1, index2, G);
    }

    // external edges
    int edge_count = 0;
    for (int i1=0; i1<m_nsegs; ++i1) {
      for (int i2=i1; i2<m_nsegs; ++i2) {
	if (m_segs[i1].seed_bin != m_segs[i2].seed_bin)
	  continue;
	bool selected = m_selected[i1][i2];
	if (selected) {
	  edge_count++;
	  add_edge(i1, i2, G);
	}
      }
    }
    cout << "number of segment sides: " << num_vertices(G) << endl;
    cout << "number of segment-segment edges: " << num_edges(G) << endl;
    massert(edge_count > 0, "no pairs of segments were associated, consider reducing p threshold");

    int num = connected_components(G, &bins[0]);
    // cout << "connected components: " << num << endl;

    return num;
  }
}
