#ifndef CLUSTER_VARIANT_H
#define CLUSTER_VARIANT_H

#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <fstream>
#include <assert.h>
#include <sstream>
#include <stdarg.h>
#include <math.h>

#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

#include <dirent.h>

#include <thread>
#include <mutex>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include "util.h"
#include "VariationSet.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ClusterVariantItem
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct ClusterVariantItem {
  string id;
  string contig;
  int coord;
  Variation var;

  int cluster_index;
  
  // library -> count
  vector < double > counts;
  
ClusterVariantItem(string _id, string _contig, int _coord, Variation _var) :
  id(_id), contig(_contig), coord(_coord), var(_var), cluster_index(-1) {};
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ClusterVariant
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct ClusterVariants {
private:
  string m_title;
  double m_p_threshold;
  vector< ClusterVariantItem >& m_vars;
  int m_nvars;

  std::mutex m_mtx;
  
  // distance matrix between all pairs of variants
  vector < vector < bool > > m_selected;

  // init all with multiple threads
  void init_matrix(int thread_count);
public:
  ClusterVariants(string title,
		  double p_threshold,
		  vector< ClusterVariantItem >& vars,
		  int nlibs,
		  int thread_count);

  // results are stored in the cluster_index field of m_vars
  int cluster_variants(vector<int>& vclusters);

  // init only slice of matrix with single threads
  void init_matrix(int index, int from_ind, int to_ind);
};

#endif
