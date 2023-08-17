#ifndef BINMATRIX_H
#define BINMATRIX_H

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

#include <boost/math/distributions/chi_squared.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include "util.h"
#include "VariationSet.h"

using namespace std;
using namespace boost::math;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// BinSegmentSide
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct BinSegmentSide {
  string id;
  string contig;
  
  // coords
  int start, end;

  // each segment has two sides
  bool is_left_side;

  // length of segment used to compute counts 
  int side_length;
  
  // seed starting bin
  string seed_bin;

  // total number of supporting reads across all samples
  int supporting_reads;

  // total number of supporting samples with one or more supporting reads
  int supporting_samples;
  
  // library -> count
  vector < double > counts;
  
BinSegmentSide(string _id, string _contig, int _start, int _end,
	       bool _is_left_side, int _side_length) :
  id(_id), contig(_contig), start(_start), end(_end),
  is_left_side(_is_left_side), side_length(_side_length),
  seed_bin(""), supporting_reads(0), supporting_samples(0) {};
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// BinMatrix
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct BinMatrix {
private:
  vector< BinSegmentSide >& m_segs;
  int m_nsegs;

  chi_squared m_chi_squared;

  // distance matrix between all pairs of segments (p-values)
  vector < vector < bool > > m_selected;

  double m_chi_threshold;
  bool m_debug;
  
  double get_chi_value(vector < double >& c1, vector < double >& c2);
  double get_p_value(vector < double >& c1, vector < double >& c2);
  double compute_seg_chi_distance(BinSegmentSide& seg1, BinSegmentSide& seg2);

  // linear estimate of coverage ratio
  double compute_coverage_factor(BinSegmentSide& seg1, BinSegmentSide& seg2);
  
public:
  BinMatrix(vector< BinSegmentSide >& segs, int nlibs,
	    double p_threshold, bool debug);

  // init for specific range
  void init_matrix(int index, int from_ind, int to_ind);

  // init all matrix
  void init_matrix(int thread_count);

  // cluster matrix
  int cluster_segments(vector<int>& bins, vector< pair<int,int> >& seg_pairs);
};

#endif
