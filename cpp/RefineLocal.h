#ifndef REFINE_LOCAL_H
#define REFINE_LOCAL_H

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

#include <boost/math/distributions/chi_squared.hpp>

#include "util.h"
#include "VariationSet.h"

using namespace std;
using namespace boost::math;

struct Breakpoint {

  // order of introduction 
  int index;

  // index in candidate vector
  int candidate_index;
  
  // sort by coord
  int coord;

  // side coordinates
  int side_requested_length;

  // p-value of chi-sqaure between side distribs
  double p_value;

  // reads that supported dangle
  int read_support;
  
  // side coordinates
  int left_side_start, left_side_end, right_side_start, right_side_end;

  friend bool operator<(const Breakpoint& lhs, const Breakpoint& rhs);
  friend bool operator==(const Breakpoint& lhs, const Breakpoint& rhs);
  
  Breakpoint(int _index, int _candidate_index, int _coord, int _side_requested_length,
	     double _p_value, int _read_support);
};
  
class RefineLocal {
 private:
  string m_contig;
  
  // basic matrix data
  vector < VariationSet >& m_cavs;

  // margin away from candidate breakpoint and segment sides
  int m_margin;
  
  // length of inspected segment near breakpoint
  int m_current_dist;
  int m_min_side_length;
  int m_max_side_length;

  // classify as breakpoint of under this p-value
  double m_min_p;
  
  // pseudo-count to avoid infinite probs
  double m_pseudo_count;

  // contig length
  int m_contig_length;

  chi_squared m_chi_squared;

  double m_chi_threshold;

  // candidate coords
  vector<int> m_candidates;
  
  // from coord to read support of dangle
  map<int, int> m_dangle_read_support;

  // avoid counting same test twice
  set <string> m_test_count_map;

  int m_candidate_step;
  
  bool m_debug;

  // get support of coords with dangles
  void init_dangle_support();
  
  // all coords along contig
  void init_candidate_coords(bool only_dangles);
  
  void get_breakpoint_side_coords(int start, int end, int coord, bool on_left_side,
				  int& seg_start, int& seg_end);
  
  // get counts left or right or candidate breakpoint
  vector < double > get_breakpoint_counts(int start, int end, int coord, bool on_left_side);


  void get_expected_probs(vector <double>& left_counts, vector <double>& right_counts,
			  double& left_p, double& right_p, vector <double>& lib_p, double& total);

  // recursive function, inspecting candidate breakpoints with a segment 
  int inspect_segment(int start, int end,
		      int candidate_start_index,
		      int candidate_end_index,
		      int depth,
		      set<Breakpoint>& breakpoints);
    
  // test breakpoint candidate
  double test_breakpoint(int start, int end, int coord);
  
 public:
  RefineLocal(string contig, vector < VariationSet >& cavs,
	      int margin, int min_side_length, int max_side_length,
	      double min_p, double pseudo_count, int candidate_step, bool debug);

  void get_breakpoints(bool only_dangles, set<Breakpoint>& breakpoints);

  int get_test_count() { return m_test_count_map.size(); };
  // dump stats and pvalues
  void dump_scores(ofstream& out);
};

#endif
