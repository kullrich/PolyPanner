#ifndef DISSOLVE_H
#define DISSOLVE_H

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

#include "util.h"
#include "VariationSet.h"

using namespace std;

enum WeightStyle { wsUniform, wsMarginal };

struct DissolveSegment {
  int start, end;
  bool is_center;
  DissolveSegment(int _start, int _end, bool _is_center) : start(_start), end(_end), is_center(_is_center) {};
};

struct ChiValue {
  int coord;
  double value;
  ChiValue(int _coord, double _value) : coord(_coord), value(_value) {};
};

class Dissolve {
 private:
  // basic matrix data
  vector < vector < double > > m_counts;
  vector < vector < double > > m_freqs;

  // marginal count per coordinate
  vector < double > m_marg;

  // 1 /(1+|z|) of marginals
  vector < double > m_weights;
  double m_weight_sum;

  // init once
  bool m_init;

  // dimensions
  int m_contig_length;
  int m_lib_count;

  // pseudo-count to avoid infinite probs
  double m_pseudo_count;
  
  WeightStyle m_weight_style;
  
  // center coords (i.e. not outliers)
  set <int > m_center_coords;

  // center freq
  vector <double > m_center_freq;

  void add_pseudo_count();

  // init freq from counts
  void init_freq();

  // init marg per coord
  void init_marg();

  void init_weights();

  ////////////////////////////////////////////////
  // center functions
  ////////////////////////////////////////////////

  // init center using all
  void init_center_coords();

  // update center using selected coords
  void compute_center_freq_uniform();
  void compute_center_freq();

  // update center by removing specific coord
  void remove_from_center(int coord);

  // compute weighted center frequency vector
  vector<double> get_center_values(int lib_index);
  vector<double> get_center_marginals();
  vector<double> get_center_weights();
  void compute_center_freq_marginal();

  ////////////////////////////////////////////////
  // chi_value functions
  ////////////////////////////////////////////////

  // sorting functions
  void sort_by_chi_values(vector < ChiValue >& chi_values);
  void sort_by_coords(vector < ChiValue >& chi_values);

  // compute single Chi-Square P-value
  double get_chi_value(int coord);

  // update chi_values for specific coords
  void update_chi_values(vector < ChiValue >& chi_values);

  ////////////////////////////////////////////////
  // outlier functions
  ////////////////////////////////////////////////

  void potential_outliers(double outlier_fraction, vector<ChiValue>& result);

  // compute chi_values for all coords
  void get_center_chi_values(vector < ChiValue >& chi_values);

  // extract segments without removing short center segments
  void get_segments_basic(vector<DissolveSegment>& segments);

 public:
  Dissolve(int lib_count, double pseudo_count, string weight_style);

  // matrix structure: coord -> library
  vector < vector < double > >& get_counts();
  vector < vector < double > >& get_freqs();

  void set_counts(string contig, vector < VariationSet >& cavs);

  int get_contig_length();

  // call after counts were directly set
  void init_once();

  // get center frequency
  vector < double >& get_center_freq();
  set < int >& get_center_coords();

  int reduce_center_round(double outlier_fraction, double p_threshold);
  void reduce_center(double outlier_fraction, double p_threshold);

  vector < double > get_chivalues();
  vector < double > get_pvalues();
  vector < bool > get_is_center();
  vector < double > get_weights();
  vector < double > get_marginals();

  void get_segments(vector<DissolveSegment>& segments, int min_center_segment_length);
};

#endif
