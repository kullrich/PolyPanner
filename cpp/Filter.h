#ifndef FILTER_H
#define FILTER_H

#include <boost/math/distributions/chi_squared.hpp>
using namespace boost::math;

#include <thread>
#include <mutex>

#include "VariationSet.h"

using namespace std;

enum SeqError { seNone, seSubstitute, seIndel, seRearrange, seCount };


struct FilterVariation {
  string contig;
  int coord;
  Variation var;
  int count;
  int total_count;
  double chi_value, p_value, q_value;
  FilterVariation(string _contig, int _coord, Variation _var,
		  int _count, int _total_count, double _chi_value, double _p_value) :
    contig(_contig), coord(_coord), var(_var),
    count(_count), total_count(_total_count),
    chi_value(_chi_value), p_value(_p_value), q_value(0) {};

  friend bool operator==(const FilterVariation& lhs, const FilterVariation& rhs);
  friend bool operator<(const FilterVariation& lhs, const FilterVariation& rhs);
  friend ostream& operator<<(ostream& os, const FilterVariation& fv);
};

class Filter {
 private:
  // cav object
  const VariationSet m_cav;
  
  // min cov of considered variant
  int m_min_cov;

  // min contig length
  int m_min_contig_length;
  
  // min/max/seed error rate
  double m_min_error, m_max_error, m_seed_error;

  // inferred error rates per type
  vector<double> m_rates;

  // inferred segregating sites
  vector < FilterVariation > m_sites;

  // requested false-discovery rate
  double m_alpha;

  // requested maximal p_value
  double m_p_threshold;
  
  // LL test p-values are distributed as a chi-square with df=1
  chi_squared m_chi_squared;
  // double m_chi_threshold;

  int m_round;

  // prefix of progress file
  string m_oprefix;

  // number of threads used
  int m_thread_count;
  
  // we keep contig vectors for threading activity
  mutex m_mtx; 
  vector<string> m_site_contigs;
  vector<string> m_rate_contigs;

  map < string, set < int > > m_masked;

  // used while estimating error rates
  vector<double> m_error_counts;
  double m_total_count;

  //////////////////////////////////////////////////////////////////////////////
  // rate inference
  //////////////////////////////////////////////////////////////////////////////

  SeqError get_error_type(Variation var);
  double get_error_rate(Variation src, Variation tgt);

  double get_weighted_error_rate(int true_index, Variation tgt,
				 vector < pair < Variation, int > >& var_counts);

  void count_total_count();  
  void count_var_errors();

  // infer error rates, keeping sites fixed
  void infer_rates();
  
  //////////////////////////////////////////////////////////////////////////////
  // site inference
  //////////////////////////////////////////////////////////////////////////////

  // compute hypothesis p-values
  double get_null_p_value(vector < pair < Variation, int > >& var_counts);
  double get_p_value(int H_index, vector < pair < Variation, int > >& var_counts);
  
  // compute p-values for all variants in position
  void test_sites_on_coord(string contig, int coord);

  double ll_test(int H_index, vector < pair < Variation, int > >& var_counts);
  
  // apply BH correction
  void apply_BH();

  // infer sites, keeping error rates fixed
  // returns true of sites remained the same
  bool infer_sites();

 public:
  Filter(VariationSet& cav, int min_cov, int min_contig_length,
	 double min_error, double max_error, double seed_error,
	 double alpha, int thread_count);

  // infer sites and error rates
  void resolve(int max_iterations, string oprefix);

  // thread task, made public to allow thread to call it
  void infer_sites(int index, int from_ind, int to_ind);
  void count_total_count(int index, int from_ind, int to_ind);
  void count_var_errors(int index, int from_ind, int to_ind);

  void save_error_rates(string ofn);
  void save_sites(string ofn);

};

#endif
