#include <iostream>
#include <random>
#include <numeric>
#include <cfloat>

#include <gsl/gsl_randist.h>

#include "util.h"
#include "Filter.h"

///////////////////////////////////////////////////////////////////////////////////////////////
// FilterVariation
///////////////////////////////////////////////////////////////////////////////////////////////

bool operator<(const FilterVariation& lhs, const FilterVariation& rhs)
{
  return lhs.p_value < rhs.p_value;
}

bool operator==(const FilterVariation& lhs, const FilterVariation& rhs)
{
  return (lhs.contig == rhs.contig && lhs.coord == rhs.coord && lhs.var == rhs.var);
}

ostream& operator<<(ostream& os, const FilterVariation& fv)
{
  os << fv.contig << ":" << fv.coord << ":" << fv.var;
  os << ",n=" << fv.count << ",N=" << fv.total_count;
  os << ",cs=" << fv.chi_value << ",p=" << fv.p_value << ",q=" << fv.q_value;
  return os;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// ctor and init
///////////////////////////////////////////////////////////////////////////////////////////////

Filter::Filter(VariationSet& cav, int min_cov, int min_contig_length,
	       double min_error, double max_error, double seed_error,
	       double alpha, int thread_count)
  : m_cav(cav), m_min_cov(min_cov), m_min_contig_length(min_contig_length),
    m_min_error(min_error), m_max_error(max_error), m_seed_error(seed_error),
    m_alpha(alpha), m_chi_squared(1), m_round(1), m_thread_count(thread_count)
{
  m_rates.resize(seCount);
  m_error_counts.resize(seCount);
  
  double total = 0;
  for (unsigned int i=1; i<m_rates.size(); ++i) {
    m_rates[i] = m_seed_error;
    total += m_seed_error;
  }
  m_rates[seNone] = 1 - total;

  // m_chi_threshold = quantile(chi_squared, 1-p_threshold);
}

///////////////////////////////////////////////////////////////////////////////////////////////
// error rates
///////////////////////////////////////////////////////////////////////////////////////////////

string seq_error_to_string(SeqError se)
{
  switch(se) {
  case seNone: return "seNone";
  case seSubstitute: return "seSubstitute";
  case seIndel: return "seIndel";
  case seRearrange: return "seRearrange";
  default: mexit("invalid SeqError, index=%d", se);
  }
  return "error";
}

SeqError Filter::get_error_type(Variation var)
{
  if (var.has_type(vtDangleLeft) || var.has_type(vtDangleRight))
      return seRearrange;

  if (var.has_type(vtDelete) || var.has_type(vtInsert))
    return seIndel;

  if (var.has_type(vtSubstitute))
    return seSubstitute;

  if (var.is_ref())
    return seNone;
  
  massert(0, "unhandled variant: %s", var.to_string().c_str());
  return seCount;
}

double Filter::get_error_rate(Variation src, Variation tgt)
{
  if (src == tgt)
    return m_rates[seNone];
  
  if (src.has_type(vtDangleLeft) || src.has_type(vtDangleRight) ||
      tgt.has_type(vtDangleLeft) || tgt.has_type(vtDangleRight))
    return m_rates[seRearrange];

  if (src.has_type(vtDelete) || src.has_type(vtInsert) ||
      tgt.has_type(vtDelete) || tgt.has_type(vtInsert))
    return m_rates[seIndel];

  if (src.has_type(vtSubstitute) || tgt.has_type(vtSubstitute))
    return m_rates[seSubstitute];
  
  massert(0, "unhandled variant pair: %s -> %s", src.to_string().c_str(), tgt.to_string().c_str());
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// site inference
///////////////////////////////////////////////////////////////////////////////////////////////

double Filter::get_null_p_value(vector < pair < Variation, int > >& var_counts)
{
  unsigned int K = var_counts.size();
  vector<double> p(K+1);
  vector<unsigned int> n(K+1);
  double tot_p = 0;
  for (unsigned int i=0; i<K; ++i) {
    n[i] = var_counts[i].second;
    p[i] = get_error_rate(var_counts[0].first, var_counts[i].first);
    tot_p += p[i];
  }
  n[K] = 0;
  p[K] = 1-tot_p;
  double result = gsl_ran_multinomial_lnpdf(K+1, &p[0], &n[0]);

  return result;
}

double Filter::get_p_value(int H_index, vector < pair < Variation, int > >& var_counts)
{
  // number of variants
  unsigned int K = var_counts.size();

  // total count
  //  note: this code should be moved below into the loop so that weights sum to 1
  int N = 0;
  for (unsigned int i=0; i<K; ++i)
    N += var_counts[i].second;
  
  vector<double> p(K+1);
  vector<unsigned int> n(K+1);
  double tot_p = 0;
  for (unsigned int i=0; i<K; ++i) {
    n[i] = var_counts[i].second;
    p[i] = 0;
    for (int j=0; j<=H_index; ++j) {
      double weight = (double)var_counts[j].second / N;
      double pj = get_error_rate(var_counts[j].first, var_counts[i].first);
      p[i] += weight * pj;
    }
    tot_p += p[i];
  }
  n[K] = 0;
  p[K] = 1-tot_p;
  double result = gsl_ran_multinomial_lnpdf(K+1, &p[0], &n[0]);

  return result;
}

// variants up to index are true
double Filter::get_weighted_error_rate(int true_index, Variation tgt,
				       vector < pair < Variation, int > >& var_counts)
{
  if (true_index == 0)
      return get_error_rate(var_counts[0].first, tgt);

  double N = 0;
  for (unsigned int i=0; i<var_counts.size(); ++i)
    N += var_counts[i].second;
  
  double result = 0;
  for (int i=0; i<=true_index; ++i) {
    double weight = (double)var_counts[i].second / N;
    double p = get_error_rate(var_counts[i].first, tgt);
    result += weight * p;
  }
  return result;
}

double Filter::ll_test(int H_index, vector < pair < Variation, int > >& var_counts)
{
  massert(H_index > 0, "H_index must be greater than 0");
  // cout << "get_ll_test: " << endl;
  
  // number of variants
  unsigned int K = var_counts.size();

  double result = 0;
  for (unsigned int i=0; i<K; ++i) {
    double n_i = var_counts[i].second;
    double p_null_i = get_weighted_error_rate(H_index-1, var_counts[i].first, var_counts);
    double p_i = get_weighted_error_rate(H_index, var_counts[i].first, var_counts);
    // cout << "ni=" << n_i << ", p_null_i=" << p_null_i << ", p_i=" << p_i << endl;
    result += n_i*(log(p_null_i) - log(p_i));
  }
  return -2 * result;
}

void Filter::test_sites_on_coord(string contig, int coord)
{  
  vector < pair < Variation, int > > var_counts;
  m_cav.get_var_counts(contig, coord, false, var_counts);

  // massert(var_counts.size() > 1, "expecting more than one variant in position: %s: %d", contig.c_str(), coord);

  if (var_counts.size() < 2 || var_counts[1].second < m_min_cov)
    return;

  //  if (var_counts.size() != 1)
  // cout << "contig=" << contig << ", coord=" << coord << ", N=" << var_counts.size() << endl;
  
  // double log_prev_P = get_null_p_value(var_counts);

  int total_count = 0;
  for (unsigned int i=0; i<var_counts.size(); ++i)
    total_count += var_counts[i].second;
  
  for (unsigned int i=1; i<var_counts.size(); ++i) {
    Variation var = var_counts[i].first;
    int count = var_counts[i].second;

    // stop if under minimal number of reads
    if (count < m_min_cov)
      break;

    double LL = ll_test(i, var_counts);
    
    if (LL < 0) LL = 0;
    double p_val = 1 - cdf(m_chi_squared, LL);
    
    m_mtx.lock();
    m_sites.push_back(FilterVariation(contig, coord, var, count, total_count, LL, p_val));
    m_mtx.unlock();
  }
}

void Filter::apply_BH()
{
  //  for(unsigned int i=0; i<m_sites.size(); ++i)
  //    cout << "i:" << i << " " << m_sites[i] << endl;
  
  // find largest k that satifies condition
  double N = m_sites.size();
  double K = N-1;
  while(m_sites[K].p_value > ((m_alpha*K) / N))
    K--;

  // cout << "Number of p-values: " << N << endl;
  // cout << "Max K: " << K << endl;

  cout << "Maximal P-value: " << m_sites[K].p_value << endl;
  
  // limit sites
  vector < FilterVariation > result;
  for (int i=0; i<=K; ++i) {
    m_sites[i].q_value = (m_sites[i].p_value*N)/(i+1);
    result.push_back(m_sites[i]);
  }
  m_sites = result;
}

void Filter::infer_sites(int index, int from_ind, int to_ind)
{
  const map< string, map< int, map <Variation, int> > >& vars = m_cav.get_vars();
  for (int i=from_ind; i<to_ind; ++i) {
    string contig = m_site_contigs[i];
  
    // process putative sites
    const map< int, map <Variation, int> >& vars_contig = vars.at(contig);
    for (map< int, map <Variation, int> >::const_iterator jt=vars_contig.begin(); jt != vars_contig.end(); ++jt) {
      int coord = (*jt).first;
      test_sites_on_coord(contig, coord);
    }
  }
}

void global_filter_infer_sites(Filter* filter, int index, int from_ind, int to_ind)
{
  filter->infer_sites(index, from_ind, to_ind);
}

bool Filter::infer_sites()
{
  vector < FilterVariation > prev_sites = m_sites;
  m_sites.clear();
  m_site_contigs.clear();

  // prepare for threading
  const map< string, map< int, map <Variation, int> > >& vars = m_cav.get_vars();
  for (map< string, map< int, map <Variation, int> > >::const_iterator it=vars.begin(); it!=vars.end(); ++it) {
    string contig = (*it).first;
    if (m_cav.get_contig_length(contig) < m_min_contig_length)
      continue;
    m_site_contigs.push_back(contig);
  }

  int n_contigs = m_site_contigs.size();
  int thread_count = min(m_thread_count, n_contigs);
  
  int step = floor(n_contigs / thread_count);
  // cout << "number of contigs: " << n_contigs << endl;
  // cout << "number of threads: " << thread_count << endl;

  vector<thread> threads;
  for (int i=0; i<thread_count; ++i) {
    int from_ind = i * step;
    int to_ind = (i < (thread_count-1)) ? (i+1) * step : n_contigs;
    threads.push_back(thread(global_filter_infer_sites, this, i, from_ind, to_ind));
  }

  cout << "waiting for " << thread_count << " site inference threads to finish\n";
  for (auto& th : threads) th.join();

  sort(m_sites.begin(), m_sites.end());
  
  // finally save result
  save_sites(m_oprefix + "." + to_string(m_round));
  cout << "number of sites before BH: " << m_sites.size() << endl;
  apply_BH();
  cout << "number of sites found after BH: " << m_sites.size() << endl;

  // check if sites remained the same
  set<FilterVariation> s1(m_sites.begin(), m_sites.end()), s2(prev_sites.begin(), prev_sites.end());

  // !!!
  // return (s1 == s2);
  return (m_sites.size() == prev_sites.size());
}

///////////////////////////////////////////////////////////////////////////////////////////////
// rate inference
///////////////////////////////////////////////////////////////////////////////////////////////

void Filter::count_var_errors(int index, int from_ind, int to_ind)
{
  const map< string, map< int, map <Variation, int> > >& vars = m_cav.get_vars();
  for (int i=from_ind; i<to_ind; ++i) {
    string contig = m_rate_contigs[i];
    
    const map< int, map <Variation, int> >& vars_contig = vars.at(contig);
    set < int >& masked_contig = m_masked[contig];
    for (map< int, map <Variation, int> >::const_iterator jt=vars_contig.begin(); jt != vars_contig.end(); ++jt) {
      int coord = (*jt).first;
      
      // skip if masked due to segregating sites
      if (masked_contig.find(coord) != masked_contig.end())
	continue;

      // get all counts
      vector < pair < Variation, int > > vars;
      m_cav.get_var_counts(contig, coord, false, vars);
      massert(vars.size() > 0, "expecting at least one element");

      // skip if reference isn't major allele
      if (!vars[0].first.is_ref()) {
	m_mtx.lock();
	m_masked[contig].insert(coord);
	m_mtx.unlock();
	continue;
      }

      for (unsigned int i=1; i<vars.size(); ++i) {
	SeqError se = get_error_type(vars[i].first);
	int count = vars[i].second;
	m_mtx.lock();
	m_error_counts[se] += count;
	m_mtx.unlock();
      }
    }
  }
}

void global_filter_count_var_errors(Filter* filter, int index, int from_ind, int to_ind)
{
  filter->count_var_errors(index, from_ind, to_ind);
}

void Filter::count_var_errors()
{
  for (unsigned int i=0; i<seCount; ++i)
    m_error_counts[i] = 0;
  m_rate_contigs.clear();
  
  const map< string, map< int, map <Variation, int> > >& vars = m_cav.get_vars();
  for (map< string, map< int, map <Variation, int> > >::const_iterator it=vars.begin(); it!=vars.end(); ++it) {
    string contig = (*it).first;
    if (m_cav.get_contig_length(contig) < m_min_contig_length)
      continue;
    m_rate_contigs.push_back(contig);
  }

  int n_contigs = m_rate_contigs.size();
  int thread_count = min(m_thread_count, n_contigs);
  
  int step = floor(n_contigs / thread_count);
  // cout << "number of contigs: " << n_contigs << endl;
  // cout << "number of threads: " << thread_count << endl;

  vector<thread> threads;
  for (int i=0; i<thread_count; ++i) {
    int from_ind = i * step;
    int to_ind = (i < (thread_count-1)) ? (i+1) * step : n_contigs;
    threads.push_back(thread(global_filter_count_var_errors, this, i, from_ind, to_ind));
  }

  cout << "waiting for " << thread_count << " variant errors threads to finish\n";
  for (auto& th : threads) th.join();
}

void Filter::count_total_count(int index, int from_ind, int to_ind)
{
  const map<string, vector<int> >& covs = m_cav.get_covs();
  for (int i=from_ind; i<to_ind; ++i) {
    string contig = m_rate_contigs[i];
    const vector<int>& cov_contig = covs.at(contig);
    set < int >& masked_contig = m_masked[contig];
    for (unsigned int coord=0; coord < cov_contig.size(); ++coord) {
      if (masked_contig.find(coord) != masked_contig.end())
	continue;
      m_mtx.lock();
      m_total_count += cov_contig[coord];
      m_mtx.unlock();
    }      
  }
}

void global_filter_count_total(Filter* filter, int index, int from_ind, int to_ind)
{
  filter->count_total_count(index, from_ind, to_ind);
}

void Filter::count_total_count()
{
  m_rate_contigs.clear();
  
  // go over all positions to count total coverage
  m_total_count = 0;
  const map<string, vector<int> >& covs = m_cav.get_covs();
  for (map<string, vector<int> >::const_iterator it=covs.begin(); it!=covs.end(); ++it) {
    string contig = (*it).first;
    if (m_cav.get_contig_length(contig) < m_min_contig_length)
      continue;
    m_rate_contigs.push_back(contig);
  }

  int n_contigs = m_rate_contigs.size();
  int thread_count = min(m_thread_count, n_contigs);
  
  int step = floor(n_contigs / thread_count);
  // cout << "number of contigs: " << n_contigs << endl;
  // cout << "number of threads: " << thread_count << endl;

  vector<thread> threads;
  for (int i=0; i<thread_count; ++i) {
    int from_ind = i * step;
    int to_ind = (i < (thread_count-1)) ? (i+1) * step : n_contigs;
    threads.push_back(thread(global_filter_count_total, this, i, from_ind, to_ind));
  }

  cout << "waiting for " << thread_count << " total coverage threads to finish\n";
  for (auto& th : threads) th.join();
}

void Filter::infer_rates()
{
  // init mask sites
  m_masked.clear();
  for (unsigned int i=0; i<m_sites.size(); ++i)
    m_masked[m_sites[i].contig].insert(m_sites[i].coord);

  count_var_errors();
  count_total_count();
  
  double total_p = 0;
  for (unsigned int i=1; i<seCount; ++i) {
    SeqError se = (SeqError)i;
    cout << seq_error_to_string(se) << " n=" << m_error_counts[i] << endl;
    double new_rate = m_error_counts[i] / m_total_count;
    if (new_rate < m_min_error)
      new_rate = m_min_error;
    if (new_rate > m_max_error && i != 0)
      new_rate = m_max_error;
    cout << "updating rate of " << seq_error_to_string(se) << ": " <<
      m_rates[i] << " -> " << new_rate << endl;
    m_rates[i] = new_rate;
    total_p += new_rate;
  }
    
  cout << "updating rate of " << seq_error_to_string(seNone) << ": " <<
    m_rates[seNone] << " -> " << 1 - total_p << endl;
  m_rates[seNone] = 1 - total_p;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// Save output
///////////////////////////////////////////////////////////////////////////////////////////////

void Filter::save_error_rates(string ofn)
{
  cout << "saving error rate file: " << ofn << endl;
  ofstream out(ofn.c_str(), ios::out);
  massert(out.is_open(), "could not open file %s", ofn.c_str());

  out << "type\trate" << endl;
  for (unsigned int i=0; i<seCount; ++i) {
    SeqError se = (SeqError)i;
    out << seq_error_to_string(se) << "\t" << m_rates[i] << endl;
  }
    
  out.close();
}

void Filter::save_sites(string ofn)
{
  cout << "saving site file: " << ofn << endl;
  ofstream out(ofn.c_str(), ios::out);
  massert(out.is_open(), "could not open file %s", ofn.c_str());

  out << "contig\tcoord\tvariant\tcount\ttotal_count\tchival\tmlog_pval\tmlog_qval" << endl;
  for (unsigned int i=0; i<m_sites.size(); ++i) {
    FilterVariation& fv = m_sites[i];
    out << fv.contig << "\t" << fv.coord+1 << "\t" << fv.var << "\t";
    out << fv.count << "\t" << fv.total_count << "\t";
    out << fv.chi_value << "\t" << -log10(fv.p_value) << "\t" << -log10(fv.q_value) << endl;
  }
    
  out.close();
}

///////////////////////////////////////////////////////////////////////////////////////////////
// main function
///////////////////////////////////////////////////////////////////////////////////////////////

void Filter::resolve(int max_iterations, string oprefix)
{
  m_oprefix = oprefix;
  
  bool converged = false;
  for (int i=0; i<max_iterations; ++i) {
    m_round = i+1;
    cout << "================================================================================" << endl;
    cout << "ITERATION=" << i+1 << endl;
    cout << "infering sites..." << endl;
    converged = infer_sites();
    cout << "inferring rates..." << endl;
    infer_rates();
    if (converged)
      break;
  }
  massert(converged, "did not converge after %d iterations", max_iterations);
}  
