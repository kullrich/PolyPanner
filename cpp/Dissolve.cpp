#include <iostream>
#include <random>
#include <numeric>
#include <cfloat>

#include <boost/math/distributions/chi_squared.hpp>
using namespace boost::math;

#include "Dissolve.h"

///////////////////////////////////////////////////////////////////////////////////////////////
// ctor and init
///////////////////////////////////////////////////////////////////////////////////////////////

Dissolve::Dissolve(int lib_count, double pseudo_count, string weight_style)
  : m_init(false), m_lib_count(lib_count), m_pseudo_count(pseudo_count)
{
  if (weight_style == "uniform")
    m_weight_style = wsUniform;
  else if (weight_style == "marginal")
    m_weight_style = wsMarginal;
  else
    mexit("unknown weight_style: %s", weight_style.c_str());
}

vector < vector < double > >& Dissolve::get_counts() { return m_counts; }
vector < vector < double > >& Dissolve::get_freqs() { return m_freqs; }
int Dissolve::get_contig_length() { return m_counts.size(); }
void Dissolve::add_pseudo_count()
{
  for (int i=0; i<m_contig_length; ++i)
    for (int j=0; j<m_lib_count; ++j)
      m_counts[i][j] += m_pseudo_count;
}

void Dissolve::init_marg()
{
  m_marg.resize(m_contig_length);
  for (int i=0; i<m_contig_length; ++i)
    for (int j=0; j<m_lib_count; ++j)
      m_marg[i] += m_counts[i][j];
}

void Dissolve::init_weights()
{
  // init weights to 1
  m_weights.resize(m_contig_length);
  transform(m_weights.begin(), m_weights.end(), m_weights.begin(), [](const double t){ return 1; });
  m_weight_sum = accumulate(m_weights.begin(), m_weights.end(), 0.0);

  if (m_weight_style != wsMarginal)
    return;

  // marg mean
  double sum = accumulate(m_marg.begin(), m_marg.end(), 0.0);
  double mean = sum / m_marg.size();

  // marg sd
  double sq_sum = inner_product(m_marg.begin(), m_marg.end(), m_marg.begin(), 0.0);
  double sd = sqrt(sq_sum / m_marg.size() - mean * mean);

  // if sd is positive define weights using z-scores
  if (sd > DBL_MIN) {

    // absolute z-score
    transform(m_marg.begin(), m_marg.end(), m_weights.begin(), [mean,sd](const double t){ return abs(t-mean)/sd; });

    // sanity check
    for_each(m_weights.begin(), m_weights.end(), [](const double t){ if(!std::isfinite(t) || t<0) { cout << "assert failed"; exit(-1); } });

    // weight = 1/(|z|+1)
    transform(m_weights.begin(), m_weights.end(), m_weights.begin(), [](const double t){ return 1/(1+t); });

    // update weight sum
    m_weight_sum = accumulate(m_weights.begin(), m_weights.end(), 0.0);
  }

}

void Dissolve::init_freq()
{
  m_freqs.resize(m_contig_length);
  for (int i=0; i<m_contig_length; ++i) {
    m_freqs[i].resize(m_lib_count);

    // compute total for coord
    double total = 0;
    for (int j=0; j<m_lib_count; ++j)
      total += m_counts[i][j];

    // compute freq for coord
    for (int j=0; j<m_lib_count; ++j)
      m_freqs[i][j] = m_counts[i][j] / total;
  }
}

void Dissolve::init_center_coords()
{
  for (int i=0; i<m_contig_length; ++i)
    m_center_coords.insert(i);
}

void Dissolve::set_counts(string contig, vector < VariationSet >& cavs)
{
  for (unsigned int i=0; i<cavs.size(); ++i) {
    map<string, vector<int> >& set_cov = cavs[i].get_covs();
    massert(set_cov.find(contig) != set_cov.end(), "contig not found");
    vector<int>& contig_cov = set_cov[contig];
    unsigned int contig_length = contig_cov.size();

    if (i == 0) {
      m_counts.resize(contig_length);
      for (unsigned int j=0; j<contig_length; ++j)
	m_counts[j].resize(cavs.size());
    } else {
      massert(m_counts.size() == contig_length, "contig length does not match up between libraries");
    }
    for (unsigned int j=0; j<contig_length; ++j)
      m_counts[j][i] = contig_cov[j];
  }
}

void Dissolve::init_once()
{
  // init once
  if (m_init)
    return;
  m_init = true;

  massert(m_counts.size() > 0, "contig vector is empty");
  m_contig_length = m_counts.size();

  add_pseudo_count();

  init_freq();
  init_marg();
  init_weights();
  init_center_coords();
  compute_center_freq();
}

///////////////////////////////////////////////////////////////////////////////////////////////
// center functions
///////////////////////////////////////////////////////////////////////////////////////////////

void Dissolve::compute_center_freq()
{
  switch (m_weight_style) {
  case wsUniform:
    compute_center_freq_uniform(); break;
  case wsMarginal:
    compute_center_freq_marginal(); break;
  default:
    mexit("unknown weight style");
  }
}

void Dissolve::compute_center_freq_uniform()
{
  int center_size = m_center_coords.size();
  vector < double > result(m_lib_count);
  for (set <int >::iterator it=m_center_coords.begin(); it != m_center_coords.end(); ++it) {
    int coord = *it;
    massert(coord >= 0 && coord < m_contig_length, "coord out of range");
    for (int j=0; j<m_lib_count; ++j)
      result[j] += m_freqs[coord][j];
  }
  double total = 0;
  for (int j=0; j<m_lib_count; ++j) {
    result[j] = result[j] / center_size;
    total += result[j];
  }

  // extra normalize step due to possible precision issue
  for (int j=0; j<m_lib_count; ++j) {
    result[j] = result[j] / total;
  }

  m_center_freq = result;
}

vector<double> Dissolve::get_center_marginals()
{
  vector<double> result;
  for (set <int >::iterator it=m_center_coords.begin(); it != m_center_coords.end(); ++it) {
    int coord = *it;
    massert(coord >= 0 && coord < m_contig_length, "coord out of range");
    result.push_back(m_marg[coord]);
  }
  return result;
}

vector<double> Dissolve::get_center_values(int lib_index)
{
  vector<double> result;
  for (set <int >::iterator it=m_center_coords.begin(); it != m_center_coords.end(); ++it) {
    int coord = *it;
    massert(coord >= 0 && coord < m_contig_length, "coord out of range");
    result.push_back(m_freqs[coord][lib_index]);
  }
  return result;
}

vector<double> Dissolve::get_center_weights()
{
  vector<double> result;
  for (set <int >::iterator it=m_center_coords.begin(); it != m_center_coords.end(); ++it) {
    int coord = *it;
    massert(coord >= 0 && coord < m_contig_length, "coord out of range");
    result.push_back(m_weights[coord]);
  }
  return result;
}

void Dissolve::compute_center_freq_marginal()
{
  // weight among center positions
  vector<double> weights = get_center_weights();
  double w_sum = accumulate(weights.begin(), weights.end(), 0.0);
  transform(weights.begin(), weights.end(), weights.begin(), [w_sum](const double t){ return t/w_sum; });

  vector < double > result(m_lib_count);
  for (int i=0; i<m_lib_count; ++i) {
    vector<double> values = get_center_values(i);
    result[i] = inner_product(values.begin(), values.end(), weights.begin(), 0.0);
  }
  m_center_freq = result;
}

void Dissolve::remove_from_center(int coord)
{
  // if uniform we can save some time
  if (m_weight_style == wsUniform) {
    double center_size = m_center_coords.size();
    for (int j=0; j<m_lib_count; ++j)
      m_center_freq[j] = (center_size*m_center_freq[j] - m_freqs[coord][j]) / (center_size-1);
  } else {
    double coord_weight = m_weights[coord];
    double new_weight_sum = m_weight_sum - coord_weight;
    for (int j=0; j<m_lib_count; ++j)
      m_center_freq[j] = (m_weight_sum*m_center_freq[j] - coord_weight*m_freqs[coord][j]) / new_weight_sum;
    m_weight_sum = new_weight_sum;
  }
  m_center_coords.erase(coord);
}

///////////////////////////////////////////////////////////////////////////////////////////////
// outlier functions
///////////////////////////////////////////////////////////////////////////////////////////////

void Dissolve::sort_by_chi_values(vector< ChiValue >& chi_values)
{
  sort(chi_values.begin(), chi_values.end(), []( const ChiValue& lhs, const ChiValue& rhs )
       {
	 return lhs.value > rhs.value;
       });
}

void Dissolve::sort_by_coords(vector< ChiValue >& chi_values)
{
  sort(chi_values.begin(), chi_values.end(), []( const ChiValue& lhs, const ChiValue& rhs )
       {
	 return lhs.coord < rhs.coord;
       });
}

double Dissolve::get_chi_value(int coord)
{
  double stat = 0;
  for (int j=0; j<m_lib_count; ++j) {
    double obs = m_counts[coord][j];
    double exp = m_center_freq[j] * m_marg[coord];
    stat += (obs - exp) * (obs - exp) / exp;
  }
  return stat;
}


void Dissolve::get_center_chi_values(vector < ChiValue >& chi_values)
{
  chi_values.clear();
  for (set <int >::iterator it=m_center_coords.begin(); it != m_center_coords.end(); ++it) {
    int coord = *it;
    ChiValue cv(coord, get_chi_value(coord));
    chi_values.push_back(cv);
  }
}

void Dissolve::update_chi_values(vector < ChiValue >& chi_values)
{
  for (unsigned int i=0; i<chi_values.size(); ++i)
    chi_values[i].value = get_chi_value(chi_values[i].coord);
}

///////////////////////////////////////////////////////////////////////////////////////////////
// outlier functions
///////////////////////////////////////////////////////////////////////////////////////////////

void Dissolve::potential_outliers(double outlier_fraction, vector<ChiValue>& result)
{
  massert(m_init, "init must be called");

  result.clear();
  vector < ChiValue > chi_values;
  get_center_chi_values(chi_values);
  sort_by_chi_values(chi_values);
  int size = floor((double)chi_values.size() * outlier_fraction);
  if (size < 10)
    size = min(10, (int)chi_values.size());
  copy_n(chi_values.begin(), size, std::back_inserter(result));
}

int Dissolve::reduce_center_round(double outlier_fraction, double p_threshold)
{
  massert(m_init, "init must be called");

  chi_squared chi_squared(m_lib_count-1);
  double chi_threshold = quantile(chi_squared, 1-p_threshold);

  vector<ChiValue> chi_values;
  potential_outliers(outlier_fraction, chi_values);
  if (chi_values.size() == 0)
    return 0;

  int removed_count = 0;
  for (unsigned int i=0; i<chi_values.size(); ++i) {
    if (chi_values[i].value > chi_threshold) {
      remove_from_center(chi_values[i].coord);
      update_chi_values(chi_values);
      removed_count++;
    }
  }

  return removed_count;
}

void Dissolve::reduce_center(double outlier_fraction, double p_threshold)
{
  massert(m_init, "init must be called");

  int max_rounds = 1000;
  int round_count = 0;
  int removed_count = 1;
  while (removed_count > 0) {
    round_count++;
    removed_count = reduce_center_round(outlier_fraction, p_threshold);
    massert(round_count < max_rounds, "did not converge in %d rounds, giving up", max_rounds);
  }

  // remove short center segments
}

vector < double >& Dissolve::get_center_freq()
{
  massert(m_init, "init must be called");

  return m_center_freq;
}

set < int >& Dissolve::get_center_coords()
{
  massert(m_init, "init must be called");

  return m_center_coords;
}

vector < double > Dissolve::get_chivalues()
{
  massert(m_init, "init must be called");

  vector < double > result(m_contig_length);
  for (int i=0; i<m_contig_length; ++i)
    result[i] = get_chi_value(i);
  return result;
}

vector < double > Dissolve::get_pvalues()
{
  massert(m_init, "init must be called");

  vector < double > result(m_contig_length);
  vector < double > chivalues = get_chivalues();
  chi_squared chi_squared(m_lib_count-1);
  result.resize(m_contig_length);
  for (int i=0; i<m_contig_length; ++i)
    result[i] = 1 - cdf(chi_squared, chivalues[i]);
  return result;
}

vector < double > Dissolve::get_weights()
{
  return m_weights;
}

vector < double > Dissolve::get_marginals()
{
  return m_marg;
}

vector < bool > Dissolve::get_is_center()
{
  massert(m_init, "init must be called");
  vector < bool > result(m_contig_length);
  for (int i=0; i<m_contig_length; ++i)
    result[i]  = m_center_coords.find(i) != m_center_coords.end();
  return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// segment functions
///////////////////////////////////////////////////////////////////////////////////////////////

void Dissolve::get_segments_basic(vector<DissolveSegment>& segments)
{
  int prev_center_coord = -1;
  int center_start = 0;
  for (set<int>::iterator it=m_center_coords.begin(); it != m_center_coords.end(); ++it) {
    int center_coord = *it;

    // gap found
    if (center_coord > prev_center_coord + 1) {
      // center segment
      if (prev_center_coord != -1)
	segments.push_back(DissolveSegment(center_start, prev_center_coord+1, true));
      center_start = center_coord;

      // outlier segment
      if (center_coord > 0)
	segments.push_back(DissolveSegment(prev_center_coord+1, center_coord, false));
    }
    prev_center_coord = center_coord;
  }
  // close last center segment
  if (prev_center_coord != -1)
    segments.push_back(DissolveSegment(center_start, prev_center_coord+1, true));

  // add last outlier segment
  if (prev_center_coord+1 < m_contig_length)
    segments.push_back(DissolveSegment(prev_center_coord+1, m_contig_length, false));
}

void Dissolve::get_segments(vector<DissolveSegment>& segments, int min_center_seg_len)
{
  vector<DissolveSegment> basic_segments;
  get_segments_basic(basic_segments);
  if (basic_segments.size() == 1) {
    segments = basic_segments;
    return;
  }

  // correct short center segments
  vector<bool> remove(basic_segments.size());
  for (unsigned int i=0; i<basic_segments.size(); ++i) {
    remove[i] = basic_segments[i].is_center && (basic_segments[i].end - basic_segments[i].start) < min_center_seg_len;
    if (!remove[i])
      continue;

    // first/last/middle segment
    if (i == 0) {
      massert(!basic_segments[i+1].is_center, "expecting outlier segment after short center");
      basic_segments[i+1].start = basic_segments[i].start;
    } else if (i == basic_segments.size()-1) {
      massert(!basic_segments[i-1].is_center, "expecting outlier segment before short center");
      basic_segments[i-1].end = basic_segments[i].end;
    } else {
      massert(!basic_segments[i+1].is_center, "expecting outlier segment after short center");
      massert(!basic_segments[i-1].is_center, "expecting outlier segment before short center");
      // use next outlier segment and get rid of both segments
      remove[i-1] = true;
      basic_segments[i+1].start = basic_segments[i-1].start;
    }
  }

  // collect segments
  for (unsigned int i=0; i<basic_segments.size(); ++i)
    if (!remove[i])
      segments.push_back(basic_segments[i]);

  // extend outliers into adjacent center segments
  int ext = min_center_seg_len;
  for (unsigned int i=0; i<segments.size(); ++i) {
    if (segments[i].is_center)
      continue;
    if (i > 0) {
      massert(segments[i-1].is_center && segments[i-1].start <= segments[i].start-ext, "cannot extend outlier segment");
      segments[i].start -= ext;
    }
    if (i < segments.size()-1) {
      massert(segments[i+1].is_center && segments[i].end+ext <= segments[i+1].end, "cannot extend outlier segment");
      segments[i].end += ext;
    }
  }
}
