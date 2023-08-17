#include <iostream>
#include <random>
#include <numeric>
#include <cfloat>

#include "util.h"
#include "RefineLocal.h"

using namespace boost::math;

///////////////////////////////////////////////////////////////////////////////////////////////
// Breakpoint
///////////////////////////////////////////////////////////////////////////////////////////////

Breakpoint::Breakpoint(int _index, int _candidate_index, int _coord,
		       int _side_requested_length, double _p_value, int _read_support) :
  index(_index), candidate_index(_candidate_index), coord(_coord),
  side_requested_length(_side_requested_length), p_value(_p_value),
  read_support(_read_support)
{}

bool operator<(const Breakpoint& lhs, const Breakpoint& rhs)
{
  return lhs.coord < rhs.coord;
}

bool operator==(const Breakpoint& lhs, const Breakpoint& rhs)
{
  return lhs.coord == rhs.coord;
}

///////////////////////////////////////////////////////////////////////////////////////////////
// ctor and init
///////////////////////////////////////////////////////////////////////////////////////////////

RefineLocal::RefineLocal(string contig, vector < VariationSet >& cavs,
			 int margin, int min_side_length, int max_side_length,
			 double min_p, double pseudo_count, int candidate_step, bool debug)
  : m_contig(contig), m_cavs(cavs),
    m_margin(margin), m_min_side_length(min_side_length), m_max_side_length(max_side_length),
    m_min_p(min_p), m_pseudo_count(pseudo_count),
    m_chi_squared(cavs.size()-1), m_candidate_step(candidate_step),
    m_debug(debug)
{
  m_contig_length = m_cavs[0].get_covs()[contig].size();
  m_chi_threshold = quantile(m_chi_squared, 1-m_min_p);
  m_current_dist = m_min_side_length;  
}

void RefineLocal::init_dangle_support()
{
  for (unsigned int i=0; i<m_cavs.size(); ++i) {
    VariationSet& cav = m_cavs[i];
    map< string, map< int, map <Variation, int> > >& vars_all =  cav.get_vars();
    map< int, map <Variation, int> >& vars =  vars_all[m_contig];
    for (map< int, map <Variation, int> >::const_iterator it=vars.begin(); it!=vars.end(); ++it) {
      unsigned int coord = (*it).first;
      
      const map <Variation, int>& vars_coord = (*it).second;
      for (map <Variation, int>::const_iterator xt=vars_coord.begin(); xt!=vars_coord.end(); ++xt) {
	Variation var = (*xt).first;
	int count = (*xt).second;
	if (var.has_type(vtDangleLeft))
	  m_dangle_read_support[coord-1] += count;
	if (var.has_type(vtDangleRight))
	  m_dangle_read_support[coord] += count;
      }
    }
  }
}

void RefineLocal::init_candidate_coords(bool only_dangles)
{
  cout << "calling init_candidate_coords, N=" << m_contig_length << endl;

  // collect first in set
  set<int> coords;
  
  // add dangle coords
  for (map<int, int>::iterator it=m_dangle_read_support.begin(); it != m_dangle_read_support.end(); ++it) {
    int coord = (*it).first;
    if ((coord <= m_min_side_length+2*m_margin) || (coord >= m_contig_length-m_min_side_length-2*m_margin))
      continue;
    double chi_val = test_breakpoint(0, m_contig_length, coord);
    if (chi_val > m_chi_threshold)
      coords.insert(coord);
  }

  // sample all other coords
  if (!only_dangles) {
    int coord = m_candidate_step-1;
    while (coord < m_contig_length) {
      if ((coord <= m_min_side_length+2*m_margin) || (coord >= m_contig_length-m_min_side_length-2*m_margin)) {
	coord += m_candidate_step;
	continue;
      }
      double chi_val = test_breakpoint(0, m_contig_length, coord);
      if (chi_val > m_chi_threshold)
	coords.insert(coord);
      
      coord += m_candidate_step;
    }
  }

  copy(coords.begin(), coords.end(), back_inserter(m_candidates));
  cout << "number of candidates: " << m_candidates.size() << endl;
}

void RefineLocal::get_expected_probs(vector <double>& left_counts, vector <double>& right_counts,
				     double& left_p, double& right_p, vector <double>& lib_p, double& total)
{
  unsigned int nlibs = m_cavs.size();
  lib_p.resize(nlibs);
  
  double total_left=0, total_right=0;
  
  for (unsigned int i=0; i<nlibs; ++i) {
    total_left += left_counts[i];
    total_right += right_counts[i];
    lib_p[i] = left_counts[i] + right_counts[i];
  }
  total = (total_left + total_right);
  left_p = total_left / total;
  right_p = total_right / total;

  // normalize lib_p
  for (unsigned int i=0; i<nlibs; ++i)
    lib_p[i] = lib_p[i] / total;
}

void RefineLocal::get_breakpoint_side_coords(int start, int end, int coord, bool on_left_side,
					     int& seg_start, int& seg_end)
{
  int left_length = min(m_current_dist, coord-start-2*m_margin);
  int right_length = min(m_current_dist, end-coord-2*m_margin);
  massert(left_length >= m_min_side_length, "left-side segment length too short");
  massert(right_length >= m_min_side_length, "right-side segment length too short");
  
  seg_start = on_left_side ? coord-left_length-m_margin : coord+m_margin;
  seg_end   = on_left_side ? coord-m_margin             : coord+right_length+m_margin;
}

vector < double > RefineLocal::get_breakpoint_counts(int start, int end, int coord, bool on_left_side)
{
  int seg_start;
  int seg_end;
  get_breakpoint_side_coords(start, end, coord, on_left_side, seg_start, seg_end);
  
  vector < double > result;
  result.resize(m_cavs.size());
  for (unsigned int i=0; i<m_cavs.size(); ++i) {
    VariationSet& cav = m_cavs[i];
    result[i] = (double)cav.get_segment_coverage(m_contig, seg_start, seg_end) + m_pseudo_count;
  }
  return result;
}

double RefineLocal::test_breakpoint(int start, int end, int coord)
{
  // count tests for FDR
  // m_test_count++;
  stringstream ss;
  ss << coord << "_" << m_current_dist;
  string key = ss.str();
  m_test_count_map.insert(key);
  
  vector <double> left_counts = get_breakpoint_counts(start, end, coord, true);
  vector <double> right_counts = get_breakpoint_counts(start, end, coord, false);
  return chi_square_test(left_counts, right_counts, NULL);
}

int RefineLocal::inspect_segment(int start, int end,
				 int candidate_start_index,
				 int candidate_end_index,
				 int depth,
				 set<Breakpoint>& breakpoints)
{
  // cout << m_contig << " D=" << m_current_dist << " g=" << depth << " I=" << breakpoints.size()+1 << " start=" << start << " end=" << end << " " << "Si=" << candidate_start_index << " Ei=" << candidate_end_index << endl;
    
  // select candidate that matches criteria
  int max_coord = -1;
  double max_chi_val = 0;
  int max_index;
  int max_dangle_read_support = 0;
  
  for (int index=candidate_start_index; index<=candidate_end_index; ++index) {
    massert(index >= 0 && index < (int)m_candidates.size(),
	    "index %d out of range, n=%d", index, m_candidates.size());
    int coord = m_candidates[index];
    if ((coord <= start+m_min_side_length+2*m_margin) || (coord >= end-m_min_side_length-2*m_margin))
      continue;
    double chi_val = test_breakpoint(start, end, coord);
    if (chi_val <= m_chi_threshold)
      continue;

    // get read support
    int dangle_read_support = 0;
    if (m_dangle_read_support.find(coord) != m_dangle_read_support.end())
      dangle_read_support = m_dangle_read_support[coord];

    // select candidate with maximal suppport
    if (dangle_read_support > max_dangle_read_support) {
      max_dangle_read_support = dangle_read_support;
      max_coord = coord;
      max_chi_val = chi_val;
      max_index = index;
    }

    // for ties use chi-square test, for candidates without dangle support
    if (dangle_read_support == max_dangle_read_support && chi_val > max_chi_val) {
      max_coord = coord;
      max_chi_val = chi_val;
      max_index = index;
    }
  }

  // if none found return nothing
  if (max_coord == -1)
    return 0;

  // init breakpoint
  double p_val = 1 - cdf(m_chi_squared, max_chi_val);
  Breakpoint bp(breakpoints.size()+1, max_index, max_coord, m_current_dist, p_val, max_dangle_read_support);
  get_breakpoint_side_coords(start, end, max_coord, true, bp.left_side_start, bp.left_side_end);
  get_breakpoint_side_coords(start, end, max_coord, false, bp.right_side_start, bp.right_side_end);

  // cout << "adding: " << "c=" << max_coord << " index=" << max_index << endl;
  // add to container
  massert(breakpoints.find(bp) == breakpoints.end(), "breakpoint already in container, coord=%d", max_coord);
  breakpoints.insert(bp);

  int count = 1;
  
  // check left 
  if (max_index>candidate_start_index)
    count += inspect_segment(start, max_coord, candidate_start_index, max_index-1, depth+1, breakpoints);

  // check right 
  if (max_index<candidate_end_index)
    count += inspect_segment(max_coord, end, max_index+1, candidate_end_index, depth+1, breakpoints);

  return count;
}

void RefineLocal::get_breakpoints(bool only_dangles, set<Breakpoint>& breakpoints)
{
  init_dangle_support();
  init_candidate_coords(only_dangles);

  int last_cand_index = m_candidates.size()-1;
  
  while (m_current_dist <= m_max_side_length) {

    // no breakpoints found yet
    if (breakpoints.size() == 0) {
      inspect_segment(0, m_contig_length, 0, last_cand_index, 0, breakpoints);
      m_current_dist *= 2;
      continue;
    }

    // must travserse copy of breakpoints since we are modifying container as we go
    set<Breakpoint> prev_breakpoints(breakpoints);
    
    int prev_coord = 0;
    int prev_bp_index = -1;
    for (set<Breakpoint>::iterator it = prev_breakpoints.begin(); it != prev_breakpoints.end(); ++it) {
      int coord = (*it).coord;
      int bp_index = (*it).candidate_index;
      if (bp_index-prev_bp_index>=2)
	inspect_segment(prev_coord, coord, prev_bp_index+1, bp_index-1, 0, breakpoints);
      prev_coord = coord;
      prev_bp_index = bp_index;
    }
    // handle last fragment
    if (prev_bp_index<last_cand_index)
      inspect_segment(prev_coord, m_contig_length, prev_bp_index+1, last_cand_index, 0, breakpoints);

    m_current_dist *= 2;
  }
}

void RefineLocal::dump_scores(ofstream& out)
{
  cout << "init dangle support called" << endl;
  init_dangle_support();
  
  out << "dist\tcontig\tcoord\tdangle_support\tstat\tp_value" << endl;

  cout << "going over contig of length: " << m_contig_length << endl;
  m_current_dist = m_min_side_length;
  while (m_current_dist <= m_max_side_length) {

    cout << "current distance: " << m_current_dist << endl;
    for (int coord=0; coord<m_contig_length; ++coord) {
      if ((coord <= m_min_side_length+2*m_margin) || (coord >= m_contig_length-m_min_side_length-2*m_margin))
	continue;
      double chi_val = test_breakpoint(0, m_contig_length, coord);
      double p_val = 1 - cdf(m_chi_squared, chi_val);

      int dangle_read_support = 0;
      if (m_dangle_read_support.find(coord) != m_dangle_read_support.end())
	dangle_read_support = m_dangle_read_support[coord];
      
      out << m_current_dist << "\t" << m_contig << "\t" << coord+1 << "\t";
      out << dangle_read_support << "\t" << chi_val << "\t" << p_val << endl;
    }
    m_current_dist *= 2;
  }
}
