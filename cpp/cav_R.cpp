#include <Rcpp.h>
#include "VariationSet.h"

using namespace Rcpp;
using namespace std;

// echo PKG_CPPFLAGS=-Wall -Wno-write-strings -std=c++0x > ~/.R/Makevars

/*
library("Rcpp")
Sys.setenv(PKG_CPPFLAGS="-Wall -Wno-write-strings -std=c++0x")
sourceCpp("cav_R.cpp", verbose=T)
cav = cav_load("/relman03/work/users/eitany/tempo/subjects/AAB/assembly/megahit/cav_v4/sets_N1/sets/s1/set.cav")
df = data.frame(contig="s1", from=1, to=10)
cav_query(cav, df)
cc = cav_contigs(cav)
*/

// [[Rcpp::export]]
Rcpp::XPtr<VariationSet> cav_load(string fn)
{
  VariationSet* varset = new VariationSet;
  varset->load(fn);
  Rcpp::XPtr<VariationSet> result(varset, true);
  return result;
}

// [[Rcpp::export]]
vector<string> cav_contigs(Rcpp::XPtr<VariationSet> xptr)
{
  Rcpp::XPtr<VariationSet> varset(xptr);
  return varset->get_contigs();
}

// [[Rcpp::export]]
DataFrame cav_query(Rcpp::XPtr<VariationSet> xptr, DataFrame df)
{
  Rcpp::XPtr<VariationSet> varset(xptr);
  StringVector contig_v = df["contig"];
  IntegerVector from_v = df["from"];
  IntegerVector to_v = df["to"];

  vector<string> contig_r;
  vector<int> coord_r;
  vector<string> var_r;
  vector<int> count_r;
  vector<int> total_r;
  vector<int> cc_start_r;
  vector<int> cc_end_r;

  for (unsigned int i=0; i<contig_v.size(); ++i) {
    string contig = as<string>(contig_v[i]);
    int from = from_v[i] - 1;
    int to = to_v[i] - 1;
    varset->get_contig_data(contig, from, to, contig_r, coord_r, var_r, count_r, total_r, cc_start_r, cc_end_r);
  }

  // switch back to 1-based
  for (unsigned int i=0; i<coord_r.size(); ++i)
    coord_r[i]++;

  DataFrame result = DataFrame::create(
				       Named("contig") = contig_r,
				       Named("coord") = coord_r,
				       Named("var") = var_r,
				       Named("count") = count_r,
				       Named("total") = total_r,
				       Named("cc_start") = cc_start_r,
				       Named("cc_end") = cc_end_r);

  return result;
}

/*
library("Rcpp")
Sys.setenv(PKG_CPPFLAGS="-Wall -Wno-write-strings -std=c++0x")
sourceCpp("cav_R.cpp", verbose=T)
cav = cav_load("/relman03/work/users/eitany/tempo/subjects/AAB/assembly/megahit/cav_v4/sets_N1/sets/s1/hosts.cav")
df = data.frame(contig="s1", from=1, to=1000)
cav_summary(cav, df, 2, 0.8)
 */

// [[Rcpp::export]]
DataFrame cav_summary(Rcpp::XPtr<VariationSet> xptr, DataFrame df, int min_cov, double max_seg_freq)
{
  Rcpp::XPtr<VariationSet> varset(xptr);
  StringVector contig_v = df["contig"];
  IntegerVector from_v = df["from"];
  IntegerVector to_v = df["to"];
  StringVector vtype_v = df["vtype"];
  
  vector<string> contig_r;
  vector<int> from_r;
  vector<int> to_r;
  vector<int> seg_r;
  vector<int> cov_r;

  for (unsigned int i=0; i<contig_v.size(); ++i) {
    string contig = as<string>(contig_v[i]);
    int from = from_v[i] - 1;
    int to = to_v[i] - 1;
    VariType vtype = Variation::string_to_type(as<string>(vtype_v[i]));
  
    int seg, cov;
    varset->get_summary(contig, vtype, from, to, min_cov, max_seg_freq, seg, cov);
    if (from < 0) {
      switch (from) {
      case -1: printf("contig %s not found\n", contig.c_str()); break;
      case -2: printf("query segment [%d,%d] is invalid\n", from_v[i], to_v[i]); break;
      }
      continue;
    }
    contig_r.push_back(contig);
    from_r.push_back(from+1);
    to_r.push_back(to+1);
    seg_r.push_back(seg);
    cov_r.push_back(cov);
  }

  DataFrame result = DataFrame::create(
				       Named("contig") = contig_r,
				       Named("from") = from_r,
				       Named("to") = to_r,
				       Named("seg.count") = seg_r,
				       Named("xcov") = cov_r);
  return result;
}

// [[Rcpp::export]]
DataFrame cav_cov_summary(Rcpp::XPtr<VariationSet> xptr, DataFrame df, int min_cov)
{
  Rcpp::XPtr<VariationSet> varset(xptr);
  StringVector contig_v = df["contig"];
  IntegerVector from_v = df["from"];
  IntegerVector to_v = df["to"];

  vector<string> contig_r;
  vector<int> from_r;
  vector<int> to_r;
  vector<int> cov_r;

  for (unsigned int i=0; i<contig_v.size(); ++i) {
    string contig = as<string>(contig_v[i]);
    int from = from_v[i] - 1;
    int to = to_v[i] - 1;
    int cov;
    varset->get_cov_summary(contig, from, to, min_cov, cov);
    if (from < 0) {
      switch (from) {
      case -1: printf("contig %s not found\n", contig.c_str()); break;
      case -2: printf("query segment [%d,%d] is invalid\n", from_v[i], to_v[i]); break;
      }
      continue;
    }
    contig_r.push_back(contig);
    from_r.push_back(from+1);
    to_r.push_back(to+1);
    cov_r.push_back(cov);
  }

  DataFrame result = DataFrame::create(
				       Named("contig") = contig_r,
				       Named("from") = from_r,
				       Named("to") = to_r,
				       Named("xcov") = cov_r);
  return result;
}

// [[Rcpp::export]]
DataFrame cav_freq_summary(Rcpp::XPtr<VariationSet> xptr, DataFrame df)
{
  Rcpp::XPtr<VariationSet> varset(xptr);
  StringVector contig_v = df["contig"];
  IntegerVector from_v = df["from"];
  IntegerVector to_v = df["to"];

  vector<string> contig_r;
  vector<int> from_r;
  vector<int> to_r;
  vector<int> count_r;
  vector<double> freq_r;

  for (unsigned int i=0; i<contig_v.size(); ++i) {
    string contig = as<string>(contig_v[i]);
    int from = from_v[i] - 1;
    int to = to_v[i] - 1;
    int count;
    double freq;
    varset->get_freq_summary(contig, from, to, count, freq);
    if (from < 0) {
      switch (from) {
      case -1: printf("contig %s not found\n", contig.c_str()); break;
      case -2: printf("query segment [%d,%d] is invalid\n", from_v[i], to_v[i]); break;
      }
      continue;
    }
    contig_r.push_back(contig);
    from_r.push_back(from+1);
    to_r.push_back(to+1);
    count_r.push_back(count);
    freq_r.push_back(freq);
  }

  DataFrame result = DataFrame::create(
				       Named("contig") = contig_r,
				       Named("from") = from_r,
				       Named("to") = to_r,
				       Named("count") = count_r,
				       Named("freq") = freq_r);
  return result;
}
