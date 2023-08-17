#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <fstream>
#include <assert.h>
#include <sstream>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "util.h"
#include "Params.h"
#include "VariationSet.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void dump_init_params(const char* name, int argc, char **argv, Parameters& params)
{
  params.add_parser("ifn", new ParserFilename("POP filename"), true);
  params.add_parser("ofn_cov_cs", new ParserFilename("coverage cumsum output file"), true);
  params.add_parser("ofn_cov", new ParserFilename("coverage output file"), true);
  params.add_parser("ofn_cav", new ParserFilename("POP output file"), true);

  if (argc == 1) {
    params.usage(name);
    exit(1);
  }

  // read command line params
  params.read(argc, argv);
  params.parse();
  params.verify_mandatory();
  params.print(cout);
}


void dump_cov(string fn, map<string, vector<int> >& vv)
{
  cout << "saving vector to table: " << fn << endl;
  ofstream out(fn.c_str(), ios::out | ios::binary);
  massert(out.is_open(), "could not open file %s", fn.c_str());

  out << "contig" << "\t" << "coord" << "\t" << "count"  << endl;
  for (map<string, vector<int> >::iterator it=vv.begin(); it!=vv.end(); ++it) {
    string contig = (*it).first;
    vector<int>& coverage_contig = (*it).second;
    for (unsigned int i=0; i<coverage_contig.size(); ++i)
      out << contig << "\t" << i+1 << "\t" << coverage_contig[i] << endl;
  }
  out.close();
}

void dump_cav(VariationSet& varset, string fn)
{
  map< string, map< int, map <Variation, int> > >& vars = varset.get_vars();

  cout << "saving cav to table: " << fn << endl;
  ofstream out(fn.c_str(), ios::out | ios::binary);
  massert(out.is_open(), "could not open file %s", fn.c_str());

  // header
  out << "contig" << "\t" << "coord" << "\t" << "var" << "\t" << "count"  << endl;

  for (map< string, map< int, map <Variation, int> > >::iterator it=vars.begin(); it!=vars.end(); ++it) {
    string contig = (*it).first;
    map< int, map <Variation, int> >& table_contig = (*it).second;

    for (map< int, map <Variation, int> >::iterator jt=table_contig.begin(); jt != table_contig.end(); ++jt) {
      int coord = (*jt).first;
      map <Variation, int>& xmap = (*jt).second;
      for (map <Variation, int>::iterator xt=xmap.begin(); xt!=xmap.end(); ++xt) {
	Variation var = (*xt).first;
	int count = (*xt).second;
	out << contig << "\t" << coord+1 << "\t" << var.to_string() << "\t" << count  << endl;
      }
    }
  }
  out.close();
}

int dump_main(const char* name, int argc, char **argv)
{
  Parameters params;
  dump_init_params(name, argc, argv, params);

  string ifn = params.get_string("ifn");
  string ofn_cov = params.get_string("ofn_cov");
  string ofn_cav = params.get_string("ofn_cav");
  string ofn_cov_cs = params.get_string("ofn_cov_cs");

  VariationSet varset;
  varset.load(ifn);

  cout << "number of reads: " << varset.get_read_count() << endl;
  cout << "max read length: " << varset.get_max_read_length() << endl;
  if (ofn_cov != "NA")
    dump_cov(ofn_cov, varset.get_covs());
  if (ofn_cov_cs != "NA")
    dump_cov(ofn_cov_cs, varset.get_covs_cs());
  if (ofn_cav != "NA")
    dump_cav(varset, ofn_cav);

  return 0;
}
