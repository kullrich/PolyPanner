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
// combine functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void combine_read_library_table(string fn, string afield, map< string, string >& ifns)
{
  cout << "reading table: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  // parse header
  split_line(in, fields, delim);
  int aid_ind = get_field_index(afield, fields);
  int fn_ind = get_field_index("fn", fields);

  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;
    string aid = fields[aid_ind];
    string ifn = fields[fn_ind];
    ifns[aid] = ifn;
  }
}

int combine_read_contigs(string fn, string afield, map< string, set< string > >& contigs)
{
 cout << "reading binned segment table: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  // parse header
  split_line(in, fields, delim);
  int aid_ind = get_field_index(afield, fields);
  int contig_ind = get_field_index("contig", fields);

  int count = 0;
  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;

    string aid = fields[aid_ind];
    string contig = fields[contig_ind];
    contigs[aid].insert(contig);
    count++;
  }
  return count;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void combine_init_params(const char* name, int argc, char **argv, Parameters& params)
{
  params.add_parser("ifn_libs", new ParserFilename("table with POP files"), true);
  params.add_parser("ifn_contigs", new ParserFilename("table with contigs"), true);
  params.add_parser("assembly_field", new ParserFilename("assembly field"), true);
  params.add_parser("ofn", new ParserFilename("output site file"), true);

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

int combine_main(const char* name, int argc, char **argv)
{
  Parameters params;
  combine_init_params(name, argc, argv, params);

  string ifn_libs = params.get_string("ifn_libs");
  string ifn_contigs = params.get_string("ifn_contigs");
  string afield = params.get_string("assembly_field");
  string ofn = params.get_string("ofn");

  // aid -> contigs
  map< string, set< string > > contigs;
  int n_contigs = combine_read_contigs(ifn_contigs, afield, contigs);
  cout << "number of contigs: " << n_contigs << endl;

  // aid -> fn
  map< string, string > ifns;
  combine_read_library_table(ifn_libs, afield, ifns);
  cout << "number of assemblies: " << ifns.size() << endl;
  
  VariationSet rr;
  for (map< string, string >::iterator it = ifns.begin(); it != ifns.end(); ++it) {
    string aid = (*it).first;
    string fn = (*it).second;

    // skip entire assembly if no contigs associated with it
    if (contigs.find(aid) == contigs.end())
      continue;
    
    // load cav file
    VariationSet vs;
    vs.load(fn);

    set< string >& contigs_aid = contigs[aid];
    for (set< string >::iterator jt = contigs_aid.begin(); jt != contigs_aid.end(); ++jt) {
      string s_contig = *jt;
      string t_contig = aid + "_" + s_contig;

      rr.get_vars()[t_contig] = vs.get_vars()[s_contig];
      rr.get_covs()[t_contig] = vs.get_covs()[s_contig];
      rr.get_covs_cs()[t_contig] = vs.get_covs_cs()[s_contig];
      rr.get_covs_noise()[t_contig] = vs.get_covs_noise()[s_contig];
    }
  }

  cout << "final contig count: "<< rr.get_contigs().size() << endl;
  cout << "final first contig: "<< rr.get_contigs()[0] << endl;
  
  cout << "saving POP file: " << ofn << endl;
  rr.save(ofn);

  return 0;
}
