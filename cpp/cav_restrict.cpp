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
#include "VariationSet.h"
#include "Params.h"

int read_site_table(string fn, map< string, map< int, set <Variation> > >& sites)
{
  cout << "reading table: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  // parse header
  split_line(in, fields, delim);
  int contig_ind = get_field_index("contig", fields);
  int coord_ind = get_field_index("coord", fields);
  int var_ind = get_field_index("variant", fields);

  int count = 0;
  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;

    string contig = fields[contig_ind];
    int coord = atoi(fields[coord_ind].c_str())-1;
    Variation var(fields[var_ind]);
    
    sites[contig][coord].insert(var);
    count++;
  }
  return count;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void restrict_init_params(const char* name, int argc, char **argv, Parameters& params)
{
  params.add_parser("ifn_cav", new ParserFilename("POP filename"), true);
  params.add_parser("ifn_sites", new ParserFilename("Table with variants"), true);
  params.add_parser("ofn", new ParserFilename("Output POP file"), true);

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

int restrict_main(const char* name, int argc, char **argv)
{
  Parameters params;
  restrict_init_params(name, argc, argv, params);

  string ifn_cav = params.get_string("ifn_cav");
  string ifn_sites = params.get_string("ifn_sites");
  string ofn = params.get_string("ofn");
  
  map< string, map< int, set <Variation> > > sites;
  int n = read_site_table(ifn_sites, sites);
  cout << "number of sites: " << n << endl;
  VariationSet cav;
  cav.load(ifn_cav);
  
  cav.restrict(sites);
  cav.save(ofn);
    
  return 0;
}
