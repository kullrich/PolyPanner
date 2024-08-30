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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void info_init_params(const char* name, int argc, char **argv, Parameters& params)
{
  params.add_parser("ifn", new ParserFilename("POP filename"), true);

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

int info_main(const char* name, int argc, char **argv)
{
  Parameters params;
  info_init_params(name, argc, argv, params);

  string ifn = params.get_string("ifn");
  
  VariationSet cav;
  cav.load(ifn);

  cout << "read count: " << cav.get_read_count() << endl;
  cout << "max read length: " << cav.get_max_read_length() << endl;
  cout << "mean read length: " << cav.get_mean_read_length() << endl;
  
  vector<string> contigs = cav.get_contigs();
  cout << "contig count: " << contigs.size() << endl;

  int vcount, unique_count;
  cav.get_total_var_counts(vcount, unique_count);
  cout << "number of unique variants: " << unique_count << endl;
  cout << "number of reads supporting variants: " << vcount << endl;

  //  cout << "contigs:";
  // for (unsigned int i=0; i<contigs.size(); ++i)
  //   cout << " " << contigs[i];
  // cout << endl;
  
  return 0;
}
