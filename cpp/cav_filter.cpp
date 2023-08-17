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
#include "Filter.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void filter_init_params(const char* name, int argc, char **argv, Parameters& params)
{
  params.add_parser("ifn", new ParserFilename("POP filename"), true);
  params.add_parser("threads", new ParserInteger("Number of threads used", 40), false);
  params.add_parser("alpha", new ParserDouble("False discovery rate", 0.001), false);
  params.add_parser("min_variant_count", new ParserInteger("Minimal count of considered variant", 5), false);
  params.add_parser("min_contig_length", new ParserInteger("Minimal contig length", 1000), false);
  params.add_parser("min_error_pct", new ParserDouble("Minimal error rate (%)", 0.01), false);
  params.add_parser("max_error_pct", new ParserDouble("Maximal error rate (%)", 2), false);
  params.add_parser("seed_error_pct", new ParserDouble("Seed error rate (%) used at start of optimization", 1), false);
  params.add_parser("ofn_sites", new ParserFilename("Output table with segregating sites"), true);
  params.add_parser("ofn_error_params", new ParserFilename("Output table with segregating sites"), true);

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

int filter_main(const char* name, int argc, char **argv)
{
  Parameters params;
  filter_init_params(name, argc, argv, params);

  string ifn = params.get_string("ifn");
  int thread_count = params.get_int("threads");
  double alpha = params.get_double("alpha");
  int min_cov = params.get_int("min_variant_count");
  int min_contig_length = params.get_int("min_contig_length");
  double min_error = params.get_double("min_error_pct") / 100;
  double max_error = params.get_double("max_error_pct") / 100;
  double seed_error = params.get_double("seed_error_pct") / 100;
  string ofn_sites = params.get_string("ofn_sites");
  string ofn_errors = params.get_string("ofn_error_params");

  // no need for a user-defined parameter as we expect to always converge quickly
  const int max_iterations = 100;
  
  VariationSet cav;
  cav.load(ifn);

  cout << "number of reads: " << cav.get_read_count() << endl;
  cout << "max read length: " << cav.get_max_read_length() << endl;

  Filter filter(cav, min_cov, min_contig_length, min_error, max_error, seed_error, alpha, thread_count);
  filter.resolve(max_iterations, ofn_sites);

  filter.save_error_rates(ofn_errors);
  filter.save_sites(ofn_sites);

  return 0;
}
