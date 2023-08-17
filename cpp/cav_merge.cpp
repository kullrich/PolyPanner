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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void merge_init_params(const char* name, int argc, char **argv, string& ofn, vector<string>& ifns)
{
  if (argc <= 2) {
    fprintf(stderr, "usage: %s <ofn> [ifn1 ifn2 ...]\n", name);
    exit(1);
  }

  ofn = argv[1];
  int i = 2;
  while (i < argc) {
    ifns.push_back(argv[i++]);
  }
}

int merge_main(const char* name, int argc, char **argv)
{
  string ofn;
  vector<string> ifns;
  merge_init_params(name, argc, argv, ofn, ifns);

  VariationSet result;

  for (unsigned int i=0; i<ifns.size(); ++i) {
    VariationSet varset_i;
    varset_i.load(ifns[i]);
    // vector<string> cc = varset_i.get_contigs();
    // for (unsigned int j=0; j<cc.size(); ++j)
    //   cout  << cc[j] << "\t";
    // cout << endl;

    if (i == 0)
      result = varset_i;
    else
      result = result + varset_i;
  }
  result.save(ofn);

  return 0;
}
