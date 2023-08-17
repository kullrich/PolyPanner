#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <fstream>
#include <assert.h>
#include <sstream>
#include <stdarg.h>
#include <math.h>

#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

#include <dirent.h>

#include "util.h"
#include "Params.h"
#include "VariationSet.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// I/O
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_query_process_mapping_file(string fn, string q_contig, int q_start, int q_end, ofstream& out, bool& first_line)
{
  cout << "reading table: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  // parse header
  split_line(in, fields, delim);
  int id_ind = get_field_index("id", fields);
  int contig_ind = get_field_index("contig", fields);
  int coord_ind = get_field_index("coord", fields);
  int back_coord_ind = get_field_index("back_coord", fields);
  int strand_ind = get_field_index("strand", fields);

  // handle header
  if (first_line) {
    first_line = false;
    out << "fn";
    for (unsigned int i=0; i<fields.size(); ++i)
      out << "\t" << fields[i];
    out << endl;
  }
  
  int read_count = 1;
  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;

    if (read_count++ % 1000000 == 0) {
      cout << "progress: " << (read_count-1)/1000000 << "M reads ..." << endl;
    }
    
    // use only longest match for each read
    string id = fields[id_ind];
    string contig = fields[contig_ind];
    int front_coord = atoi(fields[coord_ind].c_str());
    int back_coord = atoi(fields[back_coord_ind].c_str());
    bool strand = fields[strand_ind] == "1";
    
    int left_coord = strand ? back_coord : front_coord;
    int right_coord = !strand ? back_coord : front_coord;

    if (q_contig != contig || right_coord < q_start || left_coord > q_end)
      continue;
    out << fn;
    for (unsigned int i=0; i<fields.size(); ++i)
      out << "\t" << fields[i];
    out << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_query_init_params(const char* name, int argc, char **argv, Parameters& params)
{
  params.add_parser("source_type", new ParserString("input source (dir | file_pair)", "file_pair"), true);
  params.add_parser("idir", new ParserFilename("input directory with parsed mapped reads (relevant when source=dir)"), false);
  params.add_parser("ifn_R1", new ParserFilename("input file R1 (relevant when source=file_pair)"), false);
  params.add_parser("ifn_R2", new ParserFilename("input file R1 (relevant when source=file_pair)"), false);
  params.add_parser("contig", new ParserString("query contig"), true);
  params.add_parser("start", new ParserInteger("query start coord"), true);
  params.add_parser("end", new ParserInteger("query end coord"), true);
  params.add_parser("ofn", new ParserFilename("output variation table"), true);
  
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

int read_query_main(const char* name, int argc, char **argv)
{
  Parameters params;
  read_query_init_params(name, argc, argv, params);
  
  string itype = params.get_string("source_type");
  string idir = params.get_string("idir");
  string ifn_R1 = params.get_string("ifn_R1");
  string ifn_R2 = params.get_string("ifn_R2");
  string contig = params.get_string("contig");
  int start = params.get_int("start");
  int end = params.get_int("end");
  string ofn = params.get_string("ofn");

  // collect files
  vector<string> ifns;
  if (itype == "dir") {
    cout << "number of fastq files found: " << get_of_files(idir) << endl;
    DIR *dir;
    struct dirent *ent;
    massert((dir = opendir (idir.c_str())) != NULL, "could not open input directory");
    // int file_count = 0;
    while ((ent = readdir (dir)) != NULL) {
      string ifn = ent->d_name;
      if (ifn.find("fastq") == string::npos || ifn.find("~") != string::npos)
	continue;
      ifns.push_back(idir + "/" + ifn);
    }
    closedir (dir);
  } else {
    ifns.push_back(ifn_R1);
    ifns.push_back(ifn_R2);
  }

  cout << "saving query reads to table: " << ofn << endl;
  ofstream out(ofn.c_str(), ios::out | ios::binary);
  massert(out.is_open(), "could not open file %s", ofn.c_str());
  
  // process files
  bool first_line = true;
  for (unsigned int i=0; i<ifns.size(); ++i) {
    if (ifns[i] == "NA")
      continue;
    read_query_process_mapping_file(ifns[i], contig, start, end, out, first_line);
  }

  return 0;
}
