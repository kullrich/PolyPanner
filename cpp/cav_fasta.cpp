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
#include <cstdlib>

#include <unistd.h>

#include "util.h"
#include "Params.h"
#include "VariationSet.h"

#include "ClusterVariants.h"
#include "Resource.h"

typedef map< string, set< VariationKey > >::iterator VkIterator;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// util functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void fasta_read_segments(string fn,
			 map<string, string>& contig_fasta,
			 map<string, map < string, string > >& bin_fasta)
{
  cout << "reading table: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  map<string, int> bin_to_index;
  
  // parse header
  split_line(in, fields, delim);
  int id_ind = get_field_index("segment", fields);
  int contig_ind = get_field_index("contig", fields);
  int start_ind = get_field_index("start", fields);
  int end_ind = get_field_index("end", fields);
  int bin_ind = get_field_index("bin", fields);
  
  string name = "";
  string seq = "";
  string prev_contig = "";
  string prev_bin = "";
  int prev_end = -1;
  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;
    string id = fields[id_ind];
    string contig = fields[contig_ind];
    int start = stoi(fields[start_ind])-1;
    int end = stoi(fields[end_ind])-1;
    string bin = fields[bin_ind];
    int len = end - start + 1;

    massert(contig_fasta.find(contig) !=contig_fasta.end(), "contig %s not found in fasta file", contig.c_str());
    string seq_t = contig_fasta[contig].substr(start, len);
    
    bool connected = (bin == prev_bin) && (contig == prev_contig) && (start == prev_end);
    // cout << bin << " " << contig << ":" << start << "-" << end << " C=" << (connected ? "T" : "F") << endl;
    if (connected) {
      seq += seq_t;
      name += "_" + id;
    } else {
      if (seq != "")
	bin_fasta[prev_bin][name] = seq;
      seq = seq_t;
      name = id;
    }
    prev_bin = bin;
    prev_contig = contig;
    prev_end = end;
  }
  // last segment
  if (seq != "")
    bin_fasta[prev_bin][name] = seq;
}

void fasta_read_bins(string fn, map<string, int>& bins)
{
  cout << "reading table: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  map<string, int> bin_to_index;
  
  // parse header
  split_line(in, fields, delim);
  int bin_ind = get_field_index("bin", fields);
  int length_ind = get_field_index("length", fields);

  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;
    string bin = fields[bin_ind];
    int length = stoi(fields[length_ind]);
    bins[bin] = length;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void fasta_init_params(const char* name, int argc, char **argv, Parameters& params)
{
  params.add_parser("ifn_fasta", new ParserFilename("input contig fasta file"), true);
  params.add_parser("ifn_bins", new ParserFilename("input bin table"), true);
  params.add_parser("ifn_segments", new ParserFilename("input segment file, with bin field"), true);
  params.add_parser("min_length", new ParserInteger("minimal bin length", 10000), false);
  params.add_parser("odir", new ParserFilename("output path"), false);
  params.add_parser("ofn", new ParserFilename("output file"), false);
  
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

int fasta_main(const char* name, int argc, char **argv)
{
  Parameters params;
  fasta_init_params(name, argc, argv, params);

  string ifn_fasta = params.get_string("ifn_fasta");
  string ifn_bins = params.get_string("ifn_bins");
  string ifn_segments = params.get_string("ifn_segments");
  int min_length = params.get_int("min_length");
  string odir = "";
  string ofn = "";
  
  massert(params.is_defined("odir") || params.is_defined("ofn"), "ofn or odir must be defined");

  bool output_to_dir = params.is_defined("odir");
  bool output_to_file = params.is_defined("ofn");

  if (output_to_dir)
    odir = params.get_string("odir");
  
  if (output_to_file)
    ofn = params.get_string("ofn");
    
  // read bin table (bin -> length)
  map<string, int> bins;
  fasta_read_bins(ifn_bins, bins);

  // read contig fasta
  map<string, string> contig_fasta;
  load_fasta(ifn_fasta, contig_fasta);
  
  // bin fastas: (bin -> segment -> sequence)
  map<string, map < string, string > > bin_fasta;
  fasta_read_segments(ifn_segments, contig_fasta, bin_fasta);

  // output each bin in separate fasta file
  if (output_to_dir) {
    cout << "saving bin fasta files to " << odir << endl;
    for (map<string, int>::iterator it=bins.begin(); it!=bins.end(); ++it) {
      string bin = (*it).first;
      int length = (*it).second;
      if (length < min_length)
	continue;
      string ofn_bin = odir + "/" + bin + ".fa";
      massert(bin_fasta.find(bin) != bin_fasta.end(), "bin %s not found", bin.c_str());
      save_fasta(ofn_bin, bin_fasta[bin]);
    }
  }
	       
  // output all bins in single fasta file
  if (output_to_file) {
    map<string, string> single_fasta;
    for (map<string, int>::iterator it=bins.begin(); it!=bins.end(); ++it) {
      string bin = (*it).first;
      int length = (*it).second;
      if (length < min_length)
	continue;
      massert(bin_fasta.find(bin) != bin_fasta.end(), "bin %s not found", bin.c_str());

      // merge segments sequences
      map<string, string>& segment_fasta = bin_fasta[bin];
      string bin_seq = "";
      for (map<string, string>::iterator it=segment_fasta.begin(); it!=segment_fasta.end(); ++it) {
	string segment = (*it).first;
	string seq = (*it).second;
	bin_seq += seq;
      }
      single_fasta[bin] = bin_seq;
    }
    save_fasta(ofn, single_fasta);
  }
  
  return 0;
}
