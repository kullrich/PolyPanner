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

#include <queue>

#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

#include <dirent.h>

#include "util.h"
#include "Params.h"
#include "VariationSet.h"

#define massert_range(v, coord, id) massert(coord >= 0 && coord < (int)v.size(), "coord %d outside vector of length %d, id: %s", coord, (int)v.size(), id.c_str());
#define INCREMENT(v, coord, var, length) \
massert(coord >=0 && coord < length, "coord out of range, coord=%d, contig_size=%d", coord, length); \
v[coord][var]++;

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// I/O
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct MapFileFieldIndex {
  int id_ind, score_ind, contig_ind, coord_ind, back_coord_ind, strand_ind,
    edit_ind, sub_ind, insert_ind, delete_ind, cigar_ind, match_ind, unique_ind;
  MapFileFieldIndex(vector<string>& fields, string suffix) {
    id_ind = get_field_index("id", fields);
    score_ind = get_field_index("score"+suffix, fields);
    contig_ind = get_field_index("contig"+suffix, fields);
    coord_ind = get_field_index("coord"+suffix, fields);
    back_coord_ind = get_field_index("back_coord"+suffix, fields);
    strand_ind = get_field_index("strand"+suffix, fields);
    edit_ind = get_field_index("edit_dist"+suffix, fields);
    sub_ind = get_field_index("substitute"+suffix, fields);
    insert_ind = get_field_index("insert"+suffix, fields);
    delete_ind = get_field_index("delete"+suffix, fields);
    cigar_ind = get_field_index("cigar"+suffix, fields);
    match_ind = get_field_index("match_length"+suffix, fields);
    unique_ind = get_field_index("unique"+suffix, fields);
  }
};

int contruct_process_read_side(MapFileFieldIndex& fi, vector<string>& fields,
			       int offset_left, int offset_right, int trim,
			       string adj_contig, bool adj_strand,
			       map<string, int>& contig_map,
			       string discard_clipped, int min_score, int min_length, int max_edit,
			       bool unique_only,
			       VariationSet& varset,
			       string q_contig, int q_coord, ofstream& q_out)
{
  map< string, map< int, map <Variation, int> > >& vars = varset.get_vars();
  map<string, vector<int> >& covs = varset.get_covs();

  vector<string> sub_fields;
  char sub_delim = ';';

  vector<string> ssub_fields;
  char ssub_delim = ',';

  // use only longest match for each read
  string id = fields[fi.id_ind];
  string contig = fields[fi.contig_ind];
  int front_coord = atoi(fields[fi.coord_ind].c_str()) - 1;
  int back_coord = atoi(fields[fi.back_coord_ind].c_str()) - 1;
  bool strand = fields[fi.strand_ind] == "1";
  int edit_dist = atoi(fields[fi.edit_ind].c_str());
  int match_length = atoi(fields[fi.match_ind].c_str());
  int score = atoi(fields[fi.score_ind].c_str());
  string sub = fields[fi.sub_ind];
  string insert = fields[fi.insert_ind];
  string del = fields[fi.delete_ind];
  string cigar_str = fields[fi.cigar_ind];
  bool unique = fields[fi.unique_ind] == "T";
  
  int left_coord = strand ? back_coord : front_coord;
  int right_coord = !strand ? back_coord : front_coord;
  int read_length = right_coord - left_coord + 1;
    
  if (q_contig != "" && q_contig != contig)
    return 0;
    
  // skip if quality or match length fall under thresholds
  if (match_length < min_length || score < min_score || edit_dist > max_edit || (unique_only && !unique))
    return 0;

  // skip contigs which are not in contig table
  if(contig_map.find(contig) == contig_map.end())
    return 0;
  
  int contig_length = contig_map[contig];

  // parse cigar format to identify dangles
  vector<pair < char, int> > cigar;
  parse_cigar(cigar_str, cigar);
  bool clip_start = (cigar.front().first == 'S' || cigar.front().first == 'H');
  bool clip_end = (cigar.back().first == 'S' || cigar.back().first == 'H');
  bool clipped = clip_start || clip_end;

  // we trim from both sides requested number, subtracting the number of clipped sequences
  int trim_left = trim - (clip_start ? cigar.front().second : 0);
  int trim_right = trim - (clip_end ? cigar.back().second : 0);
  if (offset_left < trim_left)
    offset_left = trim_left;
  if (offset_right < trim_right)
    offset_right = trim_right;
  
  if (discard_clipped == "discard_any") {
    if (clip_start || clip_end)
      return 0;
  } else if (discard_clipped == "discard_both") {
    if (clip_start && clip_end)
      return 0;
  } else if (discard_clipped != "keep") {
    cout << "unknown value for discard_clipped field: " << discard_clipped << endl;
    exit(-1);
  }

  // coord->variation
  // establish what are all the variations the read contains
  map< int, Variation> read_vars;
  
  // dangle at beginning of read
  if (strand ? clip_start : clip_end) {
    int coord = strand ? left_coord : right_coord;
    if (strand)
      read_vars[coord].add_dangle_left("BACK", 1);
    else
      read_vars[coord].add_dangle_right("BACK", 1);
  }
  
  // dangle at end of read
  if (!strand ? clip_start : clip_end) {
    int coord = strand ? right_coord : left_coord;
    if (strand)
      read_vars[coord].add_dangle_right(adj_contig, adj_strand);
    else
      read_vars[coord].add_dangle_left(adj_contig, adj_strand);
  }
  
  // substitutes
  split_string(sub, sub_fields, sub_delim);
  for (unsigned int i=0; i<sub_fields.size(); ++i) {
    split_string(sub_fields[i], ssub_fields, ssub_delim);
    if (ssub_fields.size() < 4)
      continue;
    int coord = atoi(ssub_fields[0].c_str()) - 1;
    if (coord == -1) {
      cout << "id=" << id << ", contig=" << contig << ", left_coord=" << left_coord << ", right_coord=" << right_coord
	   << ", strand=" << (strand ? "+" : "-")
	   << ", cigar=" << cigar_str << ", sub=" << sub << ", insert=" << insert << ", delete=" << del
	   << ", clipped=" << clipped << ", read_length=" << read_length << ", match_length=" << match_length << endl;
    }
    
    massert(coord >= left_coord && coord <= right_coord, "substitution, coord %d out of range, left=%d, right=%d, id=%s\n",
	    coord, left_coord, right_coord, id.c_str());
    // string nt = (strand ? ssub_fields[3] : reverse_complement(ssub_fields[3]));
    string nt = ssub_fields[3];
    if (nt == "N")
      continue;
    read_vars[coord].add_sub(nt);
  }

  // deletions
  split_string(del, sub_fields, sub_delim);
  for (unsigned int i=0; i<sub_fields.size(); ++i) {
    split_string(sub_fields[i], ssub_fields, ssub_delim);
    if (ssub_fields.size() < 2)
      continue;
    int start_coord = atoi(ssub_fields[0].c_str()) - 1;
    if (start_coord == -1) {
      cout << "id=" << id << ", contig=" << contig << ", left_coord=" << left_coord << ", right_coord=" << right_coord
	   << ", strand=" << (strand ? "+" : "-")
	   << ", cigar=" << cigar_str << ", sub=" << sub << ", insert=" << insert << ", delete=" << del
	   << ", clipped=" << clipped << ", read_length=" << read_length << ", match_length=" << match_length << endl;
    }
    massert(start_coord >= left_coord && start_coord <= right_coord, "deletion, coord %d out of range, left=%d, right=%d, id=%s\n",
	    start_coord, left_coord, right_coord, id.c_str());
    int del_length = atoi(ssub_fields[1].c_str());
    read_vars[start_coord].add_delete(del_length);
  }
  
  // insertions
  split_string(insert, sub_fields, sub_delim);
  for (unsigned int i=0; i<sub_fields.size(); ++i) {
    split_string(sub_fields[i], ssub_fields, ssub_delim);
    if (ssub_fields.size() < 3)
      continue;
    int coord = atoi(ssub_fields[0].c_str()) - 1;

    // insertion to the right of the last coord is omitted since the associated coord is not covered
    if (coord == (right_coord+1))
      continue;
    
    if (coord == -1) {
      cout << "id=" << id << ", contig=" << contig << ", left_coord=" << left_coord << ", right_coord=" << right_coord
	   << ", strand=" << (strand ? "+" : "-")
	   << ", cigar=" << cigar_str << ", sub=" << sub << ", insert=" << insert << ", delete=" << del
	   << ", clipped=" << clipped << ", read_length=" << read_length << ", match_length=" << match_length << endl;
    }
    massert(coord >= left_coord && coord <= right_coord, "insertion, coord %d out of range, left=%d, right=%d, id=%s\n",
	    coord, left_coord, right_coord, id.c_str());
    // string seq = (strand ? ssub_fields[2] : reverse_complement(ssub_fields[2]));
    string seq = ssub_fields[2];
    read_vars[coord].add_insert(seq);
  }
  
  // increment for all variations identified in read
  map< int, map <Variation, int> >& coord_map = vars[contig];
  for (map< int, Variation>::iterator it=read_vars.begin(); it != read_vars.end(); ++it) {
    int coord = (*it).first;
    int read_index = coord - left_coord;
    Variation var = (*it).second;
    if (read_index < offset_left || read_index >= (read_length - offset_right))
      continue;
    INCREMENT(coord_map, coord, var, contig_length);
  }

  // coverage along read
  vector<int>& coverage_contig = covs[contig];
  for (int read_index=0; read_index<read_length; ++read_index) {
    int coord = left_coord + read_index;
    if (read_index < offset_left || read_index >= (read_length - offset_right))
      continue;
    massert_range(coverage_contig, coord, id);
    coverage_contig[coord]++;

    if (q_contig != "") {
      if (q_contig == contig && coord == q_coord) {
	q_out << id << "\t" << contig << "\t" << left_coord+1 << "\t" << right_coord+1
	      << "\t" << (strand ? "+" : "-")
	      << "\t" << cigar_str << "\t" << sub << "\t" << insert << "\t" << del
	      << "\t" << clipped << "\t" << read_length << "\t" << match_length << endl;
      }
    }
  }
  
  // append read
  varset.add_read(contig, left_coord+offset_left, read_length-(offset_left+offset_right));
  
  // one read side added
  return 1;
}

void construct_process_mapped_pair(string fn,
				   map<string, int>& contig_map,
				   string discard_clipped, int min_score, int min_length, int max_edit,
				   int trim, bool unique_only,
				   VariationSet& varset,
				   string q_contig, int q_coord, ofstream& q_out)
{
  cout << "reading fastq: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  
  vector<string> fields;
  char delim = '\t';

  // parse header
  split_line(in, fields, delim);

  // collect field indices for efficiency
  MapFileFieldIndex fi1(fields, "1"), fi2(fields, "2");

  int read_count = 1;
  int ok_read_count = 0;
  while(1) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;
    
    if (read_count++ % 1000000 == 0) {
      cout << "paired read: " << read_count-1 << " ..." << endl;
    }

    int offset_left1 = 0;
    int offset_right1 = 0;
    int offset_left2 = 0;
    int offset_right2 = 0;

    string id = fields[fi1.id_ind];
    string contig1 = fields[fi1.contig_ind];
    string contig2 = fields[fi2.contig_ind];

    // R1
    bool strand1 = fields[fi1.strand_ind] == "1";
    int front_coord1 = atoi(fields[fi1.coord_ind].c_str()) - 1;
    int back_coord1 = atoi(fields[fi1.back_coord_ind].c_str()) - 1;
      
    // R2
    bool strand2 = fields[fi2.strand_ind] == "1";
    int front_coord2 = atoi(fields[fi2.coord_ind].c_str()) - 1;
    int back_coord2 = atoi(fields[fi2.back_coord_ind].c_str()) - 1;

    int l1 = (strand1 ? back_coord1 : front_coord1);
    int r1 = (!strand1 ? back_coord1 : front_coord1);
    int l2 = (strand2 ? back_coord2 : front_coord2);
    int r2 = (!strand2 ? back_coord2 : front_coord2);
    bool skip_R1 = false;
    bool skip_R2 = false;
    
    // handle read overlap
    bool overlap_found = (contig1 == contig2) && (r1 > l2) && (r2 > l1);
    if (overlap_found) {
      int delta = 0;
      bool R1_on_left = (l1 < l2);
      if ((l1 >= l2) && (r1 <= r2)) {
	skip_R1 = true;
      } else if ((l2 >= l1) && (r2 <= r1)) {
	skip_R2 = true;
      } else if (R1_on_left)  {
	delta = r1-l2+1;
	offset_right1 += delta/2;
	offset_left2 += delta - delta/2;
	// cout << "delta=" << delta << ", offset_right1=" << offset_right1 << ", offset_left2=" << offset_left2 << endl;
	massert((r1 - offset_right1 + 1) == (l2 + offset_left2), "internal error computing offsets, R1 on left, %d != %d",
		r1 - offset_right1 + 1, l2 + offset_left2);
      } else {
	delta = r2-l1+1;
	offset_right2 += delta/2;
	offset_left1 += delta - delta/2;
	// cout << "delta=" << delta << ", offset_right2=" << offset_right2 << ", offset_left1=" << offset_left1 << endl;
	massert((r2 - offset_right2 + 1) == (l1 + offset_left1), "internal error computing offsets, R1 on right, %d != %d",
		r2 - offset_right2 + 1, l1 + offset_left1);
      }
    } 

    // if (!overlap_found) continue;

    if (!skip_R1)
      ok_read_count += contruct_process_read_side(fi1, fields,
						  offset_left1, offset_right1, trim,
						  contig2, strand2,
						  contig_map, discard_clipped, min_score, min_length, max_edit,
						  unique_only, varset,
						  q_contig, q_coord, q_out);
    
    // R2
    if (!skip_R2)
      ok_read_count += contruct_process_read_side(fi2, fields,
						  offset_left2, offset_right2, trim,
						  contig1, strand1,
						  contig_map, discard_clipped, min_score, min_length, max_edit,
						  unique_only, varset,
						  q_contig, q_coord, q_out);

    // if (ok_read_count > 10)
    //   break;
  }
  cout << "number of reads sides added to cav: " << ok_read_count << endl;
}

void construct_process_mapped_single(string fn,
				     map<string, int>& contig_map,
				     string discard_clipped, int min_score, int min_length, int max_edit,
				     int trim, bool unique_only,
				     VariationSet& varset,
				     string q_contig, int q_coord, ofstream& q_out)
{
  cout << "reading fastq: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  
  vector<string> fields;
  char delim = '\t';

  // parse header
  split_line(in, fields, delim);

  // collect field indices for efficiency
  MapFileFieldIndex fi(fields, "");

  int read_count = 1;
  int ok_read_count = 0;
  while(1) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;
    
    if (read_count++ % 1000000 == 0)
      cout << "paired read: " << read_count-1 << " ..." << endl;
    int added = contruct_process_read_side(fi, fields,
					   0, 0, trim,
					   "SINGLE", 1,
					   contig_map, discard_clipped, min_score, min_length, max_edit,
					   unique_only, varset,
					   q_contig, q_coord, q_out);
    ok_read_count += added;
  }
  cout << "number of reads sides added to cav: " << ok_read_count << endl;
}

void construct_read_contig_table(string fn, map<string, int>& contig_map, string q_contig)
{
  cout << "reading table: " << fn << endl;
  ifstream in(fn.c_str());
  massert(in.is_open(), "could not open file %s", fn.c_str());
  vector<string> fields;
  char delim = '\t';

  // parse header
  split_line(in, fields, delim);
  int contig_ind = get_field_index("contig", fields);
  int length_ind = get_field_index("length", fields);

  while(in) {
    split_line(in, fields, delim);
    if(fields.size() == 0)
      break;

    string contig = fields[contig_ind];
    int length = (int)atof(fields[length_ind].c_str());

    if (q_contig == "" || q_contig == contig)
      contig_map[contig] = length;
  }
}

inline double mround(double v, int digits) {
  double f = pow(10,digits);
  return (round(v * f) / f);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void construct_init_params(const char* name, int argc, char **argv, Parameters& params)
{
  params.add_parser("contig_table", new ParserFilename("contig table"), true);
  params.add_parser("ifn_paired", new ParserFilename("input file table of mapped paired reads"), true);
  params.add_parser("ifn_R1", new ParserFilename("input file table of mapped non-paired reads R1"), true);
  params.add_parser("ifn_R2", new ParserFilename("input file table of mapped non-paired reads R2"), true);
  params.add_parser("discard_clipped", new ParserString("handle clipped reads", "discard_both"), false);
  params.add_parser("trim", new ParserInteger("keep trim from read sides", 30), false);
  params.add_parser("min_length", new ParserInteger("minimal match length (nts)", 30), false);
  params.add_parser("min_score", new ParserInteger("minimal score", 30), false);
  params.add_parser("unique_only", new ParserBoolean("limit to unique hits", true), false);
  params.add_parser("max_edit", new ParserInteger("maximal edit distance", 2), false);
  params.add_parser("ofn", new ParserFilename("output variation table"), true);

  params.add_parser("query_contig", new ParserString("query contig (debug)", ""), false);
  params.add_parser("query_coord", new ParserInteger("query coordinate (debug)"), false);
  params.add_parser("query_ofn", new ParserFilename("output query filename (debug)"), false);
  
  if (argc == 1) {
    params.usage(name);
    cout << "discard_clipped options: " << endl;
    cout << "- discard_any: discard reads clipped on either side" << endl;
    cout << "- discard_both: discard reads clipped on both side" << endl;
    cout << "- keep: discard reads clipped on both side" << endl;
    exit(1);
  }

  // read command line params
  params.read(argc, argv);
  params.parse();
  params.verify_mandatory();
  params.print(cout);
}

int construct_main(const char* name, int argc, char **argv)
{
  Parameters params;
  construct_init_params(name, argc, argv, params);
  
  string ifn_paired = params.get_string("ifn_paired");
  string ifn_R1 = params.get_string("ifn_R1");
  string ifn_R2 = params.get_string("ifn_R2");
  string contig_table_ifn = params.get_string("contig_table");
  
  bool unique_only = params.get_bool("unique_only");
  string discard_clipped = params.get_string("discard_clipped");
  
  int min_score = params.get_int("min_score");
  int min_length = params.get_int("min_length");
  int max_edit = params.get_int("max_edit");
  int trim = params.get_int("trim");
  
  string ofn = params.get_string("ofn");

  string q_contig = params.get_string("query_contig");
  int q_coord = params.get_int("query_coord") - 1;
  string q_ofn = params.get_string("query_ofn");
  
  // contig lengths
  map<string, int> contig_map;
  construct_read_contig_table(contig_table_ifn, contig_map, q_contig);
  cout << "number of contigs: " << contig_map.size() << endl;

  VariationSet varset(contig_map);

  ofstream q_out;
  if (q_ofn != "") {
    q_out.open(q_ofn.c_str(), ios::out | ios::binary);
    cout << "saving query reads to table: " << q_ofn << endl;
    massert(q_out.is_open(), "could not open file %s", q_ofn.c_str());
    q_out << "id" << "\tcontig" << "\tleft_coord" << "\tright_coord" << "\tstrand" << "\tcigar" << "\tsub" << "\tinsert" << "\tdelete";
    q_out << "\tclipped" << "\tread_length" << "\tmatch_length" << endl;
  }

  // mapped paired reads
  construct_process_mapped_pair(ifn_paired, contig_map, discard_clipped,
				min_score, min_length, max_edit, trim, unique_only,  varset,
				q_contig, q_coord, q_out);

  // mapped non-paired reads
  construct_process_mapped_single(ifn_R1, contig_map, discard_clipped,
				  min_score, min_length, max_edit, trim, unique_only,  varset,
				  q_contig, q_coord, q_out);
  construct_process_mapped_single(ifn_R2, contig_map, discard_clipped,
				  min_score, min_length, max_edit, trim, unique_only,  varset,
				  q_contig, q_coord, q_out);

  varset.reads_done();

  if (q_ofn == "")
    varset.save(ofn);
  else
    q_out.close();

  return 0;
}
