#ifndef __VARIATION__
#define __VARIATION__

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

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Variation
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*
vtSubstitute: Single nucleotide replaced on coordinate
vtDelete: Nucleotides deleted starting from coordinate
vtInsert: Nucleotides inserted to left of coordinate
vtDangleLeft: Read dangles to the left of coordinate
vtDangleRight: Read dangles the right of coordinate
*/

enum VariType { vtNone, vtSubstitute, vtDelete, vtInsert, vtDangleLeft, vtDangleRight, vtCount };
class Variation {
 private:
  int m_type_bitmap;

  // vtSubstitute data
  char m_sub_nt;

  // vtInsert data
  string m_insert_seq;

  // vtDelete data
  int m_delete_length;

  // vtDangleLeft/vtDangleRight data
  string m_adj_contig;
  bool m_adj_strand;
  
  void from_string(const string str);

 public:

  // ctor
  Variation();
  Variation(string str);

  void save(ofstream& out);
  void load(ifstream& in);

  void add_sub(string nt);
  void add_delete(int length);
  void add_insert(string seq);
  void add_dangle_left(string contig, bool strand);
  void add_dangle_right(string contig, bool strand);

  void set_type(VariType type);
  void get_types(vector<VariType>& types) const;
  bool has_type(VariType type);
  
  // a variation
  bool is_ref() { return m_type_bitmap == 0; };

  string to_string() const;

  friend bool operator<(const Variation& lhs, const Variation& rhs);
  friend bool operator==(const Variation& lhs, const Variation& rhs);
  friend ostream& operator<<(ostream& os, const Variation& var);

  ////////////////////////////////////////////////////////////////////
  // static
  ////////////////////////////////////////////////////////////////////

  // type to/from string
  static string type_to_string(const VariType type);
  static VariType string_to_type(const string str);
};

struct VariationKey {
  string id;
  string contig;
  int coord;
  Variation var;
  VariationKey(string _id, string _contig, int _coord, Variation _var) :
    id(_id), contig(_contig), coord(_coord), var(_var) {};
  friend bool operator<(const VariationKey& lhs, const VariationKey& rhs);
};

struct Segment {
  string id;
  string contig;
  int start;
  int end;
  Segment(string _id, string _contig, int _start, int _end) :
    id(_id), contig(_contig), start(_start), end(_end) {};
  friend bool operator<(const Segment& lhs, const Segment& rhs);
  friend bool operator==(const Segment& lhs, const Segment& rhs);
};

#endif
