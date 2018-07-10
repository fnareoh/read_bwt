#include <sdsl/bit_vectors.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <tuple>
#include "read_bwt.hpp"
#include "internal/rle_string.hpp"

using namespace std;
using namespace sdsl;

read_bwt brut(){
  ifstream reads("data/reads.in", ios::in);
  if (reads){
    cout << "reads.in loaded." << endl;
    int N; int len_ext; int min_overlap;
    stringstream s;
    reads >> N >> len_ext >> min_overlap;
    for(int i = 0; i < N; i ++) {
      string r; reads >> r;
      s << r << "$";
    }
    string t = s.str();
    //Bit vector marking where the dolars are.
    bit_vector is_ext = bit_vector(t.size(), 0);
    for(int i = 0; i < t.size(); i++){
      if(t[i] == '$') is_ext[i] = 1;
    }
    return read_bwt(t, t.size(), is_ext);
  }
  else cout << "Can't open data/reads_ext.in." << endl;
  return read_bwt();
}

read_bwt format(string file_reads, string file){
  ifstream reads_ext(file_reads, ios::in);
  if (reads_ext){
    cout << file_reads << " loaded." << endl;
    int N, l, o;
    stringstream s;
    reads_ext >> N >> l;
    int len_s = 0;
    vector<int> extension_start(N, 0);
    vector<int> extension_length(N, 0);
    for(int i = 0; i < N; i ++){
      string r,e;
      reads_ext >> r >> e;
      extension_start[i] = len_s + r.size();
      extension_length[i] = e.size();
      s << r << e;
      len_s = len_s + r.size() + e.size();
    }
    string t = s.str();
    int n = t.size();
    cout << "Reads concatenation is loaded." << endl;
    //Bit vector marking where the extensions are.
    bit_vector is_ext = bit_vector(n, 0);
    for(int i = 0; i < N; i++){
      for(int j = 0; j < extension_length[i]; j++){
        is_ext[extension_start[i]+j] = 1;
      }
    }
		cout << "Bit vector is_ext is saved" << endl;
    reads_ext.close();
		return read_bwt(t,t.size(),is_ext);
  }
  else cout << "Can't open " << file_reads << "." << endl;

  return read_bwt();
}


int main() {

  read_bwt no_ext = brut();
  float size_no_ext = no_ext.save(".out");
  no_ext.load(".out");
  cout << "Brut done." << endl;
  cout << "Size of the construction without any extensions: " << size_no_ext << endl;
  read_bwt perfect = format("data/reads_ext.in", "_ext.out");
  float size_perfect = perfect.save("_ext.out");
  perfect.load("_ext.out");
  cout << "Extracted extension done." << endl;
  cout << "Size of the construction with perfect extensions: " << size_perfect << endl;
  read_bwt gen = format("data/reads_ext_gen.in", "_gen.out");
  float size_gen = gen.save("_gen.out");
  gen.load("_gen.out");
  cout << "Generated extension done." << endl;
  cout << "Size of the construction with gen extensions: " << size_gen << endl;
  read_bwt new_gen = format("data/reads_ext_new_gen.in", "_new_gen.out");
  float size_new_gen = new_gen.save("_new_gen.out");
  new_gen.load("_new_gen.out");
  cout << "New way generated extension done." << endl;
  cout << "Size of the construction with new_gen extensions: " << size_new_gen << endl;

	ifstream pattern("data/pattern.in", ios::in);
  if (pattern){
    int nb_pattern, pattern_length;
    while (pattern >> nb_pattern >> pattern_length){
      int nb_found = 0; int nb_found_ext = 0; int nb_found_gen = 0; int nb_found_new_gen = 0;
      for(int i = 0; i < nb_pattern; i++){
        string pat; pattern >> pat;
        nb_found += no_ext.search_pattern(pat);
        nb_found_ext += perfect.search_pattern(pat);
        nb_found_gen += gen.search_pattern(pat);
        nb_found_new_gen += new_gen.search_pattern(pat);

      }
      cout << "For " << nb_pattern << " patterns of size " << pattern_length << "." << endl;
      cout << nb_found << " match where found in the reads without extension." << endl;
      cout << nb_found_ext << " match where found in the reads with perfect extension." << endl;
      cout << nb_found_gen << " match where found in the reads with generated extension." << endl;
      cout << nb_found_new_gen << " match where found in the reads with new generated extension." << endl;
      if (nb_found_ext != 0) cout << (float)(nb_found_ext-nb_found) / (float) nb_found_ext * 100 << "% of the matches are false positive with perfect extension." << endl;
      if (nb_found_gen != 0) cout << (float)(nb_found_gen-nb_found) / (float) nb_found_gen * 100 << "% of the matches are false positive with generated extension." << endl;
      if (nb_found_new_gen != 0) cout << (float)(nb_found_new_gen-nb_found) / (float) nb_found_new_gen * 100 << "% of the matches are false positive with newly generated extension." << endl;
    }
  }
  return 0;
}
