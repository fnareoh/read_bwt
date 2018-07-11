#include <iostream>
#include <fstream>
#include <sstream>
#include <time.h>
#include <streambuf>
#include <algorithm>
#include <string>
#include <vector>
using namespace std;

/*
extract.cpp creates artificial reads and extension and patterns from a dna
squence.
It takes two files in input:
- data/dna.in where the dna sequence is.
- data/conf.in where the parameters for the generation of the reads and patterns are.
And produces two file in output:
- data/reads.in where the raw reads are stored.
- data/reads_ext.in where the reads with the perfect extension are stored.
- data/pattern.in where all the extracted patterns (used later for false-positive
  analysis) are stored.
*/


/*
filter_dna() returns the dna sequence from dna.in but whith only 'A' 'C' 'G' 'T'
characters.
*/
string filter_dna(){
  stringstream T;
  ifstream dna("data/dna.in", ios::in);
  if (dna)
  {
    string t((std::istreambuf_iterator<char>(dna)),
    std::istreambuf_iterator<char>());
    t.erase(std::remove(t.begin(), t.end(), '\n'), t.end());
    for(int i = 0; i < t.size();i++){
      if (t[i] == 'A' || t[i] == 'C' || t[i] == 'G' || t[i] == 'T' || t[i] == 'N'){
        T << t[i];
      }
    }
    dna.close();
  }
  else cout << "Can't oppen data/reads.in or data/reads_ext.in file" << endl;
  return T.str();
}


int main(){

  /* Recuperation of the generating parameters form conf.in */

  int S;//Memory space that can be used.
  int c;//Coverage (usally 30).
  int n;//Length of the reads (usally 200).
  int l;//Length of the extensions (usally 20).
  int N;//Numbre of reads to generate.
  int nb_pattern;//number of pattern we want to extract for each type of pattern.
  int nb_type_pattern;//number of type of pattern.
  vector<int> pattern_length;//length of the pattern we will generate.

  ifstream conf("data/conf.in", ios::in);
  if (conf){
    /* Recuperation of the generating parameters for reads and extensions */
    conf >>  S >> c >> n >> l;
    /* Recuperation of the generating parameters for patterns */
    conf >> nb_pattern >> nb_type_pattern;
    for(int i = 0; i < nb_type_pattern; i++){
      int length; conf >> length;
      pattern_length.push_back(length);
    }
    N = S/n;
    conf.close();
  }
  else  cout << "Can't open data/conf.in file" << endl;

  /* Generation of the reads and extensions */

  cout << "DNA sequence loading..." << endl;
  string t = filter_dna(); // Get dna sequence
  cout << "DNA sequence loaded." << endl;
  int part_dna_used = min(S/c,int(t.size())); // restrict the size in order to
  //limit the space usage and respect the 30 coverage.
  cout << "DNA file size: " << t.size() << endl;
  cout << "DNA file size we extract from: " << part_dna_used << endl;
  ofstream reads("data/reads.in");
  ofstream reads_ext("data/reads_ext.in");
  if (reads.is_open() && reads_ext.is_open()){
    cout << "Reads generation started..." << endl;
    /* We put at the begining of reads_ext the number of reads and the length
    of the extensions.*/
    reads_ext << N << " " << l << endl;
    /* We put at the begining of readsthe number of reads and the length
    of the extensions, and the minimum overlap (used in extension.cpp) temporary
    defined as the same length of the extension.*/
    reads << N << " " << l  << " " << l << endl;
    /* Random extraction of reads from our dna sample */
    srand (time(NULL));
    for(int i = 0; i < N; i ++){
      if(i % 10000 == 0) cout << i << " reads generated." << endl;
      int r = rand() % (part_dna_used-n-l); //only extract from part_dna_used
      // to keep the 30 coverage.
      string read = t.substr(r, n);
      string ext = t.substr(r+n, l);
      reads << read << endl;
      reads_ext << read << " " << ext << "$" << endl;
    }
    cout << "Reads generation done." << endl;
    reads.close(); reads_ext.close();
  }
  else cout << "Can't oppen data/reads.in or data/reads_ext.in file" << endl;


  /* Generation of the patterns */

  ofstream pattern("data/pattern.in");
  if (pattern.is_open()){
    /* Random extraction of patterns from our dna sample */
    srand (time(NULL));
    for(int j = 0; j < nb_type_pattern; j++){
      pattern << nb_pattern << " " << pattern_length[j] << endl;
      for(int i = 0; i < nb_pattern; i ++){
        int r = rand() % (part_dna_used-n-l);
        pattern << t.substr(r, pattern_length[j]) << endl;
      }
    }
    cout << "pattern generation done." << endl;
    pattern.close();
  }
  else cout << "Can't open data/pattern.in file" << endl;

  return 0;
}
