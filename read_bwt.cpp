#include <sdsl/bit_vectors.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <iostream>
#include <vector>

using namespace std;
using namespace sdsl;


read_bwt::read_bwt(){}

read_bwt::read_bwt(string &T, int n, bit_vector &is_ext){
	csa_bitcompressed<> csa;
	cout << "Constructing the compressed suffix array..." << endl;
  construct_im(csa, T, 1);
  cout << "Compressed suffix array constructed." << endl;

  bit_vector::rank_0_type is_ext_rank0(&is_ext);
  int bwt_size = is_ext_rank0(is_ext.size());
	bit_vector a = bit_vector(bwt_size,0);
	bit_vector c = bit_vector(bwt_size,0);
  bit_vector g = bit_vector(bwt_size,0);
	bit_vector t = bit_vector(bwt_size,0);
  map<char, bit_vector> bp_vector = {{'A',a},{'C',c},{'G',g},{'T',t}};
  stringstream sbwt;
  vector<int> C(4,0);
  int C_index = 1;
  char last_letter_F = 'A';
  int sbwt_size = -1;
  bit_vector end_read = bit_vector(bwt_size,0);
  for(int i = 1; i < csa.size(); i++){
		//loop in the suffix array order
    int prev_i = (csa[i] - 1 + n) % n;
		//previous character in the text order
    if  ( T[prev_i]!='$' && not is_ext[prev_i] ) {
      sbwt <<  T[prev_i]; sbwt_size++;
      if (is_ext[csa[i]]) {
        end_read[sbwt_size] = 1;
      }
      else if (T[(csa[i] - 2 + n) % n]=='$') {
        bp_vector[T[prev_i]][sbwt_size] = 1;
      }
      if (T[csa[i]]!='$') {
        if (T[csa[i]]!=last_letter_F){
          C[C_index] = sbwt_size;
          C_index++;
        }
        last_letter_F = T[csa[i]];
      }
		}
	}
	wavelet_F wt_F_(C,end_read);
	wavelet_L wt_L_(sbwt.str(), sd_vector<>(bp_vector['A']), sd_vector<>(bp_vector['C']), sd_vector<>(bp_vector['G']), sd_vector<>(bp_vector['T']));
	wt_F.copy(wt_F_); //improper gestion
	wt_L.copy(wt_L_);
}

char lower_case(char a){
  switch (a) {
    case 'A': return 'a';
    case 'C': return 'c';
    case 'G': return 'g';
    case 'T': return 't';
  }
  return '$';
}

int read_bwt::search_pattern_aux(string pattern, int l, int i, int j){
  if (i > j || i < 0 || j > wt_L.size()) return 0;
  char X = pattern[l];
  int r1 = wt_L.rank(i, X) + 1;
  int r2 = wt_L.rank(j, X);
  if(r1 <= r2){
    int i2 = wt_F.select(((r1 == 0) ? 1 : r1) , X);
    int j2 = wt_F.select(r2 , X) + 1;
    if (l == pattern.size()-1) {
      int nb_pattern = max(j2 - i2,0) + wt_L.rank(j,lower_case(X)) - wt_L.rank(i,lower_case(X));
      return nb_pattern;
    }
    return search_pattern_aux(pattern, l+1, i2, j2);
  }
  else if (l == pattern.size()-1){
    int nb_pattern = wt_L.rank(j,lower_case(X)) - wt_L.rank(i,lower_case(X));
    return nb_pattern;
  }
  else
    return 0;
}

int read_bwt::search_pattern(string pattern){
	reverse(pattern.begin(), pattern.end());
	return search_pattern_aux(pattern, 0, 0, wt_L.size());
}

int read_bwt::LF(int i){
	if (wt_L[i] - 'A' >= 27) return -1;
	return wt_F.select(wt_L.rank(i,wt_L[i])+1,wt_L[i]);
}

float read_bwt::save(string file){
	return wt_F.save(file)+wt_L.save(file);
}

int read_bwt::load(string file){
	wt_L.load(file);
	wt_F.load(file);
	return 0;
}
