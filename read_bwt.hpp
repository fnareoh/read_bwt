#include <iostream>
#include <sdsl/bit_vectors.hpp>
#include "wt.hpp"

using namespace std;
using namespace sdsl;

class read_bwt {
	int search_pattern_aux(string , int, int, int);
public:
	wavelet_F wt_F;
	wavelet_L wt_L;
	read_bwt();
	read_bwt(string &, int, bit_vector &);
	float save(string);
	int load(string);
	int search_pattern(string);
	int LF(int);
};

#include "read_bwt.cpp"
