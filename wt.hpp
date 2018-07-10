#include <sdsl/bit_vectors.hpp>
#include <vector>
#include <iostream>

#include "internal/rle_string.hpp"


using namespace ri;
using namespace std;
using namespace sdsl;

class wavelet_F {
	std::vector<int> C;
	sd_vector<> end_read;
	sd_vector<>::rank_1_type end_read_rank1;
	sd_vector<>::rank_0_type end_read_rank0;
	sd_vector<>::select_0_type end_read_select0;
	sd_vector<>::select_1_type end_read_select1;
public:
	wavelet_F();
	wavelet_F(vector<int>,bit_vector);
	int select(int,char);
	int rank(int,char);
	float save(string);
	int load(string);
	int size(){return end_read.size();};
	int copy(const wavelet_F &);
	uchar operator[](ulint);
};


class wavelet_L {
	rle_string<> bwt;
	sd_vector<> a;
	sd_vector<>::rank_1_type a_rank;
	sd_vector<> c;
	sd_vector<>::rank_1_type c_rank;
	sd_vector<> g;
	sd_vector<>::rank_1_type g_rank;
	sd_vector<> t;
	sd_vector<>::rank_1_type t_rank;
public:
	wavelet_L();
	wavelet_L(string, sd_vector<>, sd_vector<>, sd_vector<>, sd_vector<>);
	int rank(int,char);
	float save(string);
	int load(string);
	int size(){return bwt.size();};
	int copy(const wavelet_L &);
	uchar operator[](ulint);
};

#include "wt.cpp"
