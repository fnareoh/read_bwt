#include <iostream>
#include <algorithm>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/lcp.hpp>
#include <sstream>

using namespace std;
using namespace sdsl;

/*
extension.cpp generates fake extensions from the raw reads.
It takes as input the file data/reads.in with the following format:
- first line: Number_of_reads length_of_the_extension minimum_overlap
- following lines: a read per line
And it outputs a file data/reads_ext_gen.in with the following format:
- first line: Number_of_reads length_of_the_extension
- following lines: read extension$
*/

/* Given the position of the end of the read (dolar) and the begining of the
	read (next_dolar) and the reverse concatenation of the reads, get_read
	returns the read */
string get_read(int dolar, int next_dolar, string &RT){
	string read = RT.substr(dolar, next_dolar - dolar);
	reverse(read.begin(),read.end());
	return read;
}

/* Given the position of the start of an extension (start) and the length of the
	extension (len_ext) and the reverse concatenation of the reads, get_extension
	returns the extension */
string get_extension(int start, int len_ext, string &RT){
	string ext = RT.substr(start-len_ext, len_ext);
	reverse(ext.begin(),ext.end());
	return ext;
}

/* Given the concatenation of the reads, the number of reads , the length of the
	extensions, and the minimum overlap, build_extensions outputs the reads with
	the generated extension in data/reads_ext_gen.in */
int build_extensions(string &T, int N, int len_ext, int min_overlap){
	/* Reverse T */
	string RT(T);
	reverse(RT.begin(), RT.end());

	/* Construction of the LongestComonPrefix and CompressedSuffixArray */
	cout << "CSA in construction..." << endl;
	csa_bitcompressed<> csa_RT;
	construct_im(csa_RT, RT, 1);
	cout << "CSA constructed" << endl;
	lcp_bitcompressed<> lcp_RT;
	cout << "LCP in construction..." << endl;
	construct_im(lcp_RT, RT, 1);
	cout << "LCP constructed" << endl;

 /* Construction of the D array that gives the distance to the next dolars */
 cout << "D array in construction..." << endl;
	bit_vector is_dolar = bit_vector(RT.size(),0);
	for(int i = 0; i < RT.size(); i++){
		if (RT[i] == '$') is_dolar[i]=1;
	}
	bit_vector::rank_1_type dolar_rank(&is_dolar);
	bit_vector::select_1_type dolar_select(&is_dolar);

	vector<int> D(csa_RT.size(),0);
	for(int i = 0; i < D.size(); i++){
		int r = dolar_rank(csa_RT[i]);
		D[i] = csa_RT[i] - dolar_select(((r == 0) ? 1 : r));
	}
	cout << "D array constructed" << endl;

	/* Construction of the extenstions */
	int nb_dolars_seen = 0;
	int last_pos = -1;
	int min_lcp_range = 0;
	vector<int> save_dolar(dolar_rank(RT.size()-1), -1);
	//saves the position of each end of read
	vector<int> save_forward_position(dolar_rank(RT.size()-1), -1);
	//saves the best position we can copy form is there is one on the foward swipe
	//over the LCP array
	vector<int> save_forward_overlap(dolar_rank(RT.size()-1), -1);
	//saves the overlap with the best position we can copy form is there is one on
	//the forward swipe over the LCP array
	cout << "Foward swipe ..." << endl;
	for(int i = 0; i < D.size(); i++){
		if (D[i] == 1){
			//if we reach the end the of a read
			if (min_lcp_range >= min_overlap){
				//if the overlap is sufficient
				save_forward_position[nb_dolars_seen] = csa_RT[last_pos];
				save_forward_overlap[nb_dolars_seen] = min_lcp_range;
				//we save the position of the extension
			}
			//else we don't
			save_dolar[nb_dolars_seen] = i;
			nb_dolars_seen++;
		}
		else if (D[i] > len_ext){
			//if the suffix has a potential extension of sufficient size
			last_pos = i;
			min_lcp_range = lcp_RT[i+1];
			//we save the overlap
		}
		min_lcp_range = min(min_lcp_range, int(lcp_RT[i+1]));
		//we update the overlap
	}
	cout << "Foward swipe done" << endl;



	cout << "Backward swipe ..." << endl;
	nb_dolars_seen = dolar_rank(RT.size()-1) -1;
	last_pos = 0;
	min_lcp_range = 0;
	vector<int> save_backward_position(dolar_rank(RT.size()-1), -1);
	//saves the best position we can copy form is there is one on the backward swipe
	//over the LCP array
	vector<int> save_backward_overlap(dolar_rank(RT.size()-1), -1);
	//saves the overlap with the best position we can copy form is there is one on
	//the backward swipe over the LCP array
	for(int i = D.size() - 1; i >= 0; i--){
		min_lcp_range = min(min_lcp_range, int(lcp_RT[i+1]));
		//we update the overlap
		if (D[i] == 1){
			//if we reach the end the of a read
			if (min_lcp_range >= min_overlap){
				//if the overlap is sufficient
				save_backward_position[nb_dolars_seen] = csa_RT[last_pos];
				save_backward_overlap[nb_dolars_seen] = min_lcp_range;
				//we save the position of the extension
			}
			//else we don't
			nb_dolars_seen--;
		}
		else if (D[i] > len_ext){
			//if the suffix has a potential extension of sufficient size
			last_pos = i;
			min_lcp_range = lcp_RT[i+1];
			//we save the overlap
		}
	}
	cout << "Backward swipe done" << endl;

	/* Saving the reads and extenstionsin reads_ext_gen.in*/
	cout << "Saving the reads and extensions to file..." << endl;
	ofstream text_with_extension_file("data/reads_ext_gen.in");
	if (text_with_extension_file.is_open()){
		text_with_extension_file << N << " " << len_ext << endl;
		for(int i = 0; i < save_forward_position.size(); i++){
			//get the original read
			string read = get_read(csa_RT[save_dolar[i]], dolar_select(dolar_rank(csa_RT[save_dolar[i]])+1),RT);
			string extension = "";
			if (save_forward_position[i] != -1 && save_backward_position[i] != -1){
				if(save_forward_overlap[i] < save_backward_overlap[i])
					extension = get_extension(save_backward_position[i], len_ext, RT);
				else
					extension = get_extension(save_forward_position[i], len_ext, RT);
			}
			else if (save_forward_position[i] != -1)
				extension = get_extension(save_forward_position[i], len_ext, RT);
			else if (save_backward_position[i] != -1)
				extension = get_extension(save_backward_position[i], len_ext, RT);
			text_with_extension_file << read << " " << extension << "$" << endl;
		}
		text_with_extension_file.close();
	}
	else cout << "Couldn't open data/reads_ext_gen.in." << endl;
	cout << "Done" << endl;
	return 0;
}

int main(){
	string t;
	int N, len_ext, min_overlap;
	ifstream reads("data/reads.in", ios::in);
  if (reads){
    stringstream s;
    reads >> N >> len_ext >> min_overlap;
    for(int i = 0; i < N; i ++) {
      string r; reads >> r;
      s << r << "$";
    }
		t = s.str();
		cout << "reads.in loaded." << endl;
		reads.close();
  }
  else cout << "Can't open data/reads.in." << endl;

	cout << "Generation of the extension..." << endl;
	build_extensions(t, N, len_ext, min_overlap);
	cout << "Generation of the extensions done." << endl;

	return 0;
}
