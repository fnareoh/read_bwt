#include <iostream>
#include <iostream>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/lcp.hpp>
#include <vector>

using namespace std;
using namespace sdsl;

/*
all_extension.cpp generates fake extensions from the raw reads.
It takes as input the file data/reads.in with the following format:
- first line: Number_of_reads length_of_the_extension minimum_overlap
- following lines: a read per line
and conf_new_ext.in which only stores the number of iteration (nb_iteration).
And it outputs a file data/reads_ext_new_gen.in with the following format:
- first line: Number_of_reads length_of_the_extension
- following lines: read extension$
*/

/* Given the number of the read (r), the array containing the position of the
	start of each reads (start) and the concatenation of the reads,
	get_read returns the read */
string get_read(int r, vector<int> &start, string &t){
	string read = t.substr(start[r],start[r+1]-start[r]-1);
	return read;
}

/* Given the number of the read (r), the array containing the position of the
	start of each reads (start), the length of the extension (len), and the
	concatenation of the reads, get_extension returns the extension */
string get_extension(int r, vector<int> &start, int len, string &t){
	string ext = t.substr(start[r+1]-1-len, len);
	return ext;
}

/* A simple binary search */
int binarysearch(int i, int j, vector<int> &t, int s){
	int m = (i+j)/2;
	if (j-i == 1) return i;
	if (s >= t[m]) return binarysearch(m,j,t,s);
	else return binarysearch(i,m,t,s);
}

/* Given a position in the read concatenation, find_predecessor returns the
	number of the read */
int find_predecessor(int s, vector<int> &start_position_read){
	return binarysearch(0,start_position_read.size(), start_position_read, s);
}

/* first_iteration builds the first iteration, the source from which we copy
  the extension for each read */
tuple<vector<int>, vector<int>, vector<int>, vector<int>> first_iteration(string &T, string &RT, csa_bitcompressed<> &csa_RT, lcp_bitcompressed<> &lcp_RT, vector<int> &D, int N){
	/* Construction of the D array that gives the distance to the next dolars */
	vector<int> start_position_read(N+1,0);
	cout << "D array in construction..." << endl;
	bit_vector is_dolar = bit_vector(RT.size(),0);
	vector<int> end_position_read(N,0);
	int cpt_read = 0;
	for(int i = 0; i < T.size(); i++){
		if (T[i] == '$') {
			is_dolar[T.size()-1-i]=1;
			end_position_read[cpt_read] = i-1;
			cpt_read++;
			start_position_read[cpt_read] = i+1;
		}
	}
	bit_vector::rank_1_type dolar_rank(&is_dolar);
	bit_vector::select_1_type dolar_select(&is_dolar);

	/* Construction of the extensions */
	vector<int> source(N,0);
	//store the source read from which the extension is copied
	vector<int> length_overlap(N,0);
	//saves the length of the overlap with the source read
	vector<int> length_extension(N,0);
	//saves the length of the length of the extension
	for(int i = 0; i < D.size()-1; i++){
		if(D[i] == 1){
			//if position i is the end of a read.
			int r = find_predecessor(T.size() - csa_RT[i], start_position_read);
			//find the number of the read coresponding to the position
			if (lcp_RT[i] < lcp_RT[i+1]){
				//if the longest overlap is the next one in the LCP array
				source[r] = find_predecessor(T.size() - csa_RT[i+1], start_position_read);
				//find the number of the read coresponding to the position of the extension
				//save the source extension
				length_overlap[r] = lcp_RT[i+1];
				length_extension[r] = D[i+1]-1;
			}
			else {
				//if the longest overlap is the previous one in the LCP array
				source[r] = find_predecessor(T.size() - csa_RT[i-1], start_position_read);
				//find the number of the read coresponding to the position of the extension
				//save the source extension
				length_overlap[r] = lcp_RT[i];
				length_extension[r] = D[i-1]-1;
			}
		}
	}
	return make_tuple(start_position_read, source, length_overlap, length_extension);
}

/* Builds the extensions returns the text that needs to go in reads_ext_new_gen.in */
string build_extensions_new_gen(string &t, string &RT, csa_bitcompressed<> &csa_RT, lcp_bitcompressed<> &lcp_RT, vector<int> &D, int N, int nb_iteration, int len_max){
	/* Geting the first iteration */
	auto res = first_iteration(t, RT, csa_RT, lcp_RT, D, N);
	vector<int> start_position_read = get<0>(res);
	vector<int> source = get<1>(res);
	vector<int> length_overlap = get<2>(res);
	vector<int> length_extension = get<3>(res);

	/* Iterate to create the extension */
	stringstream t_ext;
	for(int i = 0; i < N; i++){
		t_ext << get_read(i,start_position_read, t) << " ";
		// get read
		int total_len_ext = 0;
		int current_source = i;
		// get extension
		for(int k = 0; k < nb_iteration; k++){
			string ext= "";
			if (length_extension[current_source] > 0 && total_len_ext < len_max) ext = get_extension(source[current_source],start_position_read,length_extension[current_source],t) ;
			t_ext << ext;
			total_len_ext += ext.size();
			current_source = source[current_source];
		}
		t_ext << "$" << endl;
	}
	return t_ext.str();
}


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
int build_extensions_gen(string &T, string &RT, csa_bitcompressed<> &csa_RT, lcp_bitcompressed<> &lcp_RT, vector<int> &D, int N, int len_ext, int min_overlap){
	bit_vector is_dolar = bit_vector(RT.size(),0);
 	for(int i = 0; i < RT.size(); i++){
 		if (RT[i] == '$') is_dolar[i]=1;
 	}
 	bit_vector::rank_1_type dolar_rank(&is_dolar);
 	bit_vector::select_1_type dolar_select(&is_dolar);

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
	string T;
	int N = 0; // number of read
	int len_ext, min_overlap;
	int nb_iteration;
	ifstream conf("data/conf_new_ext.in", ios::in);
	if  (conf) {
		conf >> nb_iteration;
		conf.close();
	}
	else cout << "Can't open data/conf_new_ext.in." << endl;

	ifstream reads("data/reads.in", ios::in);
  if (reads){
    cout << "reads.in loaded." << endl;
    stringstream s;
    reads >> N >> len_ext >> min_overlap;
    for(int i = 0; i < N; i ++) {
      string r; reads >> r;
      s << r << '$';
    }
    cout << "evrey parameter loaded, starting generation" << endl;
		T = s.str();
    reads.close();
  }
  else cout << "Can't open data/reads.in." << endl;

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


	cout << "Generation of the extension..." << endl;
	build_extensions_gen(T, RT, csa_RT, lcp_RT, D, N, len_ext, min_overlap);
	cout << "Generation of the extensions done." << endl;


	ofstream reads_ext_gen("data/reads_ext_new_gen.in");
	if (reads_ext_gen.is_open()){
		reads_ext_gen << N << " " <<  len_ext << endl << build_extensions_new_gen(T, RT, csa_RT, lcp_RT, D, N, nb_iteration, len_ext);
		reads_ext_gen.close();
	}
	else cout << "Can't open data/reads_ext_new_gen.in." << endl;

	return 0;
}
