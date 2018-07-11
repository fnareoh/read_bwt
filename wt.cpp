
using namespace std;
using namespace sdsl;

wavelet_F::wavelet_F(){}

wavelet_F::wavelet_F(vector<int> C_, bit_vector end_read_){
	C = C_;
	end_read = sd_vector<>(end_read_);
	end_read_rank1 = sd_vector<>::rank_1_type(&end_read);
	end_read_rank0 = sd_vector<>::rank_0_type(&end_read);
	end_read_select0 = sd_vector<>::select_0_type(&end_read);
	end_read_select1 = sd_vector<>::select_1_type(&end_read);
}

char int_for_bp(char a){
  switch (a) {
    case 'A': return 0;
    case 'C': return 2;
    case 'G': return 4;
		case 'N': return 6;
		case 'T': return 8;
  }
  return -1;
}

int wavelet_F::select(int i, char X){
	int intX = int_for_bp(X);
	if (intX%2 == 0){
		return end_read_select0(end_read_rank0(C[intX/2])+i);
	}
	else {
		return end_read_select1(end_read_rank1(C[intX/2])+i);
	}
}


int wavelet_F::rank(int i, char X){
	int intX = int_for_bp(X);
	if (intX%2 == 0)
		return end_read_rank0(i)-end_read_rank0(C[intX/2]);
	else
		return end_read_rank1(i)-end_read_rank1(C[intX/2]);
}

int store_vector_C(vector<int> C, string file_C){
	ofstream C_file(file_C);
	if (C_file.is_open() ){
		C_file << C[0] << " " << C[1] << " " << C[2] << " " << C[3] << " " << C[4] << endl;
		C_file.close();
		return 0;
	}
	else cout << "Couldn't open " << file_C << "." << endl;
	return 1;
}

int load_vector_C(vector<int> &C, string file_C){
	ifstream C_file(file_C, ios::in);
	if (C_file){
		C_file >> C[0] >> C[1] >> C[2] >> C[3] >> C[4];
		C_file.close();
		return 0;
	}
	else cout << "Can't open " << file_C << " file" << endl;
	return 1;
}

float wavelet_F::save(string file){
	store_vector_C(C,"data/C"+file);
	store_to_file(end_read, "data/end_read"+file);
	return (float) size_in_mega_bytes(end_read);
}

int wavelet_F::load(string file){
	C = vector<int>(5,0);
	load_vector_C(C,"data/C"+file);
	load_from_file(end_read,"data/end_read"+file);
	end_read_rank1 = sd_vector<>::rank_1_type(&end_read);
	end_read_rank0 = sd_vector<>::rank_0_type(&end_read);
	end_read_select0 = sd_vector<>::select_0_type(&end_read);
	end_read_select1 = sd_vector<>::select_1_type(&end_read);
	return 0;
}

int wavelet_F::copy(const wavelet_F &r){
			C = r.C;
			end_read = r.end_read;
			end_read_rank1 = sd_vector<>::rank_1_type(&end_read);
			end_read_rank0 = sd_vector<>::rank_0_type(&end_read);
			end_read_select0 = sd_vector<>::select_0_type(&end_read);
			end_read_select1 = sd_vector<>::select_1_type(&end_read);
      return 0;
}

uchar wavelet_F::operator[](ulint i){
	if (i < C[1]) return 'A';
	if (i < C[2]) return 'C';
	if (i < C[3]) return 'G';
	if (i < C[4]) return 'N';
	else return 'T';
}

wavelet_L::wavelet_L(){}

wavelet_L::wavelet_L(string bwt_, sd_vector<> a_, sd_vector<> c_, sd_vector<> g_, sd_vector<> n_, sd_vector<> t_):
	bwt(bwt_), a(a_), c(c_), g(g_), n(n_), t(t_)
{
	a_rank = sd_vector<>::rank_1_type(&a);
	c_rank = sd_vector<>::rank_1_type(&c);
	g_rank = sd_vector<>::rank_1_type(&g);
	n_rank = sd_vector<>::rank_1_type(&n);
	t_rank = sd_vector<>::rank_1_type(&t);
}

char higher_case(char a){
  switch (a) {
    case 'a': return 'A';
    case 'c': return 'C';
		case 'g': return 'G';
    case 'n': return 'N';
    case 't': return 'T';
  }
  return '$';
}

int wavelet_L::rank(int i, char X){
	map<char, sd_vector<>::rank_1_type> bp_rank = {{'A',a_rank},{'C',c_rank},{'G',g_rank},{'N',n_rank},{'T',t_rank}};
	if (X - 'A' < 27) {
		return bwt.rank(i,X) - bp_rank[X](i);
	}
	return bp_rank[higher_case(X)](i);
}


float wavelet_L::save(string file){
	filebuf fb;
	if (fb.open("data/bwt"+file, ios::out)) {
		ostream bwt_file(&fb);
		float mb_size = (float) bwt.serialize(bwt_file)/ (float) 1000000;
		mb_size += size_in_mega_bytes(a)+size_in_mega_bytes(c)+size_in_mega_bytes(g)+size_in_mega_bytes(n)+size_in_mega_bytes(t);
		store_to_file(a, "data/a"+file);
		store_to_file(c, "data/c"+file);
		store_to_file(g, "data/g"+file);
		store_to_file(n, "data/n"+file);
		store_to_file(t, "data/t"+file);
		fb.close();
		return mb_size;
	}
	else cout << "Can't open " << "data/bwt"+file << " file" << endl;
	return 0;
}

int wavelet_L::load(string file){
	filebuf fb;
	if (fb.open("data/bwt"+file, ios::in)) {
		istream bwt_file(&fb);
		bwt.load(bwt_file);
		load_from_file(a, "data/a"+file);
		load_from_file(c, "data/c"+file);
		load_from_file(g, "data/g"+file);
		load_from_file(n, "data/n"+file);
		load_from_file(t, "data/t"+file);
		a_rank = sd_vector<>::rank_1_type(&a);
		c_rank = sd_vector<>::rank_1_type(&c);
		g_rank = sd_vector<>::rank_1_type(&g);
		n_rank = sd_vector<>::rank_1_type(&n);
		t_rank = sd_vector<>::rank_1_type(&t);
		fb.close();
	}
	return 0;
}

int wavelet_L::copy(const wavelet_L &r){
      wavelet_L res;
			bwt = r.bwt;
			a = r.a;
			c = r.c;
			g = r.g;
			n = r.n;
			t = r.t;
			a_rank = sd_vector<>::rank_1_type(&a);
			c_rank = sd_vector<>::rank_1_type(&c);
			g_rank = sd_vector<>::rank_1_type(&g);
			t_rank = sd_vector<>::rank_1_type(&t);
			n_rank = sd_vector<>::rank_1_type(&n);
			t_rank = sd_vector<>::rank_1_type(&t);
      return 0;
}

uchar wavelet_L::operator[](ulint i){
	assert(i<bwt.size());
	return bwt[i];
}
