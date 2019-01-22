#include"SparseMatrix/SparseMatrix.h"
#include<vector>
using namespace std;
using namespace H5;

int main(int argc,const char** argv) {
	SparseMatrix m;
	vector<string> paths;
	for(int i=1;i<argc;i++){
	    paths.push_back(argv[i]);
	}
	m.mergeDate(paths);
	m.h5Compressed("./data/test.h5", "s", 1000, 1);
}
