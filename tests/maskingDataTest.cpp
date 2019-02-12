#include"SparseMatrix/SparseMatrix.h"
#include<vector>
using namespace std;
using namespace H5;

int main(int argc, const char ** argv) {
	SparseMatrix m;
	m.readHDF5File(argv[1], "original");
	m.maskingData(argv[2], argv[3]);
}
