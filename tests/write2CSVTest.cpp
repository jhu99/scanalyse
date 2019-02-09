#include"SparseMatrix/SparseMatrix.h"
using namespace std;

int main(int argc, const char ** argv) {
	SparseMatrix m;
	m.readHDF5File(argv[1], argv[3]);
	m.write2CSV(argv[2],argv[3]);
}
