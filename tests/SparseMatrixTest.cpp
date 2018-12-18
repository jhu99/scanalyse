#include"qqNorm/qqNorm.h"
#include<iostream> 
#include"SparseMatrix/SparseMatrix.h"
#include"rank/rankNormalize.h"
#include <string>  
#include"qqNorm/caculateInterface.h"
using namespace std;

int main(int argc, const char ** argv)
{
	qqNorm qq;
	SparseMatrix sm;
	string path_read = argv[1];
	sm.readHDF5File(path_read);
	rankNormalize rn(sm);
	rn.ranks(5);
	unsigned short *rank = rn.getRank();
	sm.set_rank(rank);
	sm.qqNormedData2HDF5Format();
	sm.write2HDF5(argv[2]);
	sm.deleteSparseMatrix();
	cin.get();
	cin.get();
	return 0;
}