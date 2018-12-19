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
	string type = "original";
	sm.readHDF5File(path_read, type);
	sm.createZeroPosPerCell();
	rankNormalize rn(sm);
	rn.ranks(5);
	unsigned short *rank = rn.getRank();
	sm.set_rank(rank);
	sm.qqNormedData2HDF5Format();
	sm.write2HDF5(argv[2]);
	string deleteType = "write_qqnorm_data";
	sm.deleteSparseMatrix(deleteType);
	cin.get();
	cin.get();
	return 0;
}