#include"qqNorm/qqNorm.h"
#include<iostream> 
#include"SparseMatrix/SparseMatrix.h"
#include"rankNormalize/rankNormalize.h"
#include <string>  
#include"qqNorm/caculateInterface.h"
using namespace std;

int main(int argc, const char ** argv)
{
	qqNorm qq;
	SparseMatrix sm;
	string path_read = argv[1];
	sm.readHDF5File(path_read);
	unordered_map<int, string> numToCell = sm.get_numToCell();
	int columnCount = sm.get_gene_count();
	rankNormalize rn(sm);
	rn.ranks(5);
	unsigned short *rank= rn.getRank();
	sm.set_rank(rank);
	unsigned short *test = sm.createCellVectorByName(numToCell[0]);

	qqNorm q;
	double *result = q.caculateTheoryQuantiles(test, columnCount);
	for (int i = 0; i < columnCount; i++)
	{
		cout << result[i] << " ";
	}
	cin.get();
	cin.get();
	return 0;
}