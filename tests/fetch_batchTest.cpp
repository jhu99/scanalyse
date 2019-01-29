#include<iostream> 
#include"SparseMatrix/SparseMatrix.h"
#include <string>  
using namespace std;

int main(int argc, const char ** argv)
{
	SparseMatrix sm;
	sm.readHDF5File(argv[1], argv[2]);
	int gene_count = sm.get_gene_count();
	double** inputMtrix = sm.fetch_batch(1,argv[2]);
	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 100; j++)
		{
			cout << inputMtrix[i][j] << " ";
		}
		cout << endl;
	}
	sm.deleteSparseMatrix(argv[2]);
	cin.get();
	cin.get();
	return 0;
}