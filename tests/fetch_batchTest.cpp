#include<iostream> 
#include"SparseMatrix/SparseMatrix.h"
#include <string>  
using namespace std;

int main(int argc, const char ** argv)
{
	SparseMatrix sm;
	string path_read = argv[1];
	string type = "qqnorm";
	sm.readHDF5File(path_read, type);
	int gene_count = sm.get_gene_count();
	double** inputMtrix = sm.fetch_batch(1);
	for (int i = 0; i < 128; i++)
	{
		for (int j = 0; j < 100; j++)
		{
			cout << inputMtrix[i][j] << " ";
		}
		cout << endl;
	}
	sm.deleteSparseMatrix();
	cin.get();
	cin.get();
	return 0;
}