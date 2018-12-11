#include "SparseMatrix/SparseMatrix.h"
#include<iostream>
using namespace std;

int main(int argc, const char ** argv)
{
	int gene_count;
	int cell_count;
	int data_count;
	SparseMatrix hr;
	string pathin = argv[1];
	string pathout = argv[2];
	hr.readHDF5File(pathin);
	cout << "there are " << hr.get_gene_count() << " genes" << endl;
	cout << "there are " << hr.get_data_count() << " data" << endl;
	cout << "Before filtration, there are " << hr.get_cell_count() << " cells." << endl;
	unordered_map<int,int>c=hr.cellFiltration();
	cout <<"After filtration, there are " << hr.get_cell_count()-c.size()<<" cells left." << endl;
	int count = 0;
	for (int i = 0; i < hr.get_cell_count(); i++)
	{
		if (c[i] == 1)
		{
			count++;
		}
	}
	cout << count << endl;
	char**a = hr.get_barcodes();
	int* data = hr.get_data();
	long long* indptr = hr.get_indptr();
	long long* indices = hr.get_indices();
	for (int i = 0; i < 10; i++)
	{
		cout << data[i] << " ";
	}
	cout << endl;
	for (int i = 0; i < 10; i++)
	{
		cout << indptr[i] << " ";
	}
	cout << endl;
	for (int i = 0; i < 10; i++)
	{
		cout << indices[i] << " ";
	}
	cout << endl;
	hr.createCellnameMap();
	int* paraCellVector1=hr.createCellVectorByName(a[0]);
	hr.write2HDF5(pathout);
	hr.deleteSparseMatrix();
	
	cin.get();
	cin.get();
	return 0;
}