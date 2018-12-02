#include "HDF5Reader/HDF5Reader.h"
#include<iostream>
using namespace std;

int main(int argc, const char ** argv)
{
	int gene_count;
	int cell_count;
	int data_count;
	HDF5reader hr;
	string path = argv[1];
	hr.readHDF5File(path);
	cout << "there are " << hr.get_gene_count() << " genes" << endl;
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
	hr.createCellnameMap();
	int* paraCellVector1=hr.createCellVectorByName(a[0]);

	
	cin.get();
	cin.get();
	return 0;
}