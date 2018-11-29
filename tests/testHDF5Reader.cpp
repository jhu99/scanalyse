#include "HDF5Reader.h"

#include<iostream>
using namespace std;
int main()
{
	int gene_count = 33694;
	int cell_count = 384000;
	int data_count = 260473471;
	HDF5reader hr(gene_count,cell_count,data_count);
	string path = "E:\\graduated\\Hdf5data\\ica_cord_blood_h5.h5";
	hr.readHDF5File(path);
	char**a = hr.get_barcodes();
	hr.createCellnameMap();
	int* paraCellVector=hr.createCellVectorByName(a[0]);
	cout << endl;
	
	/*for (int i = 0; i < 1000; i++)
	{
		cout << paraCellVector1[i]<<" ";
	}*/
	cin.get();
	cin.get();
	return 0;
}