#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include "H5Cpp.h"
#include<iostream>
#include<unordered_map>
using namespace std;
class HDF5reader
{
private:
	int *data;
	long long *indptr,*indices;
	char** barcodes;
	char** gene_names;
	int *startPos;
	unordered_map<int, string> numToCell;
	unordered_map<string, int> cellToNum;
	int gene_count;
	int cell_count, data_count;

public:
	HDF5reader() {
	}
	~HDF5reader()
	{
		delete[] indices;
		delete[] data;
		delete[] indptr;
		delete[] startPos;
		for (int i = 0; i < cell_count; i++)
		{
			delete[] barcodes[i];
		}
		delete[] barcodes;
		for (int i = 0; i < gene_count; i++)
		{
			delete[] gene_names[i];
		}
		delete[] gene_names;
	}
	char** get_barcodes();
	char** get_gene_names();
	long long* get_indptr();
	int readHDF5File(string path);
	void createCellnameMap();
	int* createCellVectorByName(string cellname);
};