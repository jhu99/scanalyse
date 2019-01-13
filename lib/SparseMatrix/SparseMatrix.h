#pragma once
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include "H5Cpp.h"
#include<iostream>
#include<unordered_map>
#include <algorithm>
#include "qqNorm/qqNorm.h"
using namespace std;
class SparseMatrix
{
private:
	int *data;
	long long *indptr, *indices;
	char** barcodes;
	char** genes;
	unsigned short *rankData;
	double *qqNormedData;
	double *qqNormedZero;
	long *zeroPosPerCell;
	unordered_map<int, string> numToCell;
	unordered_map<string, int> cellToNum;
	int gene_count;
	int cell_count, data_count;
	int str_barcodes_len;
	int str_genes_length;

public:
	SparseMatrix() {
	}
	~SparseMatrix()
	{
		
	}
	void set_rank(unsigned short *rank);
	char** get_barcodes();
	char** get_genes();
	long long* get_indices();
	int* get_data();
	int get_cell_count();
	int get_gene_count();
	int get_data_count();
	long long* get_indptr();
	int get_str_barcodes_len();
	int get_str_genes_length();
	unordered_map<int, string> get_numToCell();
	int readHDF5File(string path, string type);
	void createCellnameMap();
	unsigned short* createCellVectorByName(string cellname);
	unordered_map<int,int> cellFiltration();
	void createZeroPosPerCell();
	void qqNormedData2HDF5Format();
	void write2HDF5(string path);
	void deleteSparseMatrix(string type);
	double** fetch_batch(int batch_index,int batch_size=128);
};
