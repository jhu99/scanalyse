#pragma once
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include "H5Cpp.h"
#include <hdf5.h>
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include "qqNorm/qqNorm.h"
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;
class SparseMatrix
{
private:
	int *data;
	long long *indptr, *indices;
	char** barcodes;
	char** genes;
	char** gene_names;
	unsigned short *rankData;
	double *log_normalize_data;
	double *qqNormedData;
	double *qqNormedZero;
	long *zeroPosPerCell;
	unordered_map<int, string> numToCell;
	unordered_map<string, int> cellToNum;
	int gene_count;
	int cell_count;
	long long data_count;
	int str_barcodes_len;
	int str_genes_length;
	int str_gene_names_len;
	long long  indptrIndex;
	long long  indicesIndex;
	int geneIndex;
	long long cellIndex;
	int dataIndex;
	int geneNamesIndex;
public:
	SparseMatrix() {
	}
	SparseMatrix(char** barcodes, char** gene_names, char** genes, 
		long long* indptr, long long* indices, int* data, 
		int cell_count, int gene_count, int data_count, 
		int str_barcodes_len, int str_genes_len, int str_gene_names_len) {
		this->barcodes = barcodes;
		this->gene_names = gene_names;
		this->genes = genes;
		this->indptr = indptr;
		this->indices = indices;
		this->data = data;
		this->cell_count = cell_count;
		this->gene_count = gene_count;
		this->data_count = data_count;
		this->str_barcodes_len = str_barcodes_len;
		this->str_genes_length = str_genes_len;
		this->str_gene_names_len = str_gene_names_len;
	}
	~SparseMatrix()
	{

	}
	void set_rank(unsigned short *rank);
	void set_log_data(double *log_normalize_data);
	char** get_barcodes();
	char** get_genes();
	char** get_gene_names();
	long long* get_indices();
	int* get_data();
	int get_cell_count();
	int get_gene_count();
	long long get_data_count();
	long long* get_indptr();
	int get_str_barcodes_len();
	int get_str_genes_length();
	int get_str_gene_names_len();
	unordered_map<int, string> get_numToCell();
	int readHDF5File(string path, string type);
	void createCellnameMap();
	unsigned short* createCellVectorByName(string cellname);
	unordered_map<int, int> cellFiltration();
	void createZeroPosPerCell();
	void qqNormedData2HDF5Format();
	void write2HDF5(string path);
	void write2CSV(string path,string type);
	void deleteSparseMatrix(string type);
	double** fetch_batch(int batch_index, string norm_type,int batch_size = 128);
	void readMtxFile(string read_path);
	void readTsvFile(string read_path);
	void read_10x_h5(string read_path);
	void read_10x_mtx(string read_path);
	void mergeDate(vector<string> paths);
	void h5Compressed(string aimFilePath, string method, int chunk, int rank);
	void write_norm_data(string write_path, string norm_type, int chunk, string method);
};
