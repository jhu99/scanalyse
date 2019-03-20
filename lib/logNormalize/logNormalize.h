#include"SparseMatrix/SparseMatrix.h"
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <gsl/gsl_statistics.h>
#include <thread>
#include <math.h>
#include <mutex>

class logNormalize
{
	int *size_factor_per_cell;
	SparseMatrix sm;
	int *size_factor;
	int *data;
	long long *indptr;
	long long *indices;
	int cell_count;
	int gene_count;
	double *normalize_data;
	int data_count;
	int thread_count;
public:
	logNormalize(SparseMatrix sm);
	~logNormalize();
	double* get_log_data();
	double *get_normalize_data();
	int* get_size_factor();
	void cacuSizeFactor(int thread_index);
	void logData(int thread_index);
	void scaleData(int thread_index);
	void logNormalizeData(int thread_count);
};