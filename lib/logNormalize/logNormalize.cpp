#include "logNormalize.h"

logNormalize::logNormalize(SparseMatrix sm)
{
	this->sm = sm;
	this->gene_count = sm.get_gene_count();
	this->cell_count = sm.get_cell_count();
	this->indptr = sm.get_indptr();
	this->indices = sm.get_indices();
	this->data = sm.get_data();
	this->data_count = sm.get_data_count();
}

logNormalize::~logNormalize()
{
}

double * logNormalize::get_log_data()
{
	return normalize_data;
}

double * logNormalize::get_normalize_data()
{
	return normalize_data;
}

void logNormalize::cacuSizeFactor(int thread_index)
{
	int start = cell_count / thread_count * thread_index;
	int end;
	if (thread_index == thread_count-1)
	{
		end = cell_count ;
	}
	else
	{
		end = cell_count / thread_count * (thread_index + 1);
	}
	for (int cell_index = start; cell_index < end; cell_index++)
	{
		for (int para_pos = indptr[cell_index]; para_pos < indptr[cell_index + 1]; para_pos++)
		{
			size_factor[cell_index] += data[para_pos];
		}
	}
}

void logNormalize::logData(int thread_index)
{
	for (int cell_index = 0; cell_index < cell_count; cell_index++)
	{
		for (int para_pos = indptr[cell_index]; para_pos < indptr[cell_index + 1]; para_pos++)
		{
			normalize_data[para_pos] = log10(data[para_pos] * 10 ^ 6 / size_factor[cell_index] + 1);
		}
	}
}

void logNormalize::scaleData(int thread_index)
{
	double standard_deviation=gsl_stats_sd(normalize_data, 1, data_count);
	double mean = gsl_stats_mean(normalize_data, 1, data_count);
	for(int i=0;i<data_count;i++)
	{
		normalize_data[i] = (normalize_data[i] - mean) / standard_deviation;
	}
	
}

void logNormalize::logNormalizeData(int thread_count)
{
	normalize_data = new double[data_count];
	size_factor = new int[cell_count];
	fill(size_factor, size_factor + cell_count, 0);
	this->thread_count = thread_count;
	thread *threads;
	threads = new thread[thread_count];
	for (int i = 0; i < thread_count; i++)
	{
		threads[i]=thread(&logNormalize::cacuSizeFactor,this, i);
		threads[i].join();
	}
	for (int i = 0; i < thread_count; i++)
	{
		threads[i] = thread(&logNormalize::logData, this, i);
		threads[i].join();
	}
	for (int i = 0; i < thread_count; i++)
	{
		threads[i] = thread(&logNormalize::scaleData, this, i);
		threads[i].join();
	}
	cout <<"after normalize mean:" << gsl_stats_mean(normalize_data, 1, data_count) << endl;
	cout << "after normalize standard deviation:" << gsl_stats_sd(normalize_data, 1, data_count) << endl;
	cout << "end" << endl;
	delete []threads;
}
