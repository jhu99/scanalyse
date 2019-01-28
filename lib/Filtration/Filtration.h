#include "SparseMatrix/SparseMatrix.h"
#include "geneTop/geneExpressionTop.h"
#include "geneVariationTop/geneVariationTop.h"
#include "H5Cpp.h"

class Filtration
{
private:
	SparseMatrix sm;
	unordered_map<int, int> filtered_cell;
	unordered_map<long, int> filtered_genes;
	unordered_map<string, int> genesToNum;
	unordered_map<long long, long long>nonfiltToFilt;

	int cell_count;
	int filt_cell_count;
	int data_count;
	int filt_data_count;
	int filt_gene_data_count;
	int gene_count;
	int filt_gene_count;
	int str_barcodes_len;
	int str_genes_len;
	int str_gene_names_len;

	long long * indptr;
	long long * indices;
	int *data;
	char** genes;
	char** barcodes;
	char** gene_names;

	long long *filt_indptr;
	long long *filt_indices;
	int *filt_data;
	char **filt_barcodes;
	char **filt_genes;
	char **filt_gene_names;
	long long *filt_gene_indptr;
	long long *filt_gene_indices;
	int *filt_gene_data;


public:
	Filtration(SparseMatrix sm);
	~Filtration();
	int* get_filt_gene_data();
	char** get_filt_barcodes();
	char** get_filt_genes();
	long long* get_filt_gene_indptr();
	long long* get_filt_gene_indices();

	void cellFiltration(int min_count = 1);
	void dataFiltrationByCell();
	void geneFiltration(int topNum, int method = 1);
	void dataFiltrationByGenes();
	SparseMatrix filtGeneAndCell(int cell_min_count, int gene_top_num, int gene_filt_method);
	void printFiltResult();
	void writeFiltH5File(string write_path);
};