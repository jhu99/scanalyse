#include"SparseMatrix/SparseMatrix.h"
class geneFilter{
	string filt_gene_file_path;
	SparseMatrix template_sm;
	char** genes;
	char** barcodes;
	char** gene_names;
	long long *indptr;
	long long *indices;
	int* data;
	int cell_count;
	int gene_count;
	int str_genes_len;
	int str_gene_names_len;
	int str_barcodes_len;

	char** filt_genes;
	char** filt_gene_names;
	int *filt_data;
	long long* filt_indptr;
	long long* filt_indices;
	int filt_gene_count;
	int filt_data_count;
	unordered_map<string, int> filt_genes_template;
	unordered_map<long long, long long> orig2filt_index;
public:
	geneFilter();
	~geneFilter();
	void createFiltGeneTemplate(string read_path);
	void out2CsvFile(string write_path);
	void createFiltGeneMap(string gene_path);
	void filtGene(string file_path);
	SparseMatrix filtDataByTemplateGene();
	void filt2H5File(string template_path, string read_path, string write_path);
};