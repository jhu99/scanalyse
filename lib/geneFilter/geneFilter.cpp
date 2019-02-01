#include "geneFilter.h"

geneFilter::geneFilter()
{
}

geneFilter::~geneFilter()
{
}

void geneFilter::createFiltGeneTemplate(string read_path)
{
	template_sm.read_10x_h5(read_path);
	gene_count = template_sm.get_gene_count();
	cell_count = template_sm.get_cell_count();
	long long* template_indptr = template_sm.get_indptr();
	long long* template_indices = template_sm.get_indices();
	char** template_genes = template_sm.get_genes();
	int *gene_express_count;
	gene_express_count = new int[gene_count];
	fill(gene_express_count, gene_express_count + gene_count, 0);
	int column_pos;
	for (int cell_index = 0; cell_index < cell_count; cell_index++)
	{
		for (int para_index = template_indptr[cell_index]; para_index < template_indptr[cell_index+1]; para_index++)
		{
			column_pos = template_indices[para_index];
			gene_express_count[column_pos]++;
		}
	}
	filt_gene_count = 0;
	for (int i = 0; i < gene_count; i++)
	{
		if (gene_express_count[i] != 0)
		{
			filt_gene_count++;
		}
	}
	filt_genes = new char*[filt_gene_count];
	int filt_index = 0;
	for (int i = 0; i < gene_count; i++)
	{
		if (gene_express_count[i] != 0)
		{
			filt_genes[filt_index] = template_genes[i];
			filt_index++;
		}
	}
	cout << "filt gene count:" << filt_gene_count << endl;
}

void geneFilter::out2CsvFile(string write_path)
{
	cout << "start out file" << endl;
	ofstream ofn(write_path);
	for (int i = 0; i < filt_gene_count; i++)
	{
		ofn << filt_genes[i] << "\n";
	}
	ofn.close();
	cout << "end out file" << endl;
}

void geneFilter::createFiltGeneMap(string gene_path)
{
	string para_genes;
	cout << "start create filt_genes_template map" << endl;
	ifstream inFile(gene_path, ios::in);
	while (getline(inFile, para_genes))
	{
		filt_genes_template[para_genes] = 1;
	}
	cout << "end create filt_genes_template map" << endl;
	filt_gene_count = filt_genes_template.size();
}

void geneFilter::filtGene(string file_path)
{
	SparseMatrix sm;
	sm.read_10x_mtx(file_path);

	gene_count = sm.get_gene_count();
	cell_count = sm.get_cell_count();
	indptr = sm.get_indptr();
	indices = sm.get_indices();
	data = sm.get_data();
	gene_names = sm.get_gene_names();
	genes = sm.get_genes();
	barcodes = sm.get_barcodes();
	str_genes_len = sm.get_str_genes_length();
	str_gene_names_len = sm.get_str_gene_names_len();
	str_barcodes_len = sm.get_str_barcodes_len();
	filt_genes = new char*[filt_gene_count];
	filt_gene_names = new char*[filt_gene_count];
	for (int i = 0; i < filt_gene_count; i++)
	{
		filt_genes[i] = new char[str_genes_len];
		filt_gene_names[i] = new char[str_gene_names_len];
	}
	int filt_gene_names_index = 0;
	for (int i = 0; i < gene_count; i++)
	{
		if (filt_genes_template[genes[i]] == 1)
		{
			filt_gene_names[filt_gene_names_index] = gene_names[i];
			filt_genes[filt_gene_names_index] = genes[i];
			orig2filt_index[i] = filt_gene_names_index;
			filt_gene_names_index++;
		}
	}
}

SparseMatrix geneFilter::filtDataByTemplateGene()
{
	int paraPos;
	filt_data_count = 0;
	for (int cell_index = 0; cell_index < cell_count; cell_index++)
	{
		for (paraPos = indptr[cell_index]; paraPos < indptr[cell_index + 1]; paraPos++)
		{
			if (filt_genes_template[genes[indices[paraPos]]] == 1)
			{
				filt_data_count++;
			}
		}
	}
	filt_data = new int[filt_data_count];
	filt_indices = new long long[filt_data_count];
	filt_indptr = new long long[cell_count + 1];
	int filt_data_index = 0;
	filt_indptr[0] = 0;
	for (int cell_index = 0; cell_index < cell_count; cell_index++)
	{
		for (paraPos = indptr[cell_index]; paraPos < indptr[cell_index + 1]; paraPos++)
		{
			if (filt_genes_template[genes[indices[paraPos]]]==1)
			{
				filt_data[filt_data_index] = data[paraPos];
				filt_indices[filt_data_index] = orig2filt_index[indices[paraPos]];
				filt_data_index++;
			}
		}
		filt_indptr[cell_index + 1] = filt_data_index;
	}
	SparseMatrix filt_sm(barcodes, filt_gene_names, filt_genes,
		filt_indptr, filt_indices, filt_data,
		cell_count, filt_gene_count, filt_data_count,
		str_barcodes_len, str_genes_len, str_gene_names_len);
	return filt_sm;
}

void geneFilter::filt2H5File(string template_path, string read_path, string write_path)
{
	cout << "start filt gene" << endl;
	createFiltGeneMap(template_path);
	filtGene(read_path);
	SparseMatrix result_sm = filtDataByTemplateGene();
	cout << "end filt gene" << endl;
	cout << "start write" << endl;
	result_sm.write_norm_data(write_path, "none", 500, "s");
	cout << "end write" << endl;
	result_sm.deleteSparseMatrix("original");
}

