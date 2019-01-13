#include"Filtration.h"

Filtration::Filtration(SparseMatrix sm)
{
	this->sm = sm;
	this->cell_count = sm.get_cell_count();
	this->indptr = sm.get_indptr();
	this->indices = sm.get_indices();
	this->data_count = sm.get_data_count();
	this->gene_count = sm.get_gene_count();
	this->data = sm.get_data();
	this->genes = sm.get_genes();
	this->barcodes = sm.get_barcodes();
	this->str_barcodes_len = sm.get_str_barcodes_len();
	this->str_genes_len = sm.get_str_genes_length();
}

Filtration::~Filtration()
{
}

int * Filtration::get_filt_gene_data()
{
	return filt_gene_data;
}

char ** Filtration::get_filt_barcodes()
{
	return filt_barcodes;
}

char ** Filtration::get_filt_genes()
{
	return filt_genes;
}

long long * Filtration::get_filt_gene_indptr()
{
	return filt_gene_indptr;
}

long long * Filtration::get_filt_gene_indices()
{
	return filt_gene_indices;
}

void Filtration::cellFiltration(int min_count)
{
	filt_data_count = data_count;
	int* geneCount_preCell;
	geneCount_preCell = new int[cell_count];
	for (int i = 0; i < cell_count; i++)
	{
		geneCount_preCell[i] = indptr[i + 1] - indptr[i];
		if (geneCount_preCell[i] <min_count)
		{
			filtered_cell[i] = 1;
			filt_data_count -= geneCount_preCell[i];
		}
	}
	filt_cell_count = cell_count - filtered_cell.size();
}

void Filtration::dataFiltrationByCell()
{
	filt_data = new int[filt_data_count];
	filt_indices = new long long[filt_data_count];
	filt_indptr = new long long[filt_cell_count + 1];
	filt_barcodes = new char*[filt_cell_count];
	for (int i = 0; i < filt_cell_count; i++)
	{
		filt_barcodes[i] = new char[str_barcodes_len];
	}
	int f_cell_index = 0;
	int f_data_index = 0;
	int f_indptr_index = 0;
	filt_indptr[0] = 0;
	for (int cell_index = 0; cell_index < cell_count; cell_index++)
	{
		if (filtered_cell[cell_index] != 1)
		{
			filt_barcodes[f_cell_index] = barcodes[cell_index];
			f_cell_index++;
			for (long long paraPos = indptr[cell_index]; paraPos < indptr[cell_index + 1]; paraPos++)
			{
				filt_data[f_data_index] = data[paraPos];
				filt_indices[f_data_index] = indices[paraPos];
				f_data_index++;
			}
			filt_indptr[f_cell_index] = f_data_index;
		}
	}
	cout << f_cell_index << endl;
	cout << filt_cell_count << endl;
	cout << "end filt cell" << endl;
}

void Filtration::geneFiltration(int topNum, int method)
{
	filt_gene_count = topNum;
	string *top_genes_str;
	long long *top_genes_index;
	if (method == 2)
	{
		geneExpressionTop geneETop(sm, topNum);
		geneETop.geneSort();
		top_genes_str = geneETop.getTop();
		top_genes_index = geneETop.get_top_index();
	}
	else
	{
		geneVariationTop GVT(sm, topNum);
		GVT.geneSort();
		top_genes_str = GVT.get_Top();
		top_genes_index = GVT.get_top_index();
	}
	for (int i = 0; i < topNum; i++)
	{
		nonfiltToFilt[top_genes_index[i]] = i;
	}
	for (int i = 0; i < topNum;i++)
	{
		filtered_genes[top_genes_index[i]] = 1;
	}
	filt_genes = new char*[topNum];
	for (int i = 0; i < topNum; i++)
	{
		filt_genes[i] = new char[str_genes_len];
	}
	for (int i = 0; i < topNum; i++)
	{
		filt_genes[i] = genes[top_genes_index[i]];
	}

}

void Filtration::dataFiltrationByGenes()
{
	int paraPos;
	filt_gene_data_count=0;
	for (int cell_index = 0; cell_index < filt_cell_count; cell_index++)
	{
		for (paraPos = filt_indptr[cell_index]; paraPos < filt_indptr[cell_index + 1]; paraPos++)
		{
			if (filtered_genes[filt_indices[paraPos]] == 1)
			{
				filt_gene_data_count++;
			}
		}
	}
	filt_gene_data = new int[filt_gene_data_count];
	filt_gene_indices = new long long[filt_gene_data_count];
	filt_gene_indptr = new long long[filt_cell_count];
	int filt_gene_data_index = 0;
	filt_gene_indptr[0] = 0;
	for (int cell_index = 0; cell_index < filt_cell_count; cell_index++)
	{
		for (paraPos = filt_indptr[cell_index]; paraPos < filt_indptr[cell_index + 1]; paraPos++)
		{
			if (filtered_genes[filt_indices[paraPos]] == 1)
			{
				filt_gene_data[filt_gene_data_index] = filt_data[paraPos];
				filt_gene_indices[filt_gene_data_index] = nonfiltToFilt[filt_indices[paraPos]];
				filt_gene_data_index++;
			}
		}
		filt_gene_indptr[cell_index + 1] = filt_gene_data_index;
	}
	cout << "end filt gene" << endl;
}

void Filtration::filtGeneAndCell(int cell_min_count, int gene_top_num, int gene_filt_method)
{
	cellFiltration(cell_min_count);
	dataFiltrationByCell();
	geneFiltration(gene_top_num, gene_filt_method);
	dataFiltrationByGenes();
}

void Filtration::printFiltResult()
{
	for (int i = 0; i < filt_gene_count; i++)
	{
		cout << filt_genes[i] << " ";
	}
	cout << endl;
	/*for (int i = 0; i < filt_cell_count; i++)
	{
		cout << filt_barcodes[i] << " ";
	}
	cout << endl;
	for (int i = 0; i < filt_gene_data_count; i++)
	{
		cout << filt_gene_data[i] << " ";
	}
	cout << endl;
	for (int i = 0; i < filt_gene_data_count; i++)
	{
		cout << filt_gene_indices[i] << " ";
	}*/
	cout << filt_data_count << endl;
	cout << filt_indptr[filt_cell_count];
	cout << endl;
	cout << filt_gene_data_count << endl;
	cout << filt_gene_indptr[filt_cell_count ];
	cout << endl;
}



