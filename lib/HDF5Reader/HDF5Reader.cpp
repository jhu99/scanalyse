#include "HDF5Reader.h"
#define _CRT_SECURE_NO_WARNINGS
using namespace std;
using namespace H5;

int HDF5reader::readHDF5File(string path)
{
	hid_t       fileId, space, datasetId, type, groupId, attriId;
	herr_t      status;
	DataSet ds;
	H5T_class_t type_class;
	string paraCellName;
	DataType datatype;
	//open file
	const char*p = path.data();
	fileId = H5Fopen(p, H5F_ACC_RDONLY, H5P_DEFAULT);

	//open group
	string groupName = "/GRCh38";
	const char*group_name = groupName.data();
	status = H5Lget_info(fileId, group_name, NULL, H5P_DEFAULT);
	printf("/GRCh38': ");
	if (status == 0)
	{
		printf("The group exists.\n");
		groupId = H5Gopen2(fileId, group_name, H5P_DEFAULT);
		Group group = Group(groupId);
	}
	else
	{
		printf("The group either does NOT exist or some other error occurred.\n");
	}

	//open dataset-gene_names
	status = H5Lget_info(fileId, "/GRCh38/gene_names", NULL, H5P_DEFAULT);
	printf("/GRCh38/gene_names: ");
	if (status == 0)
	{
		printf("The dataset exists.\n");
		datasetId = H5Dopen(fileId, "/GRCh38/gene_names", H5P_DEFAULT);
		ds = DataSet(datasetId);
		datatype = ds.getDataType();
		size_t str_gene_names_len = datatype.getSize();
		size_t s = ds.getInMemDataSize();
		DataSpace dataspace = ds.getSpace();
		hsize_t genenum;
		dataspace.getSimpleExtentDims(&genenum, NULL);
		gene_count = genenum;
		gene_names = new char*[gene_count];
		for (int i = 0; i < gene_count; i++)
		{
			gene_names[i] = new char[str_gene_names_len];
		}
		char*para_genenames;
		para_genenames = new char[str_gene_names_len * gene_count];
		hid_t str_genenames = H5Tcopy(H5T_C_S1);
		H5Tset_size(str_genenames, str_gene_names_len);
		status = H5Dread(datasetId, str_genenames, H5S_ALL, H5S_ALL, H5P_DEFAULT, para_genenames);
		for (int i = 0; i < gene_count; i++)
		{
			strncpy(gene_names[i], para_genenames + i * str_gene_names_len, str_gene_names_len);
		}
		for (int i = 0; i < 10; i++)
		{
			cout << gene_names[i] << endl;
		}
		cout << "finish read gene_names" << endl;
	}
	else
	{
		printf("The dataset either does NOT exist or some other error occurred.\n");
	}
	//read dataset-genes
	status = H5Lget_info(fileId, "/GRCh38/genes", NULL, H5P_DEFAULT);
	printf("/GRCh38/genes: ");
	if (status == 0)
	{
		printf("The dataset exists.\n");
		datasetId = H5Dopen(fileId, "/GRCh38/genes", H5P_DEFAULT);
		ds = DataSet(datasetId);
		datatype = ds.getDataType();
		size_t str_genes_len = datatype.getSize();
		size_t s = ds.getInMemDataSize();

		genes = new char*[gene_count];
		for (int i = 0; i < gene_count; i++)
		{
			genes[i] = new char[str_genes_len];
		}
		char*para_genes;
		para_genes = new char[str_genes_len * gene_count];
		hid_t str_genes = H5Tcopy(H5T_C_S1);
		H5Tset_size(str_genes, str_genes_len);
		status = H5Dread(datasetId, str_genes, H5S_ALL, H5S_ALL, H5P_DEFAULT, para_genes);
		for (int i = 0; i < gene_count; i++)
		{
			strncpy(genes[i], para_genes + i * str_genes_len, str_genes_len);
		}
		for (int i = 0; i < 10; i++)
		{
			cout << genes[i] << endl;
		}
		cout << "finish read genes" << endl;
	}
	else
	{
		printf("The dataset either does NOT exist or some other error occurred.\n");
	}
	//read dataset-cellname
	status = H5Lget_info(fileId, "/GRCh38/barcodes", NULL, H5P_DEFAULT);
	printf("/GRCh38/barcodes': ");
	if (status == 0)
	{
		printf("The dataset exists.\n");
		datasetId = H5Dopen(fileId, "/GRCh38/barcodes", H5P_DEFAULT);
		ds = DataSet(datasetId);
		datatype = ds.getDataType();
		size_t str_cell_names_len = datatype.getSize();
		size_t s = ds.getInMemDataSize();
		DataSpace dataspace = ds.getSpace();
		hsize_t cellnum;
		dataspace.getSimpleExtentDims(&cellnum, NULL);
		cell_count = cellnum;
		char* para_barcodes;
		para_barcodes = new char[str_cell_names_len * cell_count];
		barcodes = new char*[cell_count];
		for (int i = 0; i < cell_count; i++)
		{
			barcodes[i] = new char[str_cell_names_len];
		}
		hid_t str_cell_names = H5Tcopy(H5T_C_S1);
		H5Tset_size(str_cell_names, str_cell_names_len);
		status = H5Dread(datasetId, str_cell_names, H5S_ALL, H5S_ALL, H5P_DEFAULT, para_barcodes);
		for (int i = 0; i < cell_count; i++)
		{
			strncpy(barcodes[i], para_barcodes + i * str_cell_names_len, str_cell_names_len);
		}
		cout << "finish read barcodes" << endl;
	}
	else
	{
		printf("The dataset either does NOT exist or some other error occurred.\n");
	}

	//read data
	status = H5Lget_info(fileId, "/GRCh38/data", NULL, H5P_DEFAULT);
	printf("/GRCh38/data ");
	if (status == 0)
	{
		printf("The dataset-data exists.\n");
		datasetId = H5Dopen(fileId, "/GRCh38/data", H5P_DEFAULT);
		ds = DataSet(datasetId);
		datatype = ds.getDataType();

		size_t s = ds.getInMemDataSize();
		DataSpace dataspace = ds.getSpace();
		hsize_t datanum;
		dataspace.getSimpleExtentDims(&datanum, NULL);
		data_count = datanum;
		data = new int[data_count];
		status = H5Dread(datasetId, datatype.getId(), H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
		cout << "finish read data" << endl;
	}
	else
	{
		printf("The dataset either does NOT exist or some other error occurred.\n");
	}

	//read indices
	status = H5Lget_info(fileId, "/GRCh38/indices", NULL, H5P_DEFAULT);
	printf("/GRCh38/indices: ");
	if (status == 0)
	{
		indices = new long long[data_count];
		datasetId = H5Dopen(fileId, "/GRCh38/indices", H5P_DEFAULT);
		ds = DataSet(datasetId);
		datatype = ds.getDataType();
		status = H5Dread(datasetId, datatype.getId(), H5S_ALL, H5S_ALL, H5P_DEFAULT, indices);
		cout << "finish read indices" << endl;
	}
	else
	{
		printf("The dataset either does NOT exist or some other error occurred.\n");
	}
	//read indptr
	status = H5Lget_info(fileId, "/GRCh38/indptr", NULL, H5P_DEFAULT);
	printf("/GRCh38/indptr': ");
	if (status == 0)
	{
		indptr = new long long[cell_count + 1];
		datasetId = H5Dopen(fileId, "/GRCh38/indptr", H5P_DEFAULT);
		ds = DataSet(datasetId);
		datatype = ds.getDataType();
		status = H5Dread(datasetId, datatype.getId(), H5S_ALL, H5S_ALL, H5P_DEFAULT, indptr);
		cout << "finish read indptr" << endl;
	}
	else
	{
		printf("The dataset either does NOT exist or some other error occurred.\n");
	}
	//close
	H5Fclose(fileId);
	H5Fclose(datasetId);
	H5Fclose(groupId);
	H5Fclose(status);
	return 0;
}



void HDF5reader::createCellnameMap()
{
	for (int i = 0; i < cell_count; i++)
	{
		numToCell[i] = barcodes[i];
		cellToNum[barcodes[i]] = i;
	}
}

int* HDF5reader::createCellVectorByName(string cellname)
{
	int cellPos = cellToNum[cellname];
	int *singleCellVector;
	singleCellVector = new int[gene_count];
	int columPos;
	int paraStart = 0;
	for (int i = 0; i < gene_count; i++)
	{
		singleCellVector[i] = 0;
	}
	for (long long paraPos = indptr[cellPos]; paraPos < indptr[cellPos + 1]; paraPos++)
	{
		columPos = indices[paraPos];
		singleCellVector[columPos] = data[paraPos];
	}
	return singleCellVector;
}

unordered_map<int, int> HDF5reader::cellFiltration()
{
	unordered_map<int, int> cellIdToNum;
	long long *genecount;
	genecount = new long long[cell_count];
	int min_count = gene_count / 100;
	for (int i = 0; i < cell_count; i++)
	{
		genecount[i] = indptr[i + 1] - indptr[i];
		if (genecount[i] <min_count)
		{
			cellIdToNum[i] = 1;
		}
	}
	return cellIdToNum;
}

char** HDF5reader::get_barcodes()
{
	return barcodes;
}

char ** HDF5reader::get_gene_names()
{
	return gene_names;
}
char ** HDF5reader::get_genes()
{
	return genes;
}
long long* HDF5reader::get_indices() {
	return indices;
}
int* HDF5reader::get_data() {
	return data;
}
int HDF5reader::get_cell_count()
{
	return cell_count;
}

int HDF5reader::get_gene_count()
{
	return gene_count;
}

int HDF5reader::get_data_count()
{
	return data_count;
}

long long * HDF5reader::get_indptr()
{
	return indptr;
}

void HDF5reader::deleteHDF(){
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
