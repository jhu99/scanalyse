#include "SparseMatrix.h"
#define _CRT_SECURE_NO_WARNINGS
using namespace std;
using namespace H5;

void SparseMatrix::set_rank(unsigned short *rank)
{
	rankData = rank;
}

char** SparseMatrix::get_barcodes()
{
	return barcodes;
}


char ** SparseMatrix::get_genes()
{
	return genes;
}

long long* SparseMatrix::get_indices() {
	return indices;
}

int* SparseMatrix::get_data() {
	return data;
}

int SparseMatrix::get_cell_count()
{
	return cell_count;
}

int SparseMatrix::get_gene_count()
{
	return gene_count;
}

int SparseMatrix::get_data_count()
{
	return data_count;
}

long long* SparseMatrix::get_indptr()
{
	return indptr;
}

int SparseMatrix::get_str_barcodes_len()
{
	return str_barcodes_len;
}

int SparseMatrix::get_str_genes_length()
{
	return str_genes_length;
}


unordered_map<int, string> SparseMatrix::get_numToCell()
{
	return numToCell;
}

int SparseMatrix::readHDF5File(string path, string type)
{
	hid_t       fileId, datasetId, groupId;
	herr_t      status;
	DataSet ds;
	H5T_class_t type_class;
	string paraCellName;
	DataType datatype;
	DataSpace dataspace;
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

	//read dataset-genes
	status = H5Lget_info(fileId, "/GRCh38/genes", NULL, H5P_DEFAULT);
	printf("/GRCh38/genes: ");
	if (status == 0)
	{
		printf("The dataset exists.\n");
		datasetId = H5Dopen(fileId, "/GRCh38/genes", H5P_DEFAULT);
		ds = DataSet(datasetId);
		dataspace = ds.getSpace();
		hsize_t genenum;
		dataspace.getSimpleExtentDims(&genenum, NULL);
		gene_count = genenum;

		datatype = ds.getDataType();
		size_t str_genes_len = datatype.getSize();
		str_genes_length = str_genes_len;
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
		str_barcodes_len = str_cell_names_len;
		size_t s = ds.getInMemDataSize();
		dataspace = ds.getSpace();
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
		dataspace = ds.getSpace();
		hsize_t datanum;
		dataspace.getSimpleExtentDims(&datanum, NULL);
		data_count = datanum;
		if (type == "original")
		{
			data = new int[data_count];
			status = H5Dread(datasetId, datatype.getId(), H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
			cout << "finish read data" << endl;
		}
		else
		{
			qqNormedData = new double[data_count];
			status = H5Dread(datasetId, datatype.getId(), H5S_ALL, H5S_ALL, H5P_DEFAULT, qqNormedData);
			cout << "finish read qqnorm_data" << endl;
		}
		
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
	//if it is qqnorm_data read dataset zero_value

	status = H5Lget_info(fileId, "/GRCh38/zero_value", NULL, H5P_DEFAULT);
	printf("/GRCh38/zero_value': ");
	if (status == 0)
	{
		qqNormedZero = new double[cell_count];
		datasetId = H5Dopen(fileId, "/GRCh38/zero_value", H5P_DEFAULT);
		ds = DataSet(datasetId);
		datatype = ds.getDataType();
		status = H5Dread(datasetId, datatype.getId(), H5S_ALL, H5S_ALL, H5P_DEFAULT, qqNormedZero);
		cout << "finish read zero_value" << endl;
	}
	else
	{
		printf("The dataset zero_value does NOT exist or some other error occurred.\n");
	}
	//close
	H5Fclose(fileId);
	H5Dclose(datasetId);
	H5Gclose(groupId);
	ds.close();
	datatype.close();
	dataspace.close();
	return 0;
}

void SparseMatrix::createCellnameMap()
{
	for (int i = 0; i < cell_count; i++)
	{
		numToCell[i] = barcodes[i];
		cellToNum[barcodes[i]] = i;
	}
}

unsigned short* SparseMatrix::createCellVectorByName(string cellname)
{
	int cellPos = cellToNum[cellname];
	unsigned short *singleCellVector;
	singleCellVector = new unsigned short[gene_count];
	int columPos;
	int paraStart = 0;
	int zeroCount = indptr[cellPos + 1]- indptr[cellPos];

	for (int i = 0; i < gene_count; i++)
	{
		singleCellVector[i] = (unsigned short)zeroCount/2;
	}
	for (long paraPos = indptr[cellPos]; paraPos < indptr[cellPos + 1]; paraPos++)
	{
		columPos = indices[paraPos];
		singleCellVector[columPos] = rankData[paraPos];
	}
	return singleCellVector;
}

unordered_map<int, int> SparseMatrix::cellFiltration()
{
	unordered_map<int, int> cellIdToNum;
	long *genecount;
	genecount = new long[cell_count];
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

void SparseMatrix::createZeroPosPerCell()
{
	cout << "start find zero_gene position for each cell." << endl;
	int startPos; 
	int randPos;
	zeroPosPerCell = new long[cell_count];
	for (int i = 0; i < cell_count; i++)
	{
		while (true)
		{
			randPos = rand() % gene_count;
			if (find(indices + indptr[i], indices + indptr[i + 1], randPos)== indices + indptr[i + 1])
			{
				zeroPosPerCell[i] = randPos;
				break;
			}
		}
	}
	
	cout << "end find zero_gene position for each cell." << endl;
}

void SparseMatrix::qqNormedData2HDF5Format()
{
	qqNorm qq;
	unsigned short *test;
	double *result;
	int paraPos;
	qqNormedData = new double[data_count];
	qqNormedZero = new double[cell_count];
	createZeroPosPerCell();
	cout << "start use rank to qqnorm." << endl;
	for (int i = 0; i < cell_count; i++)
	{
		test = createCellVectorByName(numToCell[i]);
		result = qq.caculateTheoryQuantiles(test, gene_count);
		for (paraPos = indptr[i]; paraPos < indptr[i + 1]; paraPos++)
		{
			qqNormedData[paraPos] = result[indices[paraPos]];
		}
		qqNormedZero[i] = result[zeroPosPerCell[i]];
	}
	cout << "end use rank to qqnorm." << endl;
}

void SparseMatrix::write2HDF5(string path)
{

	const H5std_string FILE_NAME(path);
	const H5std_string GROUP_NAME("GRCh38");
	const int RANK = 1;
	H5File file(FILE_NAME, H5F_ACC_TRUNC);
	Group group = file.createGroup(GROUP_NAME);

	//write Dataset
	hsize_t dims[RANK], chunk_size[1];
	//write barcodes
	chunk_size[0] = 3276;
	char *para_barcodes;
	int cell_name_str_len = 40;
	para_barcodes = new char[cell_count * cell_name_str_len];
	for (int i = 0; i < cell_count; i++)
	{
		strncpy(para_barcodes + i * cell_name_str_len, barcodes[i], cell_name_str_len);
	}
	dims[0] = cell_count;
	DataSpace *dataspace_barcodes = new DataSpace(RANK, dims);
	size_t str_cell_names_len = cell_name_str_len;
	hid_t barcodes_type = H5Tcopy(H5T_C_S1);
	H5Tset_size(barcodes_type, str_cell_names_len);
	DataSet *dataset_barcodes = new DataSet(group.createDataSet("barcodes", barcodes_type, *dataspace_barcodes));
	dataset_barcodes->write(para_barcodes, barcodes_type);
	cout << "finish write barcodes" << endl;
	//write data
	dims[0] = data_count;
	DataSpace *dataspace_data = new DataSpace(RANK, dims);
	DataSet *dataset_data = new DataSet(group.createDataSet("data", H5T_NATIVE_DOUBLE, *dataspace_data));
	H5Dwrite(dataset_data->getId(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, qqNormedData);
	cout << "finish write data" << endl;
	//write genes
	char *para_genes;
	int genes_str_len = 15;
	para_genes = new char[gene_count * genes_str_len];
	for (int i = 0; i < gene_count; i++)
	{
		strncpy(para_genes + i * genes_str_len, genes[i], genes_str_len);
	}
	dims[0] = gene_count;
	DataSpace *dataspace_genes = new DataSpace(RANK, dims);
	size_t str_genes_len = genes_str_len;
	hid_t genes_type = H5Tcopy(H5T_C_S1);
	H5Tset_size(genes_type, str_genes_len);
	DataSet *dataset_genes = new DataSet(group.createDataSet("genes", genes_type, *dataspace_genes));
	dataset_genes->write(para_genes, genes_type);
	cout << "finish write genes" << endl;
	//write indices
	dims[0] = data_count;
	DataSpace *dataspace_indices = new DataSpace(RANK, dims);
	DataSet *dataset_indices = new DataSet(group.createDataSet("indices", H5T_NATIVE_LONG, *dataspace_indices));
	H5Dwrite(dataset_indices->getId(), H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, indices);
	cout << "finish write indices" << endl;
	//write indptr
	dims[0] = cell_count+1;
	DataSpace *dataspace_indptr = new DataSpace(RANK, dims);
	DataSet *dataset_indptr = new DataSet(group.createDataSet("indptr", H5T_NATIVE_LONG, *dataspace_indptr));
	H5Dwrite(dataset_indptr->getId(), H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, indptr);
	cout << "finish write indptr" << endl;
	//write qqNormedZero
	dims[0] = cell_count;
	DataSpace *dataspace_zero_value = new DataSpace(RANK, dims);
	DataSet *dataset_zero_value = new DataSet(group.createDataSet("zero_value", H5T_NATIVE_DOUBLE, *dataspace_zero_value));
	H5Dwrite(dataset_zero_value->getId(), H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, qqNormedZero);
	cout << "finish write zero_value" << endl;

	delete dataspace_data;
	delete dataspace_barcodes;
	delete dataspace_genes;
	delete dataspace_indices;
	delete dataspace_indptr;
	delete dataspace_zero_value;
	delete dataset_barcodes;
	delete dataset_data;
	delete dataset_genes;
	delete dataset_indices;
	delete dataset_indptr;
	delete dataset_zero_value;
	group.close();
	file.close();
}

void SparseMatrix::deleteSparseMatrix(string type){
	delete[] indices;
	delete[] indptr;
	if (type == "orignal")
	{
		delete[] data;
	}
	else if (type == "fetch_batch")
	{
		delete[] qqNormedData;
		delete[] qqNormedZero;
	}
	else if(type=="write_qqnorm_data")
	{
		delete[] data;
		delete[] qqNormedData;
		delete[] rankData;
		delete[] zeroPosPerCell;
		delete[] qqNormedZero;
		delete[] rankData;
	}
	for (int i = 0; i < cell_count; i++)
	{
		delete[] barcodes[i];
	}
	delete[] barcodes;
	for (int i = 0; i < gene_count; i++)
	{
		delete[] genes[i];
	}
	delete[] genes;
	cout << "finish delete" << endl;
}

double ** SparseMatrix::fetch_batch(int batch_index, int batch_size)
{
	double **inputMatrix;
	long columnPos;
	inputMatrix = new double*[batch_size];
	for (int i = 0; i < batch_size; i++)
	{
		inputMatrix[i] = new double[gene_count];
	}
	int startPos = (batch_index - 1)*batch_size;
	for (int i = 0; i < batch_size; i++)
	{
		for (int j = 0; j < gene_count; j++)
		{
			inputMatrix[i][j] = qqNormedZero[i + startPos];
		}
		for (long paraPos = indptr[startPos + i]; paraPos < indptr[startPos + i + 1]; paraPos++)
		{
			columnPos = indices[paraPos];
			inputMatrix[i][columnPos] = qqNormedData[paraPos];
		}
	}

	return inputMatrix;
}
