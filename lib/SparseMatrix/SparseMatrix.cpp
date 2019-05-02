#include "SparseMatrix.h"
#define _CRT_SECURE_NO_WARNINGS
using namespace std;
using namespace H5;


void SparseMatrix::set_rank(unsigned short *rank)
{
	rankData = rank;
}

void SparseMatrix::set_log_data(double * log_normalize_data)
{
	this->log_normalize_data = log_normalize_data;
}

void SparseMatrix::set_size_factor(int * size_factor)
{
	this->size_factor = size_factor;
}

char** SparseMatrix::get_barcodes()
{
	return barcodes;
}


char ** SparseMatrix::get_genes()
{
	return genes;
}

char ** SparseMatrix::get_gene_names()
{
	return gene_names;
}

long long* SparseMatrix::get_indices() {
	return indices;
}

int* SparseMatrix::get_data() {
	return data;
}

int * SparseMatrix::get_rank_zero()
{
	return rank_zero;
}

int SparseMatrix::get_cell_count()
{
	return cell_count;
}

int SparseMatrix::get_gene_count()
{
	return gene_count;
}

long long SparseMatrix::get_data_count()
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

int SparseMatrix::get_str_gene_names_len()
{
	return str_gene_names_len;
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

	//read dataset-gene_names
	status = H5Lget_info(fileId, "/GRCh38/gene_names", NULL, H5P_DEFAULT);
	printf("/GRCh38/gene_names: ");
	if (status == 0)
	{
		printf("The dataset exists.\n");
		datasetId = H5Dopen(fileId, "/GRCh38/gene_names", H5P_DEFAULT);
		ds = DataSet(datasetId);
		dataspace = ds.getSpace();
		datatype = ds.getDataType();
		size_t gene_names_len = datatype.getSize();
		str_gene_names_len = gene_names_len;

		gene_names = new char*[gene_count];
		for (int i = 0; i < gene_count; i++)
		{
			gene_names[i] = new char[str_gene_names_len];
		}
		char* para_gene_names;
		para_gene_names = new char[str_gene_names_len * gene_count];
		hid_t str_gene_names = H5Tcopy(H5T_C_S1);
		H5Tset_size(str_gene_names, str_gene_names_len);
		status = H5Dread(datasetId, str_gene_names, H5S_ALL, H5S_ALL, H5P_DEFAULT, para_gene_names);
		for (int i = 0; i < gene_count; i++)
		{
			strncpy(gene_names[i], para_gene_names + i * str_gene_names_len, str_gene_names_len);
		}
		for (int i = 0; i < 5; i++)
		{
			cout << gene_names[i] << endl;
		}
		cout << "finish read gene_names" << endl;
	}
	else
	{
		printf("The dataset either does NOT exist or some other error occurred.\n");
	}

	//read dataset-barcodes
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
		cout << datanum << endl;
		if (type == "original")
		{
			data = new int[data_count];
			status = H5Dread(datasetId, datatype.getId(), H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
			cout << "finish read data" << endl;
		}
		else if (type == "qqNorm")
		{
			qqNormedData = new double[data_count];
			status = H5Dread(datasetId, datatype.getId(), H5S_ALL, H5S_ALL, H5P_DEFAULT, qqNormedData);
			cout << "finish read qqnorm_data" << endl;
		}
		else if (type == "rank")
		{
			rankData = new unsigned short[data_count];
			status = H5Dread(datasetId, datatype.getId(), H5S_ALL, H5S_ALL, H5P_DEFAULT, rankData);
			cout << "finish read rank_data" << endl;
		}
		else if (type == "log")
		{
			log_normalize_data = new double[data_count];
			status = H5Dread(datasetId, datatype.getId(), H5S_ALL, H5S_ALL, H5P_DEFAULT, log_normalize_data);
			cout << "finish read log_data" << endl;
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
	if (type == "qqNorm")
	{
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
	}
	if (type == "rank")
	{
		status = H5Lget_info(fileId, "/GRCh38/zero_value", NULL, H5P_DEFAULT);
		printf("/GRCh38/zero_value': ");
		if (status == 0)
		{
			rank_zero = new int[cell_count];
			datasetId = H5Dopen(fileId, "/GRCh38/zero_value", H5P_DEFAULT);
			ds = DataSet(datasetId);
			datatype = ds.getDataType();
			status = H5Dread(datasetId, datatype.getId(), H5S_ALL, H5S_ALL, H5P_DEFAULT, rank_zero);
			cout << "finish read zero_value" << endl;
		}
		else
		{
			printf("The dataset zero_value does NOT exist or some other error occurred.\n");
		}
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

void SparseMatrix::cacuRankZero()
{
	rank_zero = new int[cell_count];
	for (int i = 0; i < cell_count; i++)
	{
		rank_zero[i] = (gene_count - (indptr[i + 1] - indptr[i])) / 2;
	}
}

unsigned short* SparseMatrix::createCellVectorByName(string cellname)
{
	int cellPos = cellToNum[cellname];
	unsigned short *singleCellVector;
	singleCellVector = new unsigned short[gene_count];
	int columPos;
	int paraStart = 0;
	int zeroCount = indptr[cellPos + 1] - indptr[cellPos];

	for (int i = 0; i < gene_count; i++)
	{
		singleCellVector[i] = (unsigned short)zeroCount / 2;
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
		if (genecount[i] < min_count)
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
			if (find(indices + indptr[i], indices + indptr[i + 1], randPos) == indices + indptr[i + 1])
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
	hsize_t dims[RANK];
	//write barcodes
	char *para_barcodes;
	para_barcodes = new char[cell_count * str_barcodes_len];
	for (int i = 0; i < cell_count; i++)
	{
		strncpy(para_barcodes + i * str_barcodes_len, barcodes[i], str_barcodes_len);
	}
	dims[0] = cell_count;
	DataSpace *dataspace_barcodes = new DataSpace(RANK, dims);
	size_t str_cell_names_len = str_barcodes_len;
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
	para_genes = new char[gene_count * str_genes_length];
	for (int i = 0; i < gene_count; i++)
	{
		strncpy(para_genes + i * str_genes_length, genes[i], str_genes_length);
	}
	dims[0] = gene_count;
	DataSpace *dataspace_genes = new DataSpace(RANK, dims);
	size_t str_genes_len = str_genes_length;
	hid_t genes_type = H5Tcopy(H5T_C_S1);
	H5Tset_size(genes_type, str_genes_len);
	DataSet *dataset_genes = new DataSet(group.createDataSet("genes", genes_type, *dataspace_genes));
	dataset_genes->write(para_genes, genes_type);
	cout << "finish write genes" << endl;
	//write gene_names
	char *para_gene_names;
	para_gene_names = new char[gene_count * str_gene_names_len];
	for (int i = 0; i < gene_count; i++)
	{
		strncpy(para_genes + i * str_gene_names_len, gene_names[i], str_gene_names_len);
	}
	dims[0] = gene_count;
	DataSpace *dataspace_gene_names = new DataSpace(RANK, dims);
	size_t str_gene_names_len_t = str_gene_names_len;
	hid_t gene_names_type = H5Tcopy(H5T_C_S1);
	H5Tset_size(gene_names_type, str_gene_names_len_t);
	DataSet *dataset_gene_names = new DataSet(group.createDataSet("gene_names", gene_names_type, *dataspace_gene_names));
	dataset_gene_names->write(para_gene_names, gene_names_type);
	cout << "finish write gene_names" << endl;
	//write indices
	dims[0] = data_count;
	DataSpace *dataspace_indices = new DataSpace(RANK, dims);
	DataSet *dataset_indices = new DataSet(group.createDataSet("indices", H5T_NATIVE_LONG, *dataspace_indices));
	H5Dwrite(dataset_indices->getId(), H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, indices);
	cout << "finish write indices" << endl;
	//write indptr
	dims[0] = cell_count + 1;
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
	delete dataspace_gene_names;
	delete dataset_barcodes;
	delete dataset_data;
	delete dataset_genes;
	delete dataset_indices;
	delete dataset_indptr;
	delete dataset_zero_value;
	delete dataset_gene_names;
	group.close();
	file.close();
}

void SparseMatrix::write2CSV(string path, string type) {
	cout << "start writing" << endl;
	long long start, end;
	ofstream ofn(path);
	ofn << "";
	for (int i = 0; i < gene_count; i++) {
		ofn << "," << genes[i];
	}

	if (type == "original") {
		int * cell;
		cell = new int[gene_count];
		for (int i = 0; i < cell_count; i++) {
			ofn << "\n";
			ofn << barcodes[i];
			fill(cell, cell + gene_count, 0);
			start = indptr[i];
			end = indptr[i + 1];
			for (long long j = start; j < end; j++) {
				cell[indices[j]] = data[j];
			}
			for (int j = 0; j < gene_count; j++) {
				ofn << "," << cell[j];
			}
			if (i % 100 == 0)
				cout << "write " << i << " line" << endl;
		}
		delete cell;
	}
	else if (type == "qqNorm") {
		double * cell;
		cell = new double[gene_count];
		for (int i = 0; i < cell_count; i++) {
			ofn << "\n";
			ofn << barcodes[i];
			fill(cell, cell + gene_count, 0);
			start = indptr[i];
			end = indptr[i + 1];
			for (long long j = start; j < end; j++) {
				cell[indices[j]] = qqNormedData[j];
			}
			for (int j = 0; j < gene_count; j++) {
				ofn << "," << cell[j];
			}
			if (i % 100 == 0)
				cout << "write " << i << " line" << endl;
		}
		delete cell;
	}
	else if (type == "rank") {
		double * cell;
		cell = new double[gene_count];
		for (int i = 0; i < cell_count; i++) {
			ofn << "\n";
			ofn << barcodes[i];
			fill(cell, cell + gene_count, rank_zero[i]);
			start = indptr[i];
			end = indptr[i + 1];
			for (long long j = start; j < end; j++) {
				cell[indices[j]] = rankData[j];
			}
			for (int j = 0; j < gene_count; j++) {
				ofn << "," << cell[j];
			}
			if (i % 100 == 0)
				cout << "write " << i << " line" << endl;
		}
		delete cell;
	}
	else if (type == "log") {
		double * cell;
		cell = new double[gene_count];
		for (int i = 0; i < cell_count; i++) {
			ofn << "\n";
			ofn << barcodes[i];
			fill(cell, cell + gene_count, 0);
			start = indptr[i];
			end = indptr[i + 1];
			for (long long j = start; j < end; j++) {
				cell[indices[j]] = log_normalize_data[j];
			}
			for (int j = 0; j < gene_count; j++) {
				ofn << "," << cell[j];
			}
			if (i % 100 == 0)
				cout << "write " << i << " line" << endl;
		}
		delete cell;
	}

	ofn.close();
	cout << "end writing" << endl;
}

void SparseMatrix::deleteSparseMatrix(string type) {
	delete[] indices;
	delete[] indptr;
	if (type == "original")
	{
		delete[] data;
	}
	else if (type == "fetch_batch")
	{
		delete[] qqNormedData;
		delete[] qqNormedZero;
	}
	else if (type == "write_qqnorm_data")
	{
		delete[] data;
		delete[] qqNormedData;
		delete[] rankData;
		delete[] zeroPosPerCell;
		delete[] qqNormedZero;
	}
	else if (type == "log")
	{
		delete[] log_normalize_data;
	}
	else if (type == "rank")
	{
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
	for (int i = 0; i < gene_count; i++)
	{
		delete[] gene_names[i];
	}
	delete[] gene_names;
	cout << "finish delete" << endl;
}

double ** SparseMatrix::fetch_batch(int batch_index, string norm_type, int batch_size)
{
	double **input_matrix;
	long long column_pos;

	int start_pos = batch_index * batch_size;
	int read_size = batch_size;
	if (cell_count - start_pos < batch_size)
	{
		read_size = cell_count - start_pos;
	}
	input_matrix = new double*[read_size];
	for (int i = 0; i < read_size; i++)
	{
		input_matrix[i] = new double[gene_count];
		fill(input_matrix[i], input_matrix[i] + gene_count, 0);
	}
	if (norm_type == "qqNorm")
	{
		for (int i = 0; i < read_size; i++)
		{
			for (int j = 0; j < gene_count; j++)
			{
				input_matrix[i][j] = qqNormedZero[i + start_pos];
			}
			for (long long para_pos = indptr[start_pos + i]; para_pos < indptr[start_pos + i + 1]; para_pos++)
			{
				column_pos = indices[para_pos];
				input_matrix[i][column_pos] = qqNormedData[para_pos];
			}
		}
	}
	else if (norm_type == "rank")
	{
		for (int i = 0; i < read_size; i++)
		{
			for (long long para_pos = indptr[start_pos + i]; para_pos < indptr[start_pos + i + 1]; para_pos++)
			{
				column_pos = indices[para_pos];
				input_matrix[i][column_pos] = rankData[para_pos];
			}
		}
	}
	else if (norm_type == "log")
	{
		for (int i = 0; i < read_size; i++)
		{
			for (long long para_pos = indptr[start_pos + i]; para_pos < indptr[start_pos + i + 1]; para_pos++)
			{
				column_pos = indices[para_pos];
				input_matrix[i][column_pos] = log_normalize_data[para_pos];
			}
		}
	}
	return input_matrix;
}

void SparseMatrix::readMtxFile(string read_path)
{

	string mtx_path = read_path + "matrix.mtx";
	std::ifstream fin(mtx_path);
	// Ignore headers and comments:
	while (fin.peek() == '%') fin.ignore(2048, '\n');

	// Read defining parameters:
	fin >> gene_count >> cell_count >> data_count;
	cout << "gene_count:" << gene_count << endl;
	cout << "cell_count:" << cell_count << endl;
	cout << "data_count:" << data_count << endl;
	data = new int[data_count];
	indices = new long long[data_count];
	indptr = new long long[cell_count + 1];
	int *gene_count_per_cell;
	gene_count_per_cell = new int[cell_count];
	fill(gene_count_per_cell, gene_count_per_cell + cell_count, 0);

	// Read the data
	int cell_index;
	for (int data_index = 0; data_index < data_count; data_index++)
	{
		int para_data;
		long long para_indptr;
		long long para_indices;
		fin >> indices[data_index] >> cell_index >> data[data_index];
		gene_count_per_cell[cell_index - 1]++;
	}
	indptr[0] = 0;
	for (int cell_index = 0; cell_index < cell_count; cell_index++)
	{
		indptr[cell_index + 1] = indptr[cell_index] + gene_count_per_cell[cell_index];
	}
	cout << "indptr[cell_count] " << indptr[cell_count] << endl;
	fin.close();
}

void SparseMatrix::readTsvFile(string read_path)
{
	string barcodes_path = read_path + "barcodes.tsv";
	string genes_path = read_path + "genes.tsv";
	char separator = '\t';
	string lineStr;
	barcodes = new char*[cell_count];
	genes = new char*[gene_count];
	gene_names = new char*[gene_count];
	//read barcodes
	ifstream barcodes_infile(barcodes_path, ios::in);
	getline(barcodes_infile, lineStr);
	str_barcodes_len = lineStr.length();
	for (int i = 0; i < cell_count; i++)
	{
		barcodes[i] = new char[str_barcodes_len];
	}
	lineStr.copy(barcodes[0], str_barcodes_len, 0);
	int barcodes_index = 1;
	while (getline(barcodes_infile, lineStr))
	{
		lineStr.copy(barcodes[barcodes_index], str_barcodes_len, 0);
		barcodes_index++;
	}
	barcodes_infile.close();

	//read genes and gene_names
	ifstream genes_infile(genes_path, ios::in);
	str_genes_length = 15;
	str_gene_names_len = 19;
	string para_str;
	int is_genes = 1;
	for (int i = 0; i < gene_count; i++)
	{
		genes[i] = new char[str_genes_length];
		gene_names[i] = new char[str_gene_names_len];
	}
	int genes_index = 0;
	while (getline(genes_infile, lineStr))
	{
		stringstream ss(lineStr);
		getline(ss, para_str, separator);
		para_str.copy(genes[genes_index], str_genes_length, 0);
		getline(ss, para_str, '\n');
		para_str.copy(gene_names[genes_index], para_str.length(), 0);
		genes_index++;
	}
	genes_infile.close();
	for (int i = 0; i < 5; i++)
	{
		cout << genes[i] << endl;
		cout << gene_names[i] << endl;
	}
}

void SparseMatrix::read_10x_h5(string read_path)
{
	readHDF5File(read_path, "original");
}

void SparseMatrix::read_10x_mtx(string read_path)
{
	readMtxFile(read_path);
	readTsvFile(read_path);
}

void SparseMatrix::mergeDate(std::vector<std::string> paths,bool log,string log_path) {
	data_count = 0;
	cell_count = 0;
	gene_count = 0; 
	int path_len = 0;
	SparseMatrix *sm = new SparseMatrix[paths.size()];
	for (int i = 0; i < paths.size(); i++) {
		sm[i].readHDF5File(paths[i], "original");
		data_count += sm[i].get_data_count();
		cell_count += sm[i].get_cell_count();
		gene_count = sm[i].get_gene_count();
		str_genes_length = sm[i].get_str_genes_length();
		str_barcodes_len = sm[i].get_str_barcodes_len();
		path_len = max(path_len,(int)paths[i].size());
	}
	cout << "data_num" << data_count << endl;
	cout << "barcodes_num" << cell_count << endl;
	data = new int[data_count];
	genes = new char *[gene_count];
	for (int i = 0; i < gene_count; i++) {
		genes[i] = new char[str_genes_length];
	}
	gene_names = new char*[gene_count];
	for (int i = 0; i < gene_count; i++) {
		gene_names[i] = new char[str_gene_names_len];
	}

	barcodes = new char *[cell_count];
	for (int i = 0; i < cell_count; i++) {
		barcodes[i] = new char[str_barcodes_len];
	}
	indptr = new long long[cell_count + 1];
	indices = new long long[data_count];

	dataIndex = 0;
	cellIndex = 0;
	indptrIndex = 1;
	indicesIndex = 0;
	geneIndex = 0;
	geneNamesIndex = 0;
	indptr[0] = 0;
	for (int i = 0; i < paths.size(); i++) {
		for (long long j = 0; j < sm[i].get_data_count(); j++, dataIndex++) {
			data[dataIndex] = sm[i].get_data()[j];
			indices[dataIndex] = sm[i].get_indices()[j];
		}

		for (long long j = 1; j < sm[i].get_cell_count() + 1; j++, indptrIndex++) {
			indptr[indptrIndex] = indptr[indptrIndex - 1] + sm[i].get_indptr()[j] - sm[i].get_indptr()[j - 1];
		}
		for (long long j = 0; j < sm[i].get_cell_count(); j++, cellIndex++) {
			strcpy(barcodes[cellIndex], sm[i].get_barcodes()[j]);
		}

	}

	for (long long j = 0; j < sm[0].get_gene_count(); j++, geneIndex++) {
		strcpy(genes[geneIndex], sm[0].get_genes()[j]);
	}
	for (long long j = 0; j < sm[0].get_gene_count(); j++, geneNamesIndex++) {
		strcpy(gene_names[geneNamesIndex], sm[0].get_gene_names()[j]);
	}
	cout << "merge finish" << endl;

	if(log){
		char** cell_type;
		cell_type = new char *[cell_count];
		for (int i = 0; i < cell_count; i++) {
			cell_type[i] = new char[path_len];
		}
		int k = 0;
		string type;
		for(int i = 0;i < paths.size(); i++) {
			std::string::size_type nPos1 = std::string::npos;
			nPos1 = paths[i].find_last_of("/");
			if(nPos1 !=-1) {
    				type = paths[i].substr(nPos1+1, paths[i].size()-nPos1-4);
			}

			for(int j = 0;j < sm[i].get_cell_count();j++,k++){
				strcpy(cell_type[k], type.c_str());
			}		
		}
		ofstream ofn(log_path);
		bool flag = 0;
		for(long long i = 0;i < cellIndex;i++){
			if(flag){
				ofn<<"\n"<<barcodes[i]<<","<<cell_type[i];
			}
			else{
				ofn<<barcodes[i]<<","<<cell_type[i];
				flag = 1;
			}
			
		}
		ofn.close();
		cout<<"log finish"<<endl;
	}

}

void SparseMatrix::h5Compressed(string aimFilePath, string method, int chunk, int rank) {
	hid_t    file_id, dataset_id, group_id, dataspace_id;
	hid_t    plist_id;
	hid_t    strGene, strCell, strGeneName;
	size_t   nelmts;
	unsigned flags, filter_info;
	H5Z_filter_t filter_type;

	herr_t   status;
	hsize_t  dims[1];
	hsize_t  cdims[1];

	int      idx;
	int      i, j, numfilt;
	int shape[] = { gene_count,cell_count };

	string groupPath = "/GRCh38";
	string dataPath = "/GRCh38/data";
	string genesPath = "/GRCh38/genes";
	string geneNamesPath = "/GRCh38/gene_names";
	string indicesPath = "/GRCh38/indices";
	string indptrPath = "/GRCh38/indptr";
	string barcodesPath = "/GRCh38/barcodes";
	string shapePath = "/GRCh38/shape";

	// Uncomment these variables to use SZIP compression
	unsigned szip_options_mask;
	unsigned szip_pixels_per_block;

	//create a file
	file_id = H5Fcreate(aimFilePath.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	//create a group
	group_id = H5Gcreate(file_id, groupPath.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	plist_id = H5Pcreate(H5P_DATASET_CREATE);
	cdims[0] = chunk;
	status = H5Pset_chunk(plist_id, rank, cdims);

	if (method == "s") {//szip
		szip_options_mask = H5_SZIP_NN_OPTION_MASK;
		szip_pixels_per_block = 16;
		status = H5Pset_szip(plist_id, szip_options_mask, szip_pixels_per_block);
	}
	else {//zlib
		status = H5Pset_deflate(plist_id, 6);
	}

	//write data
	dims[0] = data_count;
	dataspace_id = H5Screate_simple(rank, dims, NULL);
	dataset_id = H5Dcreate2(file_id, dataPath.c_str(), H5T_STD_I32LE,
		dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	cout << "data writed" << endl;

	//write indptr
	dims[0] = cell_count + 1;
	dataspace_id = H5Screate_simple(rank, dims, NULL);
	dataset_id = H5Dcreate2(file_id, indptrPath.c_str(), H5T_STD_I64LE,
		dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, indptr);
	cout << "indptr writed" << endl;

	//write indices
	dims[0] = data_count;
	dataspace_id = H5Screate_simple(rank, dims, NULL);
	dataset_id = H5Dcreate2(file_id, indicesPath.c_str(), H5T_STD_I64LE,
		dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, indices);
	cout << "indices writed" << endl;

	//write genes
	strGene = H5Tcopy(H5T_C_S1);
	H5Tset_size(strGene, str_genes_length);
	dims[0] = gene_count;
	char *para_genes;
	para_genes = new char[gene_count * str_genes_length];
	for (int i = 0; i < gene_count; i++) {
		strncpy(para_genes + i * str_genes_length, genes[i], str_genes_length);
	}
	//printf("%s", para_genes);
	dataspace_id = H5Screate_simple(rank, dims, NULL);
	dataset_id = H5Dcreate2(file_id, genesPath.c_str(), strGene,
		dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, strGene, H5S_ALL, H5S_ALL, H5P_DEFAULT, para_genes);
	delete[] para_genes;
	cout << "genes writed" << endl;

	//write gene_names
	strGeneName = H5Tcopy(H5T_C_S1);
	str_gene_names_len = 20;
	H5Tset_size(strGeneName, str_gene_names_len);
	dims[0] = gene_count;
	char *para_gene_names;
	para_gene_names = new char[gene_count * str_gene_names_len];
	for (int i = 0; i < gene_count; i++) {
		strncpy(para_gene_names + i * str_gene_names_len, gene_names[i], str_gene_names_len);
	}

	dataspace_id = H5Screate_simple(rank, dims, NULL);
	dataset_id = H5Dcreate2(file_id, geneNamesPath.c_str(), strGeneName,
		dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, strGeneName, H5S_ALL, H5S_ALL, H5P_DEFAULT, para_gene_names);
	delete[] para_gene_names;
	cout << "gene_names writed" << endl;

	//write barcodes
	strCell = H5Tcopy(H5T_C_S1);
	H5Tset_size(strCell, str_barcodes_len);

	dims[0] = cell_count;
	char *para_barcodes;
	para_barcodes = new char[cell_count * str_barcodes_len];
	for (int i = 0; i < cell_count; i++) {
		strncpy(para_barcodes + i * str_barcodes_len, barcodes[i], str_barcodes_len);
	}
	dataspace_id = H5Screate_simple(1, dims, NULL);
	dataset_id = H5Dcreate2(file_id, barcodesPath.c_str(), strCell,
		dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, strCell
		, H5S_ALL, H5S_ALL, H5P_DEFAULT, para_barcodes);
	delete[] para_barcodes;
	cout << "barcodes writed" << endl;

	//write shape
	dims[0] = 2;
	dataspace_id = H5Screate_simple(rank, dims, NULL);
	dataset_id = H5Dcreate2(file_id, shapePath.c_str(), H5T_STD_I32LE,
		dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, shape);
	cout << "shape writed" << endl;

	H5Fclose(file_id);

}

void SparseMatrix::write_norm_data(string write_path, string norm_type, int chunk, string method)
{
	hid_t    file_id, dataset_id, group_id, dataspace_id;
	hid_t    plist_id;
	hid_t    strGene, strCell, strGeneName;
	size_t   nelmts;
	unsigned flags, filter_info;
	H5Z_filter_t filter_type;

	herr_t   status;
	hsize_t  dims[1];
	hsize_t  cdims[1];

	int      idx;
	int      i, j, numfilt;
	int shape[] = { gene_count,cell_count };

	string groupPath = "/GRCh38";
	string dataPath = "/GRCh38/data";
	string genesPath = "/GRCh38/genes";
	string geneNamesPath = "/GRCh38/gene_names";
	string indicesPath = "/GRCh38/indices";
	string indptrPath = "/GRCh38/indptr";
	string barcodesPath = "/GRCh38/barcodes";
	string shapePath = "/GRCh38/shape";

	// Uncomment these variables to use SZIP compression
	unsigned szip_options_mask;
	unsigned szip_pixels_per_block;
	int rank = 1;
	//create a file
	file_id = H5Fcreate(write_path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	//create a group
	group_id = H5Gcreate(file_id, groupPath.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	plist_id = H5Pcreate(H5P_DATASET_CREATE);
	cdims[0] = chunk;
	status = H5Pset_chunk(plist_id, rank, cdims);

	if (method == "s") {//szip
		szip_options_mask = H5_SZIP_NN_OPTION_MASK;
		szip_pixels_per_block = 16;
		status = H5Pset_szip(plist_id, szip_options_mask, szip_pixels_per_block);
	}
	else {//zlib
		status = H5Pset_deflate(plist_id, 6);
	}

	//write data
	if (norm_type == "rank")
	{
		dims[0] = data_count;
		dataspace_id = H5Screate_simple(rank, dims, NULL);
		dataset_id = H5Dcreate2(file_id, dataPath.c_str(), H5T_STD_I16LE,
			dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
		status = H5Dwrite(dataset_id, H5T_NATIVE_USHORT, H5S_ALL, H5S_ALL, H5P_DEFAULT, rankData);
		cout << "data writed" << endl;
	}
	else if (norm_type == "log")
	{
		dims[0] = data_count;
		dataspace_id = H5Screate_simple(rank, dims, NULL);
		dataset_id = H5Dcreate2(file_id, dataPath.c_str(), H5T_IEEE_F64LE,
			dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
		status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, log_normalize_data);
		cout << "data writed" << endl;
	}
	else if (norm_type == "none")
	{
		dims[0] = data_count;
		dataspace_id = H5Screate_simple(rank, dims, NULL);
		dataset_id = H5Dcreate2(file_id, dataPath.c_str(), H5T_STD_I32LE,
			dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
		status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
		cout << "data writed" << endl;
	}
	else if (norm_type == "qqnorm")
	{
		dims[0] = data_count;
		dataspace_id = H5Screate_simple(rank, dims, NULL);
		dataset_id = H5Dcreate2(file_id, dataPath.c_str(), H5T_IEEE_F64LE,
			dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
		status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, qqNormedData);
		cout << "data writed" << endl;
	}
	//write indices
	dims[0] = data_count;
	dataspace_id = H5Screate_simple(rank, dims, NULL);
	dataset_id = H5Dcreate2(file_id, indicesPath.c_str(), H5T_STD_I64LE,
		dataspace_id, H5P_DEFAULT, plist_id, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, indices);
	cout << "indices writed" << endl;

	//write indptr
	dims[0] = cell_count + 1;
	dataspace_id = H5Screate_simple(rank, dims, NULL);
	dataset_id = H5Dcreate2(file_id, indptrPath.c_str(), H5T_STD_I64LE,
		dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, indptr);
	cout << "indptr writed" << endl;

	//write genes
	strGene = H5Tcopy(H5T_C_S1);
	H5Tset_size(strGene, str_genes_length);
	dims[0] = gene_count;
	char *para_genes;
	para_genes = new char[gene_count * str_genes_length];
	for (int i = 0; i < gene_count; i++) {
		strncpy(para_genes + i * str_genes_length, genes[i], str_genes_length);
	}
	//printf("%s", para_genes);
	dataspace_id = H5Screate_simple(rank, dims, NULL);
	dataset_id = H5Dcreate2(file_id, genesPath.c_str(), strGene,
		dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, strGene, H5S_ALL, H5S_ALL, H5P_DEFAULT, para_genes);
	delete[] para_genes;
	cout << "genes writed" << endl;

	//write gene_names
	strGeneName = H5Tcopy(H5T_C_S1);
	str_gene_names_len = 20;
	H5Tset_size(strGeneName, str_gene_names_len);
	dims[0] = gene_count;
	char *para_gene_names;
	para_gene_names = new char[gene_count * str_gene_names_len];
	for (int i = 0; i < gene_count; i++) {
		strncpy(para_gene_names + i * str_gene_names_len, gene_names[i], str_gene_names_len);
	}

	dataspace_id = H5Screate_simple(rank, dims, NULL);
	dataset_id = H5Dcreate2(file_id, geneNamesPath.c_str(), strGeneName,
		dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, strGeneName, H5S_ALL, H5S_ALL, H5P_DEFAULT, para_gene_names);
	delete[] para_gene_names;
	cout << "gene_names writed" << endl;

	//write barcodes
	strCell = H5Tcopy(H5T_C_S1);
	H5Tset_size(strCell, str_barcodes_len);

	dims[0] = cell_count;
	char *para_barcodes;
	para_barcodes = new char[cell_count * str_barcodes_len];
	for (int i = 0; i < cell_count; i++) {
		strncpy(para_barcodes + i * str_barcodes_len, barcodes[i], str_barcodes_len);
	}
	dataspace_id = H5Screate_simple(1, dims, NULL);
	dataset_id = H5Dcreate2(file_id, barcodesPath.c_str(), strCell,
		dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, strCell
		, H5S_ALL, H5S_ALL, H5P_DEFAULT, para_barcodes);
	delete[] para_barcodes;
	cout << "barcodes writed" << endl;

	//write shape
	dims[0] = 2;
	dataspace_id = H5Screate_simple(rank, dims, NULL);
	dataset_id = H5Dcreate2(file_id, shapePath.c_str(), H5T_STD_I32LE,
		dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, shape);
	cout << "shape writed" << endl;

	if (norm_type == "rank")
	{
		dims[0] = cell_count;
		dataspace_id = H5Screate_simple(rank, dims, NULL);
		dataset_id = H5Dcreate2(file_id, "/GRCh38/zero_value", H5T_STD_I32LE,
			dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, rank_zero);
		cout << "zero_value writed" << endl;
	}
	if (norm_type == "log")
	{
		dims[0] = cell_count;
		dataspace_id = H5Screate_simple(rank, dims, NULL);
		dataset_id = H5Dcreate2(file_id, "/GRCh38/size_factor", H5T_STD_I32LE,
			dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, size_factor);
		cout << "size_factor writed" << endl;
	}
	if (norm_type == "qqnorm")
	{
		dims[0] = cell_count;
		dataspace_id = H5Screate_simple(rank, dims, NULL);
		dataset_id = H5Dcreate2(file_id, "/GRCh38/zero_value", H5T_STD_I32LE,
			dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, qqNormedZero);
		cout << "zero_value writed" << endl;
	}
	H5Fclose(file_id);
}


void SparseMatrix::maskingData(int mask_probability, string write_path, string log_path,int seed, string write_type) {
	long long cnt = 0;
	long long start, end;
	srand(seed);
	long long* mask_column,mask_cnt;
	mask_column = new long long[gene_count];
	bool flag = false;
	ofstream ofn(log_path);
	for (int i = 0; i < cell_count; i++) {
		mask_cnt=0;
		start = indptr[i];
		end = indptr[i + 1];
		indptr[i] = indptr[i] - cnt;
		for (long long j = start; j < end; j++) {
			if (rand() % 100 < mask_probability) {
				data[j] = 0;
				data_count--;
				cnt++;
				mask_column[mask_cnt] = indices[j];
				mask_cnt++;
			}
			else {
				data[j - cnt] = data[j];
				indices[j - cnt] = indices[j];
			}
		}
		//sort(mask_column,mask_column+mask_cnt);
		for(int j = 0;j<mask_cnt;j++){
			//if(flag)ofn<<",";
			//else flag=true;
			ofn << i * gene_count + mask_column[j] + 1 << "\n";
		}
	}
	ofn.close();
	delete mask_column;
	indptr[cell_count] = data_count;
	if (write_type == "h5") {
		h5Compressed(write_path,"s",1000,1);
		//write2HDF5(write_path);
	}
	else {
		write2CSV(write_path, "original");
	}
}


