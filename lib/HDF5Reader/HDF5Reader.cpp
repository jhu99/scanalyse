#include "HDF5Reader.h"
#define _CRT_SECURE_NO_WARNINGS
using namespace std;
using namespace H5;
int HDF5reader::readHDF5File(string path)
{
	hid_t       fileId, space, datasetId,type,groupId; 
	herr_t      status;

	data = new int[data_count];
	indices = new int[data_count];
	indptr = new int[cell_count+1];
	char* para_barcodes;
	para_barcodes = new char[40 * cell_count];
	char*para_genenames;
	para_genenames = new char[19 * gene_count];
	barcodes = new char*[cell_count];
	for (int i = 0; i < cell_count; i++)
	{
		barcodes[i] = new char[40];
	}
	gene_names = new char*[gene_count];
	for (int i = 0; i < gene_count; i++)
	{
		gene_names[i] = new char[19];
	}
	const char*p = path.data();
	string groupName = "/GRCh38";
	const char*group_name = groupName.data();
	string dataSetName = "/GRCh38/data";
	const char* dataset_name=dataSetName.data();
	string paraCellName;
	
	fileId = H5Fopen(p, H5F_ACC_RDONLY, H5P_DEFAULT);
	status = H5Lget_info(fileId, group_name, NULL, H5P_DEFAULT);
	printf("/GRCh38': ");
	if (status == 0)
	{
		printf("The group exists.\n");
		groupId= H5Gopen2(fileId, group_name, H5P_DEFAULT);
	}
	else
	{
		printf("The group either does NOT exist or some other error occurred.\n");
	}
	string info;
	status = H5Lget_info(fileId, dataset_name,NULL, H5P_DEFAULT);
	printf("/GRCh38/genes': ");
	if (status == 0)
	{
		printf("The dataset exists.\n");
		datasetId = H5Dopen(fileId, dataset_name, H5P_DEFAULT);
	}
	else
	{
		printf("The dataset either does NOT exist or some other error occurred.\n");
	}
	//read cellname
	hid_t str40 = H5Tcopy(H5T_C_S1);
	H5Tset_size(str40, 40);
	datasetId = H5Dopen(fileId, "/GRCh38/barcodes", H5P_DEFAULT);
	status = H5Dread(datasetId, str40, H5S_ALL, H5S_ALL, H5P_DEFAULT, para_barcodes);
	for (int i = 0; i < cell_count; i++)
	{
		strncpy(barcodes[i], para_barcodes+i*40, 40);
	}
	//read genenames
	
	hid_t str_genenames = H5Tcopy(H5T_C_S1);
	H5Tset_size(str_genenames, 19);
	datasetId = H5Dopen(fileId, "/GRCh38/gene_names", H5P_DEFAULT);
	status = H5Dread(datasetId, str_genenames, H5S_ALL, H5S_ALL, H5P_DEFAULT, para_genenames);
	for (int i = 0; i < gene_count; i++)
	{
		strncpy(gene_names[i], para_genenames + i * 19, 19);
	}


	cout << "finish read barcodes" << endl;
	//read data
	datasetId = H5Dopen(fileId, dataset_name, H5P_DEFAULT);
	status = H5Dread(datasetId, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	cout << "finish read data" << endl;
	//read indices
	datasetId = H5Dopen(fileId, "/GRCh38/indices", H5P_DEFAULT);
	status = H5Dread(datasetId, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, indices);
	cout << "finish read indices" << endl;
	//read indptr
	datasetId = H5Dopen(fileId, "/GRCh38/indptr", H5P_DEFAULT);
	status = H5Dread(datasetId, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, indptr);
	cout << "finish read indptr" << endl;
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
	for (int paraPos = indptr[cellPos]; paraPos < indptr[cellPos+1]; paraPos++)
	{
		columPos = indices[paraPos];
		singleCellVector[columPos] = data[paraPos];
	}
	return singleCellVector;
}

char** HDF5reader::get_barcodes()
{
	return barcodes;
}

char ** HDF5reader::get_gene_names()
{
	return gene_names;
}

int * HDF5reader::get_indptr()
{
	return indptr;
}
