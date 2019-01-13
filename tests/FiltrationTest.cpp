#include"Filtration/Filtration.h"
int main(int argc, const char ** argv)
{
	SparseMatrix sm;
	string path_read = argv[1];
	string type = "original";
	sm.readHDF5File(path_read, type);
	Filtration f(sm);
	f.filtGeneAndCell(3000, 100, 2);
	string deleteType= "orignal";
	
	f.printFiltResult();
	sm.deleteSparseMatrix(deleteType);
	cin.get();
	cin.get();

}