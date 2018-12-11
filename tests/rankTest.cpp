#include<iostream>
#include"SparseMatrix/SparseMatrix.h"
#include"rank/rankNormalize.h"
#include<ctime>
using namespace std;

int main() {
	clock_t start,end;
	HDF5reader hr;
	hr.readHDF5File("./data/ica_cord_blood_h5.h5");
	start=clock();
	rankNormalize rn(hr);
	rn.ranks(5);
	end=clock();
	double endtime=(double)(end-start)/CLOCKS_PER_SEC;
	//cout<<"Total time:"<<endtime<<endl;		
	//rn.print();
	cout<<"Total time:"<<endtime<<endl;		
	hr.deleteSparseMatrix();
}
