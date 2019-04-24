#include<vector>
#include<cstdlib>
#include<iostream>
#include"fun/fun.h"
#include"SparseMatrix/SparseMatrix.h" 
using namespace std;
using namespace Scanalyse;
using namespace std;
using namespace H5;

int main(int argc,const char** argv) {
	cout<<"merge H5 begin"<<endl;
	SparseMatrix m;
	vector<string> paths;
	bool log;
	string out_path = argv[1];
	
	
	if(strcmp(argv[argc-4],"true")==0){
		log = true;	
	}	
	else{
		log = false;	
	}
	if(log){
		for(int i=2;i< argc-4;i++){
			paths.push_back(argv[i]);
		}
		m.mergeDate(paths,1,argv[argc-3]);
	}
	else{
		for(int i=2;i< argc;i++){
			paths.push_back(argv[argc-2]);
		}
		m.mergeDate(paths);
	}
	
	m.h5Compressed(out_path+".h5", "s", 1000, 1);
	
	cout<<"H52CSV begin"<<endl;

	m.readHDF5File(out_path+".h5", argv[argc-2]);
	m.write2CSV(out_path+".csv",argv[argc-2]);
	
	cout<<"CSV_T begin"<<endl;
	Fun fun;
	fun = Fun();
	string type = argv[argc-1];
	fun.transfer_matrix(out_path+".csv", out_path+"_t.csv", type);
	return 0;
}
