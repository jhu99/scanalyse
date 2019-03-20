#include"SparseMatrix/SparseMatrix.h"
#include<vector>
#include<cstdlib>
using namespace std;
using namespace H5;

int main(int argc,const char** argv) {
	SparseMatrix m;
	vector<string> paths;
	bool log;
	printf("%s\n",argv[argc-1]);
	
	
	if(strcmp(argv[argc-2],"true")==0){
		log = true;	
	}	
	else{
		log = false;	
	}
	if(log){
		for(int i=2;i< argc-2;i++){
			paths.push_back(argv[i]);
		}
		m.mergeDate(paths,1,argv[argc-1]);
	}
	else{
		for(int i=2;i< argc;i++){
			paths.push_back(argv[i]);
		}
		m.mergeDate(paths);
	}
	
	m.h5Compressed(argv[1], "s", 1000, 1);
}
