#include"SparseMatrix/SparseMatrix.h"
#include<vector>
using namespace std;
using namespace H5;

int main() {
	SparseMatrix m;
	vector<string> paths;
	string path = "./data/ica_cord_blood_h5.h5";
	//string path = "D:/test.h5";
	paths.push_back(path);
	//paths.push_back(path);
	m.mergeDate(paths);
	m.h5Compressed("./data/test.h5", "s", 1000, 1);
}
