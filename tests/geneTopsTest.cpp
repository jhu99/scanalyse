#include<iostream>
#include"geneTop/geneExpressionTop.h"
#include"SparseMatrix/SparseMatrix.h"
using namespace std;

int main() {
	SparseMatrix sm;
	//freopen("C:/Users/ltw/Desktop/1.txt", "w", stdout);
	string path = "./data/ica_cord_blood_h5.h5";
	sm.readHDF5File(path,"original");
	/*for (int i = 0; i < hr.get_gene_count(); i++) {
		cout << i<<" "<<hr.get_genes()[i] << endl;
	}*/
	
	geneExpressionTop geneTops(sm, 500);
	geneTops.geneSort();
	geneTops.printTop();
	cin.get();
	cin.get();
	sm.deleteSparseMatrix("original");
}
