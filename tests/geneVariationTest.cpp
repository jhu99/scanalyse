
#include "geneVariationTop/geneVariationTop.h"

int main() {
	SparseMatrix sm;
	string type = "original";
	int test = sm.readHDF5File("./data/ica_cord_blood_h5.h5", type);
	geneVariationTop GVT(sm, 5000);
	GVT.geneSort();
	GVT.printTop();
	sm.deleteSparseMatrix("original");
	cin.get();
	cin.get();
	return 0;
}