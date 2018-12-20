
#include "geneVariationTop.h"

int main() {
	SparseMatrix sm;
	int test = sm.readHDF5File("ica_cord_blood_h5.h5");
	geneVariationTop GVT(sm, 5000);
	GVT.geneSort();
	GVT.printTop();
	system("pause");
	return 0;
}