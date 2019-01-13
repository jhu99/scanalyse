#ifndef GENEVARIATIONTOP
#define GENEVARIATIONTOP
#include "geneVariation.h"
#include "SparseMatrix/SparseMatrix.h"
#include "geneTop/geneCount.h"
#include <algorithm>
class geneVariationTop {
	int topNum = 0, gnum = 0;
	SparseMatrix sm;
	geneCount *gc;
	geneVariation *gv;
	string * top;
	long long* top_index;
public:
	geneVariationTop(SparseMatrix sm, int topNum);
	~geneVariationTop();
	void set_SparseMatrix(SparseMatrix sm);
	void set_TopNum(int topNum);
	void set_Num(int gnum);
	void set_geneVariation(geneVariation *gv);
	void set_Top(string *top);

	SparseMatrix get_SparseMatrix();
	int get_TopNum();
	int get_Num();
	geneVariation* get_geneVariation();
	string* get_Top();
	long long *get_top_index();

	void geneSort();
	void printTop();

};
#endif
