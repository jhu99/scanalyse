#ifndef GENEEXPRESSIONTOP
#define GENEEXPRESSIONTOP
#include"geneCount.h"
#include"SparseMatrix/SparseMatrix.h"
#include<algorithm>
class geneExpressionTop {
	int topNum = 0, num=0;
	SparseMatrix sm;
	geneCount * gc;
	string *top;
	long long *top_index;
public:
	geneExpressionTop(SparseMatrix sm,int topNum);
	~geneExpressionTop();
	void setSparseMatrix(SparseMatrix sm);
	void setTopNum(int topNum);
	void setNum(int num);
	void setGeneCount(geneCount *gc);
	void setTop(string *top);

	SparseMatrix getSparseMatrix();
	int getTopNum();
	int getNum();
	geneCount* getGeneCount();
	string* getTop();
	long long *get_top_index();

	void geneSort();
	void printTop();
};
#endif // !GENEEXPRESSIONTOP

