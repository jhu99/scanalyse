#ifndef GENEEXPRESSIONTOP
#define GENEEXPRESSIONTOP
#include"geneCount.h"
#include"HDF5Reader.h"
#include<algorithm>
class geneExpressionTop {
	int topNum = 0, num=0;
	HDF5reader hr;
	geneCount * gc;
	string *top;
public:
	geneExpressionTop(HDF5reader hr,int topNum);
	~geneExpressionTop();
	void setHDF5reader(HDF5reader hr);
	void setTopNum(int topNum);
	void setNum(int num);
	void setGeneCount(geneCount *gc);
	void setTop(string *top);

	HDF5reader getHDF5reader();
	int getTopNum();
	int getNum();
	geneCount* getGeneCount();
	string* getTop();

	void geneSort();
	void printTop();
};
#endif // !GENEEXPRESSIONTOP

