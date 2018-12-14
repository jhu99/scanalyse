#include"geneExpressionTop.h"

geneExpressionTop::geneExpressionTop(SparseMatrix sm,int topNum) {
	num = sm.get_gene_count();
	int n = sm.get_data_count();
	this->topNum = topNum;
	gc = new geneCount[num];
	top = new string[topNum];
	for (int i = 0; i < num; i++) {
		gc[i].setGeneName(sm.get_genes()[i]);
	}
	for (int i = 0; i < n; i++) {
		gc[sm.get_indices()[i]].countAdd(sm.get_data()[i]);
	}
}
geneExpressionTop::~geneExpressionTop() {
	delete[] gc;
	delete[] top;
	sm.deleteHDF5();
}

void geneExpressionTop::setHDF5reader(SparseMatrix sm) {
	this->sm = sm;
}
void geneExpressionTop::setTopNum(int topNum) {
	this->topNum = topNum;
}
void geneExpressionTop::setNum(int num) {
	this->num = num;
}
void geneExpressionTop::setGeneCount(geneCount *gc) {
	this->gc = gc;
}
void geneExpressionTop::setTop(string *top) {
	this->top = top;
}

HDF5reader geneExpressionTop::getSparseMatrix() {
	return sm;
}
int geneExpressionTop::getTopNum() {
	return topNum;
}
int geneExpressionTop::getNum() {
	return num;
}
geneCount* geneExpressionTop::getGeneCount() {
	return gc;
}
string* geneExpressionTop::getTop() {
	return top;
}

void geneExpressionTop::geneSort() {
	sort(gc,gc+num);
	for (int i = 0; i < topNum; i++) {
		top[i] = gc[i].getGeneName();
	}
}

void geneExpressionTop::printTop() {
	for (int i = 0; i < topNum; i++) {
		cout << i << " " << gc[i].getGeneName() << " " << gc[i].getCount() << endl;
	}
}
