#include"geneExpressionTop.h"

geneExpressionTop::geneExpressionTop(SparseMatrix sm,int topNum) {
	num = sm.get_gene_count();
	int n = sm.get_data_count();
	this->topNum = topNum;
	gc = new geneCount[num];
	top = new string[topNum];
	top_index = new long long[topNum];
	for (int i = 0; i < num; i++) {
		gc[i].setGeneName(sm.get_genes()[i]);
		gc[i].setIndex(i);
	}
	for (int i = 0; i < n; i++) {
		gc[sm.get_indices()[i]].countAdd(sm.get_data()[i]);
	}
}
geneExpressionTop::~geneExpressionTop() {
	delete[] gc;
	delete[] top;
}

void geneExpressionTop::setSparseMatrix(SparseMatrix sm) {
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

SparseMatrix geneExpressionTop::getSparseMatrix() {
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

long long * geneExpressionTop::get_top_index()
{
	return top_index;
}

void geneExpressionTop::geneSort() {
	sort(gc,gc+num);
	for (int i = 0; i < topNum; i++) {
		top[i] = gc[i].getGeneName();
		top_index[i] = gc[i].getIndex();
	}
}

void geneExpressionTop::printTop() {
	for (int i = 0; i < topNum; i++) {
		cout << i << " " << gc[i].getGeneName() << " " << gc[i].getCount() << endl;
	}
}
