#include"geneExpressionTop.h"

geneExpressionTop::geneExpressionTop(HDF5reader hr,int topNum) {
	num = hr.get_gene_count();
	int n = hr.get_data_count();
	this->topNum = topNum;
	gc = new geneCount[num];
	top = new string[topNum];
	for (int i = 0; i < num; i++) {
		gc[i].setGeneName(hr.get_genes()[i]);
	}
	for (int i = 0; i < n; i++) {
		gc[hr.get_indices()[i]].countAdd(hr.get_data()[i]);
	}
}
geneExpressionTop::~geneExpressionTop() {
	delete[] gc;
	delete[] top;
	hr.deleteHDF5();
}

void geneExpressionTop::setHDF5reader(HDF5reader hr) {
	this->hr = hr;
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

HDF5reader geneExpressionTop::getHDF5reader() {
	return hr;
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
