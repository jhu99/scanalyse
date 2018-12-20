#include "geneVariationTop.h"

geneVariationTop::geneVariationTop(SparseMatrix sm, int topNum) {
	gnum = sm.get_gene_count();
	int n = sm.get_data_count();
	int p = sm.get_cell_count();
	this->topNum = topNum;
	gc = new geneCount[gnum];
	gv = new geneVariation[gnum];
	top = new string[topNum];
	for (int i = 0; i < gnum; i++) {
		gc[i].setGeneName(sm.get_genes()[i]);
		gv[i].set_geneName(sm.get_genes()[i]);
	}
	for (int i = 0; i < n; i++) {
		gc[sm.get_indices()[i]].countAdd(sm.get_data()[i]);// 对同gene不同cell表达值求和
	}
	for (int i = 0; i < n; i++) {
		gv[sm.get_indices()[i]].VarAdd(sm.get_data()[i], (double)gc[sm.get_indices()[i]].getCount() / p, p);
	} //对gene求方差
}

geneVariationTop::~geneVariationTop() {
	delete[] gc;
	delete[] gv;
	delete[] top;
}
void geneVariationTop::set_SparseMatrix(SparseMatrix sm) {
	this->sm = sm;
}
void geneVariationTop::set_TopNum(int topNum) {
	this->topNum = topNum;
}
void geneVariationTop::set_Num(int gnum) {
	this->gnum = gnum;
}
void geneVariationTop::set_geneVariation(geneVariation *gv) {
	this->gv = gv;
}
void geneVariationTop::set_Top(string *top) {
	this->top = top;
}

SparseMatrix geneVariationTop::get_SparseMatrix() {
	return sm;
}
int geneVariationTop::get_TopNum() {
	return topNum;
}
int geneVariationTop::get_Num() {
	return gnum;
}
geneVariation* geneVariationTop::get_geneVariation() {
	return gv;
}
string* geneVariationTop::get_Top() {
	return top;
}

void geneVariationTop::geneSort() {
	sort(gv, gv + gnum);
	for (int i = 0; i < topNum; i++) {
		top[i] = gv[i].get_geneName();
	}
}

void geneVariationTop::printTop() {
	for (int i = 0; i < topNum; i++) {
		cout << i << '\t' << top[i] << '\t' << gv[i].get_Variation() << endl;
	}
}