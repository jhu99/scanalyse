#include"geneInfo.h"

geneInfo::geneInfo(){}
geneInfo::geneInfo(int data, long long indices) {
	this->data = data;
	this->indices = indices;
}

int  geneInfo::getData() {
	return data;
}
long long geneInfo::getIndices() {
	return indices;
}
double geneInfo::getRank() {
	return rank;
}

void geneInfo::setData(int data) {
	this->data = data;
}
void geneInfo::setIndices(long long indices) {
	this->indices = indices;
}
void geneInfo::setRank(double rank) {
	this->rank = rank;
}
