#include"geneInfo.h"

geneInfo::geneInfo(){}
geneInfo::geneInfo(int data, long long indices) {
	this->data = data;
	this->indices = indices;
}
geneInfo::~geneInfo(){}

int  geneInfo::getData() {
	return data;
}
long long geneInfo::getIndices() {
	return indices;
}
unsigned short geneInfo::getRank() {
	return rank;
}

void geneInfo::setData(int data) {
	this->data = data;
}
void geneInfo::setIndices(long long indices) {
	this->indices = indices;
}
void geneInfo::setRank(unsigned short rank) {
	this->rank = rank;
}
