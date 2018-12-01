#include"geneCount.h"

void geneCount::setGeneName(string geneName) {
	this->geneName = geneName;
}
void geneCount::setCount(int count) {
	this->count = count;
}

string geneCount::getGeneName() {
	return geneName;
}
int geneCount::getCount() {
	return count;
}
void geneCount::countAdd(int num) {
	count += num;
}



bool geneCount::operator < (const geneCount &x)const
{
	return x.count < count;
}
