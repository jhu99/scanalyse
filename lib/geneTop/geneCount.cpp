#include"geneCount.h"

void geneCount::setGeneName(string geneName) {
	this->geneName = geneName;
}
void geneCount::setCount(int count) {
	this->count = count;
}

void geneCount::setIndex(long long index)
{
	this->index = index;
}

string geneCount::getGeneName() {
	return geneName;
}
int geneCount::getCount() {
	return count;
}
long long geneCount::getIndex()
{
	return index;
}
void geneCount::countAdd(int num) {
	count += num;
}



bool geneCount::operator < (const geneCount &x)const
{
	return x.count < count;
}
