#include "geneVariation.h"
#include<cmath>

void geneVariation::set_geneName(string geneName) {
	this->geneName = geneName;
}

void geneVariation::set_Variation(double Variation) {
	this->Variation = Variation;
}

void geneVariation::set_Index(long long index)
{
	this->index = index;
}

string geneVariation::get_geneName() {
	return geneName;
}

double geneVariation::get_Variation() {
	return Variation;
}

long long geneVariation::get_index()
{
	return index;
}

void geneVariation::VarAdd(int num, double avg,int n) {
	Variation += pow((num - avg), 2)/n;
}

bool geneVariation::operator<(const geneVariation &a)const {
	return a.Variation < Variation;
}