#include"linearRegressionParameter.h"

LinearRegressionParameter::LinearRegressionParameter() {

}
LinearRegressionParameter::LinearRegressionParameter(double c0, double c1, double cov00, double cov01, double cov11, double sumsq, double sumtot, int n) {
	this->c0 = c0;
	this->c1 = c1;
	this->cov00 = cov00;
	this->cov01 = cov01;
	this->cov11 = cov11;
	this->sumsq = sumsq;
	this->sumtot = sumtot;
	this->n = n;
	r2 = 1 - (sumsq / sumtot);
	t = sqrt(r2 / (1 - r2)*(n - 2));
	pvalue = gsl_cdf_tdist_Q(t, n - 1) * 2;

}
LinearRegressionParameter::~LinearRegressionParameter() {

}

void LinearRegressionParameter::setC0(double c0) {
	this->c0 = c0;
}
void LinearRegressionParameter::setC1(double c1) {
	this->c1 = c1;
}
void LinearRegressionParameter::setCov00(double cov00) {
	this->cov00 = cov00;
}
void LinearRegressionParameter::setCov01(double cov01) {
	this->cov01 = cov01;
}
void LinearRegressionParameter::setCov11(double cov11) {
	this->cov11 = cov11;
}
void LinearRegressionParameter::setSumsq(double sumsq) {
	this->sumsq = sumsq;
}
void LinearRegressionParameter::setSumtot(double sumtot) {
	this->sumtot = sumtot;
}
void LinearRegressionParameter::setN(int n) {
	this->n = n;
}
void LinearRegressionParameter::setR2(double r2) {
	this->r2 = r2;
}
void LinearRegressionParameter::setT(double t) {
	this->t = t;
}
void LinearRegressionParameter::setPvalue(double pvalue) {
	this->pvalue = pvalue;
}

double LinearRegressionParameter::getC0() {
	return c0;
}
double LinearRegressionParameter::getC1() {
	return c1;
}
double LinearRegressionParameter::getCov00() {
	return cov00;
}
double LinearRegressionParameter::getCov01() {
	return cov01;
}
double LinearRegressionParameter::getCov11() {
	return cov11;
}
double LinearRegressionParameter::getSumsq() {
	return sumsq;
}
double LinearRegressionParameter::getSumtot() {
	return sumtot;
}
double LinearRegressionParameter::getN() {
	return n;
}
double LinearRegressionParameter::getR2() {
	return r2;
}
double LinearRegressionParameter::getT() {
	return t;
}
double LinearRegressionParameter::getPvalue() {
	return pvalue;
}

void LinearRegressionParameter::print() {
	cout << "c0:" << c0 << endl;
	cout << "c1:" << c1 << endl;
	cout << "cov00:" << cov00 << endl;
	cout << "cov01:" << cov01 << endl;
	cout << "cov11:" << cov11<< endl;
	cout << "sumsq:" << sumsq << endl;
	cout << "sumtot:" << sumtot << endl;
	cout << "n:" << n << endl;
	cout << "r2:" << r2 << endl;
	cout << "t:" << t << endl;
	cout << "pvalue:" << pvalue << endl;
}
