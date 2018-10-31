#include"linearRegressionParameter.h"

LinearRegressionParameter::LinearRegressionParameter() {

}
LinearRegressionParameter::LinearRegressionParameter(double c0, double c1, double cov00, double cov01, double cov11, double sumsq) {
	this->c0 = c0;
	this->c1 = c1;
	this->cov00 = cov00;
	this->cov01 = cov01;
	this->cov11 = cov11;
	this->sumsq = sumsq;

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

void LinearRegressionParameter::print() {
	cout << c0 << " " << c1 << endl;
}
