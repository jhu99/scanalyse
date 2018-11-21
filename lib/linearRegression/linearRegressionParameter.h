#ifndef LINEARREGRESSIONOARAMETERS
#define LINEARREGRESSIONOARAMETERS
#include<iostream>
#include<cmath>
#include<gsl/gsl_cdf.h>
using namespace std;
class LinearRegressionParameter {
	double c0;
	double c1;
	double cov00;
	double cov01;
	double cov11;
	double sumsq;
	double sumtot;
	int n;
	double r2;
	double t;
	double pvalue;
public:
	LinearRegressionParameter();
	LinearRegressionParameter(double c0, double c1, double cov00, double cov01, double cov11, double sumsq,double sumtot, int n);
	~LinearRegressionParameter();

	void setC0(double c0);
	void setC1(double c1);
	void setCov00(double cov00);
	void setCov01(double cov01);
	void setCov11(double cov11);
	void setSumsq(double sumsq);
	void setSumtot(double sumtot);
	void setN(int n);
	void setR2(double r2);
	void setT(double t);
	void setPvalue(double pvalue);

	double getC0();
	double getC1();
	double getCov00();
	double getCov01();
	double getCov11();
	double getSumsq();
	double getSumtot();
	double getN();
	double getR2();
	double getT();
	double getPvalue();

	void print();
};
#endif // !LINEARREGRESSIONOARAMETERS

