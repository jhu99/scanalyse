#ifndef LINEARREGRESSIONOARAMETERS
#define LINEARREGRESSIONOARAMETERS
#include<iostream>
using namespace std;
class LinearRegressionParameter {
	double c0;
	double c1;
	double cov00;
	double cov01;
	double cov11;
	double sumsq;
public:
	LinearRegressionParameter();
	LinearRegressionParameter(double c0, double c1, double cov00, double cov01, double cov11, double sumsq);
	~LinearRegressionParameter();

	void setC0(double c0);
	void setC1(double c1);
	void setCov00(double cov00);
	void setCov01(double cov01);
	void setCov11(double cov11);
	void setSumsq(double sumsq);

	double getC0();
	double getC1();
	double getCov00();
	double getCov01();
	double getCov11();
	double getSumsq();

	void print();
};
#endif // !LINEARREGRESSIONOARAMETERS
