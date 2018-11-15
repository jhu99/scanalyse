#include<iostream>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_cdf.h>
#include<cmath>
#include <cstdio>
#include <cstdlib>
#include <chrono>    
#include <thread>  
#include<gsl/gsl_statistics_double.h>

#include<gsl/gsl_fit.h>
#include"linearRegression/linearRegressionParameter.h"
using namespace std;

int main() {
	
	double x[] = { 4,8,9,8,7,12,6,10,6,9 };
	double y[] = { 9,20,22,15,17,23,18,25,10,20 };
	double c0, c1, cov00,cov01, cov11, sumsq, sumtot=0;
	//sumtot = gsl_stats_tss(y, 1, 10);
	sumtot = gsl_stats_variance(y,1,10)*9;
	cout << sumtot<<":"<< sumtot << endl;
	gsl_fit_linear(x, 1, y, 1, 10, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);
	LinearRegressionParameter lpr = LinearRegressionParameter(c0, c1, cov00, cov01, cov11, sumsq,sumtot, 10);
	lpr.print();
	
	
	cin.get();
	cin.get();
	return 0;
}
