#ifndef LINEARREGRESSION_H
#define LINEARREGRESSION_H

#include<iostream>
#include<cstdio>
#include<cstdlib>
#include<chrono>    
#include<thread>  
#include<mutex>
#include<algorithm>
#include<unordered_map>
#include<gsl/gsl_fit.h>
#include<gsl/gsl_statistics_double.h>
#include<fstream>
#include"cell/cell.h"
#include"linearRegressionParameter.h"
#
using namespace std;

template <class T> class LinearRegression {
	int n, p;
	mutex mut;
	Cells<T> cells;
	LinearRegressionParameter **linearRegressionParameters;

public:
	LinearRegression();
	LinearRegression(Cells<T> cells);
	~LinearRegression();

	void setCell(Cells<T> cells);
	void setLinearRegressionParameter(LinearRegressionParameter **linearRegressionParameters);

	Cells<T> getCell();
	LinearRegressionParameter **getLinearRegressionParameter();

	void calculate(int i);
	void calculateLinearRegression(int n);
	void writeFile(string path);


};



#endif // !LINEARREGRESSION_H
