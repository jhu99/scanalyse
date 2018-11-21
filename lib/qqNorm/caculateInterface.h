#pragma once
#include <iostream>

using namespace std;

// base class
class caculateInterface
{
public:
	// provide virtual function for interface
	virtual double** cacuQByColumn(short**rankMatrix, int rowCount, int columnCount)=0;
protected:

};

//subClass
class cacuTheoryQ :caculateInterface
{
public:
	cacuTheoryQ() {

	}
	~cacuTheoryQ() {
	}
	double** cacuQByColumn(short**rankMatrix,int rowCount,int columnCount);
};
