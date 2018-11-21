#pragma once
#include<list>
#include <iostream>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include "caculateInterface.h"

class qqNorm
{
public:
	double **theoryQ;
	qqNorm() {

	}
	~qqNorm()
	{

	}

	double * caculateTheoryQuantiles(short *p,int size);

};