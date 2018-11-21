#include"qqNorm.h"
using namespace std;

double* qqNorm::caculateTheoryQuantiles(short* p,int size)
{
	typedef std::numeric_limits<double> Info;
	double const NAN_d = Info::quiet_NaN();
	list<double> result;
	float *inputP = new float[size];
	double *theoryQuantiles = new double[size];
	if (size > 10)
	{
		for (int i = 0; i < size; i++)
		{
			inputP[i] = (p[i]-0.5) / size;
		}

	}
	else
	{
		for (int i = 0; i < size; i++)
		{
			inputP[i] = (p[i] - 0.375) / (size +0.25);
		}
	}
	for (int i = 0; i < size; i++)
	{
		if (inputP[i] < 0 || inputP[i]>1)
		{
			p[i] = NAN_d;
		}
	}
	for (int i = 0; i < size; i++)
	{
		theoryQuantiles[i] = gsl_cdf_ugaussian_Pinv(inputP[i]);
	}
	
	return theoryQuantiles;
}