#include "caculateInterface.h"
#include "qqNorm.h"
double ** cacuTheoryQ::cacuQByColumn(short ** rankMatrix, int rowCount, int columnCount)
{
	double** theoryQMatix;
	theoryQMatix = new double*[rowCount];
	for (int i = 0; i < rowCount; i++)
	{
		theoryQMatix[i] = new double[columnCount];
	}
	qqNorm q;
	double* paraRow;
	for (int i = 0; i < rowCount; i++)
	{
		theoryQMatix[i]=q.caculateTheoryQuantiles(rankMatrix[i],columnCount);
	}

	for (int i = 0; i < rowCount; i++)
	{
		for (int j = 0; j < columnCount; j++)
		{
			cout << theoryQMatix[i][j]<<" ";
		}
		cout << endl;
	}

	return theoryQMatix;
}
