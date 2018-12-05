#include"qqNorm/qqNorm.h"
#include<iostream> 
#include <string>  
#include"qqNorm/caculateInterface.h"
using namespace std;

int main()
{
	qqNorm qq;
	unsigned short **test;
	int columCount=12, rowCount=2;
	test = new unsigned short*[rowCount];
	for (int i = 0; i < rowCount; i++)
	{
		test[i] =new unsigned short [columCount];
	}
	for (int i = 0; i < rowCount; i++)
	{
		for (int j = 0; j < columCount; j++)
		{
			test[i][j] = (short)(j+1);
		}
	}
	cacuTheoryQ c;
	double** result = c.cacuQByColumn(test,rowCount,columCount);
	cin.get();
	cin.get();
	return 0;
}