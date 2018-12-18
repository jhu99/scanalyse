#include"qqNorm/qqNorm.h"
#include<iostream> 
#include <string>  
#include"qqNorm/caculateInterface.h"
using namespace std;

int main(int argc, const char ** argv)
{
	qqNorm qq;
	unsigned short *test;
	test = new unsigned short[10];
	for (int i = 0; i <8; i++)
	{
		test[i] = i;
	}
	test[8] = 0; test[9] = 0;
	double *result = qq.caculateTheoryQuantiles(test, 10);
	for (int i = 0; i < 10; i++)
	{
		cout << result[i] << " ";
	}
	cin.get();
	cin.get();
	return 0;
}