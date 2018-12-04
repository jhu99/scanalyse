#ifndef RANKNORMALIZE
#define RANKNORMALIZE
#include"HDF5Reader.h"
#include"geneInfo.h"
#include <iostream>
#include <algorithm>
#include <functional>
#include <chrono>    
#include <thread> 

using namespace std;
using std::placeholders::_1;
using std::placeholders::_2;

class rankNormalize
{
	long long n;
	HDF5reader hr;
	geneInfo *geneInfos;
	double *rank;
public:
	rankNormalize();
	rankNormalize(HDF5reader hr);
	~rankNormalize();

	void setN(long long n);
	void setHr(HDF5reader hr);
	void setGeneInfos(geneInfo *geneInfos);
	void setRank(double *rank);

	long long getN();
	HDF5reader getHr();
	geneInfo* getGeneInfos();
	double* getRank();

	void setGeneInfoByHr();
	bool cmp1(geneInfo x, geneInfo y);
	bool cmp2(geneInfo x, geneInfo y);
	void sortByData(int begin, int len);
	void sortByIndices(int begin, int len);
	//void ranksThread(int i);
	void ranks();
	void print();
};

#endif // !RANKNORMALIZE



