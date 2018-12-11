#ifndef RANKNORMALIZE
#define RANKNORMALIZE
#include"SparseMatrix/SparseMatrix.h"
#include"geneInfo.h"
#include <iostream>
#include <algorithm>
#include <functional>
#include <chrono>    
#include <thread> 
#include<mutex>

using namespace std;
using std::placeholders::_1;
using std::placeholders::_2;

class rankNormalize
{
	mutex mut;
	long long n;
	SparseMatrix hr;
	geneInfo *geneInfos;
	unsigned short *rank;
public:
	rankNormalize();
	rankNormalize(SparseMatrix hr);
	~rankNormalize();

	void setN(long long n);
	void setHr(SparseMatrix hr);
	void setGeneInfos(geneInfo *geneInfos);
	void setRank(unsigned short *rank);

	long long getN();
	SparseMatrix getHr();
	geneInfo* getGeneInfos();
	unsigned short* getRank();

	void setGeneInfoByHr();
	bool cmp1(geneInfo x, geneInfo y);
	bool cmp2(geneInfo x, geneInfo y);
	void sortByData(int begin, int len);
	void sortByIndices(int begin, int len);
	void ranksThread(int i);
	void ranks(int nt);
	void print();
};

#endif // !RANKNORMALIZE


