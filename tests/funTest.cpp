#include"fun/fun.h"
#include<iostream>
#include <string>  

using namespace std;
using namespace Scanalyse;

int main()
{
	Fun fun;
	fun = Fun();
	string filePath = "../../cell";
	string outpath = "../../";
	fun.read(filePath);
	fun.CreatAllGeneMap();
	fun.initMergeMatrix();
	fun.mergeMatrixs();
	fun.outToCsvFile(outpath);

	cin.get();
	cin.get();
	return 0;
}
