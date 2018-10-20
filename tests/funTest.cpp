#include"fun.h"
#include"cell/cell.cpp"
#include<iostream>
#include <fstream>  
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
