#include"fun/fun.h"
#include<iostream>  
using namespace std;
using namespace Scanalyse;


int main(int argc, const char **argv)
{
	Fun fun;
	fun = Fun();
	string filePath = argv[1];
	string outpath = argv[2];
	fun.transfer_matrix(filePath, outpath, argv[3]);
	cin.get();
	cin.get();
	return 0;
}