#include"fun/fun.h"
#include<iostream>  
using namespace std;
using namespace Scanalyse;


int main(int argc, const char **argv)
{
	Fun fun;
	fun = Fun();
	string filepath = argv[1];
	string outpath = argv[2];
	string type = argv[3];
	fun.transfer_matrix(filepath, outpath, type);
	return 0;
}
