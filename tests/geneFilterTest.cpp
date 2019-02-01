#include"geneFilter/geneFilter.h"
int main(int argc, const char ** argv)
{
	geneFilter gf;
	gf.filt2H5File(argv[1], argv[2], argv[3]);
	cin.get();
	cin.get();
}