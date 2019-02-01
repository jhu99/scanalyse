#include"geneFilter/geneFilter.h"
int main(int argc, const char ** argv)
{
	geneFilter gf;
	gf.createFiltGeneTemplate(argv[1]);
	gf.out2CsvFile(argv[2]);
	cin.get();
	cin.get();
}