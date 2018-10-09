#include"cell.h"
#include<iostream>
using namespace std;
int main() {
	Cells cells;
	cells.readFile();
	//cells.findCell("AAACCTGAGCCACTAT.1");
	cells.findGene("ENSMUSG00000037221");
	//	cells.findCellAndGene("AAACCTGAGCCACTAT.1", "ENSMUSG00000015312");
	cin.get();
	cin.get();
	return 0;
}
