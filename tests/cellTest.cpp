#include"cell/cell.cpp"
#include<iostream>
using namespace std;

int main() {
	Cells<short> cells(5005,4001);
	cells.readFile("scimpute-0.01-1-data.csv");
	//cells.findCell("AAACCTGAGCCACTAT.1");
	cells.findGene("ENSMUSG00000037221");
	//cells.findCellAndGene("AAACCTGAGCCACTAT.1", "ENSMUSG00000015312");
	cout<<"end"<<endl;
	return 0;
}
