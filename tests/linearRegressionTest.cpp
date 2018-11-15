#include<iostream>
#include <cstdio>
#include <cstdlib>
#include <chrono>    
#include <thread>  
#include"cell/cell.cpp"
#include"linearRegression/linearRegression.cpp"
#include"linearRegression/linearRegressionParameter.h"
using namespace std;

int main() {
	Cells<short> cell(6000, 6000);
	cell.readFile("./data/t.csv");
	cout << cell.getNumToGene()[0] << endl;
	//cell.findCell("AAACCTGAGCCACTAT.1");
	//cell.findGene("ENSMUSG00000037221");
	LinearRegression<short> linearRegressions(cell);
	linearRegressions.calculateLinearRegression(5);
	linearRegressions.writeFile(".data/123.csv");
	cin.get();
	cin.get();
	return 0;
}
