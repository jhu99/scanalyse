#include"cell/cell.cpp"
#include<iostream>
#include<cstdlib>
using namespace std;

int main(int argc, const char ** argv) {
	Cells<short> cells(atoi(argv[1]),atoi(argv[2]));
	cells.readFile(argv[3]);
	cells.maskCheck(argv[4]);
	return 0;
}
