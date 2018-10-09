#ifndef CELL_H
#define CELL_H

#include<string>
#include<unordered_map>
using namespace std;


class Cells {
	unsigned short **cell;
	unordered_map<string, int> cellToNum;
	unordered_map<string, int> geneToNum;
	unordered_map<int, string> numToCell;
	unordered_map<int, string> numToGene;
public:
	Cells() {

	}
	Cells(int N,int P) {
		cell = new unsigned short *[N];
		for (int i = 0; i <= N; i++) {
			cell[i] = new unsigned short[P];
		}
	}
	~Cells() {
		delete cell;
	}
	void readFile();
	void findCell(string cellName);
	void findGene(string geneName);
	void findCellAndGene(string cellName, string geneName);
};

#endif // !CELL_H
