#ifndef CELL_H
#define CELL_H

#include<string>
#include<unordered_map>
using namespace std;

class Cells {
	int n, p;
	unsigned short **cell;
	unordered_map<string, int> cellToNum;
	unordered_map<string, int> geneToNum;
	unordered_map<int, string> numToCell;
	unordered_map<int, string> numToGene;
public:
	Cells() {

	}
	Cells(int nn,int pp):n(nn),p(pp){
		cell = new unsigned short *[n];
		for (int i = 0; i <= n; i++) {
			cell[i] = new unsigned short[p];
		}
	}
	~Cells() {
		for (int i = 0; i <= n; i++) {
			delete[] cell[i];
		}
		delete[] cell;
	}
	void readFile(string path);
	void findCell(string cellName);
	void findGene(string geneName);
	void findCellAndGene(string cellName, string geneName);
};

#endif // !CELL_H
