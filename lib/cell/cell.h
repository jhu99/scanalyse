#ifndef CELL_H
#define CELL_H

#include<string>
#include<unordered_map>
#include<iostream>
#include<cstdlib>
#include<fstream>
#include<sstream>
#include<cstdlib>

using namespace std;
template <class T> class Cells {
	int row, col;
	string type;
	T **cell;
	unordered_map<string, int> cellToNum;
	unordered_map<string, int> geneToNum;
	unordered_map<int, string> numToCell;
	unordered_map<int, string> numToGene;
public:
	Cells();
	Cells(int row, int col,string type);
	~Cells();
	void setCell(T** cell);
	void setCellToNUm(unordered_map<string, int> cellToNum);
	void setGeneToNUm(unordered_map<string, int> geneToNum);
	void setNumTocell(unordered_map<int, string> numToCell);
	void setNumToGene(unordered_map<int, string> numToGene);

	T** getCell();
	unordered_map<string, int> getCellToNum();
	unordered_map<string, int> getGeneToNum();
	unordered_map<int, string> getNumToCell();
	unordered_map<int, string> getNumToGene();

	void readFile(string path);
	void findCell(string cellName);
	void findGene(string geneName);
	void findCellAndGene(string cellName, string geneName);
	void write2CSV(string write_path);
	void releaseMemory();

};

#endif // !CELL_H


