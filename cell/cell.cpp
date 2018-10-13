#include<iostream>
#include<cstdlib>
#include<string>
#include<fstream>
#include<sstream>
#include"cell.h"
typedef unsigned short dataType;
using namespace std;

Cells::Cells() {

}

Cells::Cells(int n, int p) {
	this->n = n;
	this->p = p;
	cell = new dataType *[n];
	for (int i = 0; i < n; i++) {
		cell[i] = new dataType[p];
	}
}

Cells::~Cells() {
	/*for (int i = 0; i < n; i++) {
		delete[] cell[i];
		cout << i;
	}
	delete[] cell;
	*/
}

void Cells::setCell(dataType ** cell)
{
	this->cell = cell;
}

void Cells::setCellToNUm(unordered_map<string, int> cellToNum) {
	this->cellToNum = cellToNum;
}

void Cells::setGeneToNUm(unordered_map<string, int> geneToNum) {
	this->geneToNum = geneToNum;
}

void Cells::setNumTocell(unordered_map<int, string> numToCell) {
	this->numToCell = numToCell;
}

void Cells::setNumToGene(unordered_map<int, string> numToGene) {
	this->numToCell = numToCell;
}

dataType** Cells::getCell() {
	return cell;
}

unordered_map<string, int> Cells::getCellToNum() {
	return cellToNum;
}

unordered_map<string, int> Cells::getGeneToNum() {
	return geneToNum;
}

unordered_map<int, string> Cells::getNumToCell() {
	return numToCell;
}

unordered_map<int, string> Cells::getNumToGene() {
	return numToGene;
}

void Cells::readFile(string path) {
	ifstream inFile(path, ios::in);
	string lineStr;
	string geneName;
	string str;
	string suffixName;
	char separator;

	suffixName = path.substr(path.size() - 3, path.size() - 1);
	if (suffixName == "csv") {
		separator = ',';
	}
	else {
		separator = '\t';
	}
	int i = 0;
	if (getline(inFile, lineStr)) {
		stringstream ss(lineStr);
		while (getline(ss, str, separator)) {
			//str = str.substr(1, str.size() - 2);
			cellToNum[str] = i;
			numToCell[i] = str;
			i++;
		}
	}

	i = 0;
	while (getline(inFile, lineStr)) {
		stringstream ss(lineStr);
		int flag = 0, j = 0;
		while (getline(ss, str, separator)) {
			if (flag == 0) {
				//str = str.substr(1, str.size() - 2);
				geneToNum[str] = i;
				numToGene[i] = str;
				flag = 1;
			}
			else {
				cell[i][j] = (dataType)atoi(str.c_str());
				j++;

			}
		}
		i++;
	}
	cout << "end read" << endl;
}

void Cells::findCell(string cellName) {
	int index = cellToNum[cellName];
	int size = numToGene.size();
	dataType* cellArray = new dataType[numToGene.size()];
	for (int i = 1; i <= geneToNum.size(); i++) {
		if (cell[i][index])
			cout << numToGene[i] << " " << cell[i][index] << endl;
	}
}

void Cells::findGene(string geneName) {
	int index = geneToNum[geneName];
	for (int i = 1; i <= cellToNum.size(); i++) {
		if (cell[i][index])
			cout << numToCell[i] << " " << cell[i][index] << endl;
	}
}

void Cells::findCellAndGene(string cellName, string geneName) {
	cout << cell[geneToNum[geneName]][cellToNum[cellName]] << endl;
}

void Cells::releaseMemory()
{
	for (int i = 1; i <n; i++) {
		delete[] cell[i];
		
	}
	delete[] cell;
}
