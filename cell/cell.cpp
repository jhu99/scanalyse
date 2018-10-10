#include<iostream>
#include<cstdlib>
#include<string>
#include<fstream>
#include<sstream>
#include"cell.h"
using namespace std;

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
			str = str.substr(1, str.size() - 2);
			cellToNum[str] = i;
			numToCell[i] = str;
			i++;
		}
	}

	i = 1;
	while (getline(inFile, lineStr)) {
		stringstream ss(lineStr);
		int flag = 0, j = 1;
		while (getline(ss, str, separator)) {
			if (flag == 0) {
				str = str.substr(1, str.size() - 2);
				geneToNum[str] = i;
				numToGene[i] = str;
				flag = 1;
			}
			else {
				cell[i][j] =(unsigned short) atoi(str.c_str());
				j++;
			}
		}
		i++;
	}
}

void Cells::findCell(string cellName) {
	int index = cellToNum[cellName];
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
