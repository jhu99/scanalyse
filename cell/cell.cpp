#include<iostream>
#include<cstdlib>
#include<string>
#include<fstream>
#include<sstream>
#include"cell.h"
using namespace std;
void Cells::readFile() {
	ifstream inFile("scimpute-0.01-1-data.csv", ios::in);
	string lineStr;
	string geneName;
	string str;
	int i = 0;
	if (getline(inFile, lineStr)) {
		stringstream ss(lineStr);
		while (getline(ss, str, ',')) {
			cellToNum[str] = i;
			numToCell[i] = str;
			i++;
		}
	}

	i = 1;
	while (getline(inFile, lineStr)) {
		stringstream ss(lineStr);
		int flag = 0, j = 1;
		while (getline(ss, str, ',')) {
			if (flag == 0) {
				geneToNum[str] = i;
				numToGene[i] = str;
				flag = 1;
			}
			else {
				cell[i][j] = atoi(str.c_str());
				j++;
			}
		}
		i++;
	}
}

void Cells::findCell(string cellName) {
	cellName = "\"" + cellName + "\"";
	int index = cellToNum[cellName];
	for (int i = 1; i <= geneToNum.size(); i++) {
		if (cell[i][index])
			cout << numToGene[i] << " " << cell[i][index] << endl;
	}
}

void Cells::findGene(string geneName) {
	geneName = "\"" + geneName + "\"";
	int index = geneToNum[geneName];
	for (int i = 1; i <= cellToNum.size(); i++) {
		if (cell[i][index])
			cout << numToCell[i] << " " << cell[i][index] << endl;
	}
}

void Cells::findCellAndGene(string cellName, string geneName) {
	cellName = "\"" + cellName + "\"";
	geneName = "\"" + geneName + "\"";
	cout << cell[geneToNum[geneName]][cellToNum[cellName]] << endl;
}
