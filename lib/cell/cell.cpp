// CELL_CPP
// Created on Sep. 1, 2018
// Author: Tianwei Liu
#include "cell.h"

template <class T> Cells<T>::Cells() {

}

template <class T> Cells<T>::Cells(int n, int p) {
	this->n = n;
	this->p = p;
	cell = new T *[p];
	for (int i = 0; i < p; i++) {
		cell[i] = new T[n];
	}
}

template <class T> Cells<T>::~Cells() {
	/*for (int i = 0; i < n; i++) {
		delete[] cell[i];
	}
	delete[] cell;*/
}



template <class T> void Cells<T>::setCell(T** cell) {
	this->cell = cell;
}

template <class T> void Cells<T>::setCellToNUm(unordered_map<string, int> cellToNum) {
	this->cellToNum = cellToNum;
}

template <class T> void Cells<T>::setGeneToNUm(unordered_map<string, int> geneToNum) {
	this->geneToNum = geneToNum;
}

template <class T> void Cells<T>::setNumTocell(unordered_map<int, string> numToCell) {
	this->numToCell = numToCell;
}

template <class T> void Cells<T>::setNumToGene(unordered_map<int, string> numToGene) {
	this->numToCell = numToCell;
}

template <class T> T** Cells<T>::getCell() {
	return cell;
}

template <class T> unordered_map<string, int> Cells<T>::getCellToNum() {
	return cellToNum;
}

template <class T> unordered_map<string, int> Cells<T>::getGeneToNum() {
	return geneToNum;
}

template <class T> unordered_map<int, string> Cells<T>::getNumToCell() {
	return numToCell;
}

template <class T> unordered_map<int, string> Cells<T>::getNumToGene() {
	return numToGene;
}



template <class T> void Cells<T>::readFile(string path) {
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
			if (i > 0) {
				cellToNum[str] = i - 1;
				numToCell[i - 1] = str;
			}
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
				cell[j][i] = (T)atof(str.c_str());
				j++;
			}
		}
		i++;
	}
}

template <class T> void Cells<T>::findCell(string cellName) {
	int index = cellToNum[cellName];
	for (int i = 0; i < geneToNum.size(); i++) {
		if (cell[index][i])
			cout << numToGene[i] << " " << cell[index][i] << endl;
	}
}

template <class T> void Cells<T>::findGene(string geneName) {
	int index = geneToNum[geneName];
	for (int i = 0; i < cellToNum.size(); i++) {
		if (cell[i][index])
			cout << numToCell[i] << " " << cell[i][index] << endl;
	}
}

template <class T> void Cells<T>::findCellAndGene(string cellName, string geneName) {
	cout << cell[geneToNum[cellName]][cellToNum[geneName]] << endl;
}

template <class T> void Cells<T>::releaseMemory() {
	for (int i = 0; i < p; i++) {
		delete[] cell[i];

	}
	delete[] cell;
}
