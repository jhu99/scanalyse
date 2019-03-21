#pragma once
#include<string>
#include<unordered_map>
#include <vector>  
#include<iostream>
#include <fstream>  
#include<sstream>
#include<limits>

//#define WINDOWS
#define LINUX 
#ifdef LINUX
#include <cstring>
#include<dirent.h>
#include"cell/cell.cpp"
#include"cell/cell.h"
#endif // LINUX

#ifdef WINDOWS
#include<io.h>
#include"cell.cpp"
#include"cell.h"
#endif // WINDOWS

using namespace std;
typedef int dataType;

namespace Scanalyse
{
	class Fun
	{
	public:
		Fun() {

		}
		~Fun() {

		}
		unordered_map<long, string> allNumToGene;
		unordered_map<string, long> allGeneToNum;
		unordered_map<long, string> allNumToCell;
		unordered_map<string, long> allCellToNum;
		vector<string> files;//store path of each file
		int rowCount[100];//store row count of each file
		int columnCount[100];
		vector<Cells<dataType>> cellList;
		unordered_map<string, string> geneSymbolToEnsemblId;
		long mergeColumnCount=0;
		dataType **mergeMatrix;

		void read(string filePath);
		void transfer_matrix(string file_path,string write_path,string type);
		int CaculateRow(string path);
		int CaculateColumn(string path);
		void GetAllFiles(string path, vector<string>& files);
		void GetAllFormatFiles(string path, vector<string>& files, string format);
		void CreatAllGeneMap();
		void replaceEnsemblId(string path1, string path2);
		void mergeMatrixs();
		void initMergeMatrix();
		void outToCsvFile(string outPath);
	};
}

