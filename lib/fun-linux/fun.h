#pragma once
#include<string>
#include<vector>
#include<unordered_map>
#include <cstring>
#include"cell.h"
using namespace std;
typedef float dataType;


namespace Scanalyse
{
	class Fun
	{
	public:
		Fun() {

		}
		~Fun() {

		}
		unordered_map<int, string> allNumToGene;
		unordered_map<string, int> allGeneToNum;
		unordered_map<int, string> allNumToCell;
		unordered_map<string, int> allCellToNum;
		vector<string> files;//存放文件夹中所有文件路径
		int rowCount[100];//每个文件的行数
		int columnCount[100];//每个文件的列数
		vector<Cells> cellList;//存放每个文件的cells实例
		unordered_map<string, string> geneSymbolToEnsemblId;
		int mergeColumnCount=0;
		dataType **mergeMatrix;

		void read(string filePath);
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

