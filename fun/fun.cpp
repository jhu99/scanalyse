#include"fun.h"
#include"cell.h"
#include <io.h>
#include<iostream>
#include <fstream>  
#include <string>  
#include <vector>  
#include<queue>
#include<sstream>

typedef unsigned short dataType;

using namespace std;
using namespace Scanalyse;




void Fun::read()
{
	string filePath = "E:\cell";
	string distAll = "AllFiles.txt";
	//读取所有的文件，包括子文件的文件 
	//GetAllFiles(filePath, files);  

	//读取所有格式为csv的文件  
	string format = ".csv";
	Fun fun;
	fun = Fun();
	fun.GetAllFormatFiles(filePath, files, format);
	//ofstream ofn(distAll);

	int size = files.size();
	//ofn << size << endl;
	for (int i = 0; i<size; i++)
	{
		//ofn << files[i] << endl;
		cout << files[i] << endl;
		rowCount[i] = CaculateRow(files[i]);
		columnCount[i] = CaculateColumn(files[i]);
		mergeColumnCount += columnCount[i];
	}
	//ofn.close();
	int j = 0;
	while(j<size)
	{
		cout << "end read";
		Cells paraCells(rowCount[j], columnCount[j]);
		paraCells.readFile(files[j]);
		cellList.push_back(paraCells);
		j++;
	}
	
}

int Fun::CaculateRow(string path)
{
	ifstream inFile(path, ios::in);
	string lineStr;
	int row = 0;

	while (getline(inFile, lineStr))
	{
		row++;
	}
	return row;
}

int Fun::CaculateColumn(string path)
{
	ifstream inFile(path, ios::in);
	string lineStr;
	string geneName;
	string str;
	int column = 0;

	if (getline(inFile, lineStr)) {
		stringstream ss(lineStr);
		while (getline(ss, str, ','))
		{
			column++;
		}
	}
	return column;
}


void Fun::GetAllFiles(string path, vector<string>& files)
{
	long   hFile = 0;
	//文件信息    
	struct _finddata_t fileinfo;
	string p;
	if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1)
	{
		do
		{
			if ((fileinfo.attrib &  _A_SUBDIR))
			{
				if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
				{
					files.push_back(p.assign(path).append("\\").append(fileinfo.name));
					GetAllFiles(p.assign(path).append("\\").append(fileinfo.name), files);
				}
			}
			else
			{
				files.push_back(p.assign(path).append("\\").append(fileinfo.name));
			}
		} while (_findnext(hFile, &fileinfo) == 0);
		_findclose(hFile);
	}

}



//获取特定格式的文件名  
void Fun::GetAllFormatFiles(string path, vector<string>& files, string format)

{
	//文件句柄    
	long   hFile = 0;
	//文件信息    
	struct _finddata_t fileinfo;
	string p;
	if ((hFile = _findfirst(p.assign(path).append("\\*" + format).c_str(), &fileinfo)) != -1)
	{
		do
		{
			if ((fileinfo.attrib &  _A_SUBDIR))
			{
				if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
				{
					//files.push_back(p.assign(path).append("\\").append(fileinfo.name) );  

					GetAllFormatFiles(p.assign(path).append("\\").append(fileinfo.name), files, format);
				}

			}
			else
			{
				files.push_back(p.assign(path).append("\\").append(fileinfo.name));

			}
		} while (_findnext(hFile, &fileinfo) == 0);

		_findclose(hFile);

	}
}

void Fun::CreatAllGeneMap()
{
	unordered_map<string, int> paraGeneToNum;
	unordered_map<int, string> paraNumToGene;
	int i = 0, iter;
	Cells paraCell;
	string paraGene;
	for (iter = 0; iter <cellList.size(); ++iter) {
		paraCell = cellList[iter];
		paraGeneToNum = paraCell.getGeneToNum();
		paraNumToGene = paraCell.getNumToGene();
		for (int j = 0; j < paraGeneToNum.size(); j++)
		{
			paraGene= paraNumToGene[j];
			if (allGeneToNum.find(paraGene)==allGeneToNum.end())
			{
				allNumToGene[i] = paraGene;
				allGeneToNum[paraGene] = i;
				i++;
			}

		}
	}
}

void Fun::replaceEnsemblId()
{
	string path1 = "E:\\cell mapping\\GSE52529.csv";
	string path2 = "E:\\cell mapping\\GSE52529_mapping.csv";
	unordered_map<string, int> geneToNum;
	unordered_map<int, string> numToGene;
	string lineStr;
	string ensemblId;
	string str,str1,str2;
	int i = 0;

	ifstream inFile(path1, ios::in);
	
	while (getline(inFile, lineStr)) {
		stringstream ss(lineStr);
		while (getline(ss, str, ',')) {
			
				
				geneToNum[str] = i;
				numToGene[i] = str;
				i++;
		}
		
	}
	ifstream inFile2(path2, ios::in);
	while (getline(inFile2, lineStr)) {
		stringstream ss(lineStr);
		int  flag = 0;
		while (getline(ss, str, ',')) {
			if (flag == 0)
			{
				str1 = str;
				flag = 1;
			}
			else
			{
				str2= str;
				geneSymbolToEnsemblId[str1] = str2;
			}
		}
	}
	string paraEnsemblId;
	for (int k = 0; k < geneToNum.size(); k++)
	{
		paraEnsemblId = numToGene[k];
		
		if (geneSymbolToEnsemblId.find(paraEnsemblId) != geneSymbolToEnsemblId.end())
		{
			numToGene[k] = geneSymbolToEnsemblId.find(paraEnsemblId)->second;
			geneToNum[geneSymbolToEnsemblId.find(paraEnsemblId)->second] = k;
		}
		else
		{
			numToGene[k] = "NaN";
		}
	}

	string distAll = "E:\\cell mapping\\GSE52529.txt";
	ofstream ofn(distAll);

	for (i = 0; i<47193; i++)
	{
		ofn << numToGene[i] << endl;
		cout << numToGene[i] << endl;
	}
	ofn.close();
}

void Fun::initMergeMatrix()
{
	mergeMatrix = new dataType *[allGeneToNum.size()];
	for (int i = 0; i < allGeneToNum.size(); i++) {
		mergeMatrix[i] = new dataType[mergeColumnCount];
	}
	
	for (int i = 0; i < allGeneToNum.size(); i++)
	{
		for (int j = 0; j < mergeColumnCount; j++)
		{
			mergeMatrix[i][j] = NAN;
		}
	}
	cout << "end init" << endl;
}


void Fun::mergeMatrixs()
{
	Cells paraCell;
	unordered_map<string, int> paraGeneToNum;
	unordered_map<int, string> paraNumToGene;
	unordered_map<string, int> paraCellToNum;
	unordered_map<int, string> paraNumToCell;
	int genePosition;
	int cellCount = 0;
	dataType** paraCellMatrix;
	for (int iter = 0; iter < cellList.size(); ++iter) {
		paraCell = cellList[iter];
		paraNumToCell = paraCell.getNumToCell();
		paraNumToGene = paraCell.getNumToGene();
		paraCellMatrix = paraCell.getCell();
		for (int i = 0; i < paraNumToCell.size(); i++)
		{
			
			for (int j = 0; j < paraNumToGene.size(); j++)
			{
				genePosition = allGeneToNum[paraNumToGene[j]];
				mergeMatrix[genePosition][cellCount] = paraCellMatrix[j][i];
				
			}
			cellCount++;
		}
	}
	for (int k = 0; k < allGeneToNum.size(); k++)
	{
		for (int t = 0; t < mergeColumnCount; t++)
		{
			cout << mergeMatrix[k][t]<<" ";
		}
		cout << endl;
	}
	cout << "end merge" << endl;
	for (int iter = 0; iter < cellList.size(); ++iter)
	{
		cellList[iter].releaseMemory();
	}
}
