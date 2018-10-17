#include"fun.h"
#include"cell.h"
#include <io.h>
#include<iostream>
#include <fstream>  
#include <string>  
#include <vector>  
#include<queue>
#include<sstream>
#include"cell.cpp"
#define WINDOS
//#define LINUX 
#ifdef LINUX
#include <cstring>
#include<dirent.h>
#endif // LINUX

typedef float dataType;

using namespace std;
using namespace Scanalyse;




void Fun::read(string filePath)
{

	//read all csv file  
	string format = ".csv";
	Fun fun;
	fun = Fun();
	fun.GetAllFiles(filePath, files);

	int size = files.size();
	for (int i = 0; i<size; i++)
	{
		cout << files[i] << endl;
		rowCount[i] = CaculateRow(files[i]);
		columnCount[i] = CaculateColumn(files[i]);
		mergeColumnCount += columnCount[i];
		cout << mergeColumnCount << endl;
	}
	
	int j = 0;
	while(j<size)
	{
		cout << "end read";
		Cells<dataType> paraCells(rowCount[j], columnCount[j]);
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
#ifdef WINDOS
	long   hFile = 0;
	//file information    
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
#endif // WINDOS
#ifdef LINUX
	DIR *dir;
	struct dirent *ptr;
	char base[1000];
	string paraPath;
	if ((dir = opendir(path.c_str())) == NULL)
	{
		perror("Open dir error..."); exit(1);
	}
	while ((ptr = readdir(dir)) != NULL)
	{
		if (strcmp(ptr->d_name, ".") == 0 || strcmp(ptr->d_name, "..") == 0)
			///current dir OR parrent dir
		{
			continue;
		}
		else if (ptr->d_type == 8) ///file 
		{
			//printf("d_name:%s/%s\n",basePath,ptr->d_name);
			files.push_back(paraPath.assign(path).append("/").append(ptr->d_name));
		}
		else if (ptr->d_type == 10) ///link file 
		{
			//printf("d_name:%s/%s\n",basePath,ptr->d_name); 
			continue;
		}
		else if (ptr->d_type == 4) ///dir 
		{
			//files.push_back(ptr->d_name);
			/* memset(base,'\0',sizeof(base));
			strcpy(base,basePath);
			strcat(base,"/");
			strcat(base,ptr->d_nSame);
			readFileList(base);
			*/
			continue;
		}
	}
	closedir(dir);
#endif // LINUX

	

}



//obtain specific format file
void Fun::GetAllFormatFiles(string path, vector<string>& files, string format)

{ 
	/*long   hFile = 0;
	//file information    
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

					GetAllFormatFiles(p.assign(path).append("\\").append(fileinfo.name), files, format);
				}

			}
			else
			{
				files.push_back(p.assign(path).append("\\").append(fileinfo.name));

			}
		} while (_findnext(hFile, &fileinfo) == 0);

		_findclose(hFile);

	}*/
}

void Fun::CreatAllGeneMap()
{
	unordered_map<string, int> paraGeneToNum;
	unordered_map<int, string> paraNumToGene;
	int i = 0, iter;
	Cells<dataType> paraCell;
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

void Fun::replaceEnsemblId(string path1,string path2)
{
	
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

	string distAll = "E:\\cell mapping\\GSE66507.txt";
	ofstream ofn(distAll);
	int size = CaculateRow(path1);
	for (i = 0; i<size; i++)
	{
		ofn << numToGene[i] << endl;
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
			mergeMatrix[i][j] = 0;
		}
	}
	cout << "end init" << endl;
}

void Scanalyse::Fun::outToCsvFile(string outPath)
{
	ofstream ofn(outPath);
	for (int k = 0; k < allNumToCell.size(); k++)
	{
		ofn << allNumToCell[k];
		ofn << ",";
	}
	ofn << "\n";
	for (int i = 0; i < allGeneToNum.size(); i++)
	{
		ofn << allNumToGene[i];
		ofn << ",";
		for (int j = 0; j < mergeColumnCount; j++)
		{
			if (j != mergeColumnCount - 1)
			{
				ofn << mergeMatrix[i][j];
					ofn<<"," ;
			}
			
		}
		ofn << "\n";
	}
	ofn.close();
	cout << "end out file" << endl;


}


void Fun::mergeMatrixs()
{
	Cells<dataType> paraCell;
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
			allNumToCell[cellCount] = paraNumToCell[i];
			allCellToNum[paraNumToCell[i]] = cellCount;
			for (int j = 0; j < paraNumToGene.size(); j++)
			{
				genePosition = allGeneToNum[paraNumToGene[j]];
				mergeMatrix[genePosition][cellCount] = paraCellMatrix[j][i];
				
			}
			cellCount++;
		}
		cout << "finished No." << iter << endl;
	}
	/*for (int k = 0; k < allGeneToNum.size(); k++)
	{
		for (int t = 0; t < mergeColumnCount; t++)
		{
			cout << mergeMatrix[k][t]<<" ";
		}
		cout << endl;
	}
	*/
	cout << allGeneToNum.size() << endl;
	cout << mergeColumnCount << endl;
	cout << "end merge" << endl;
	for (int iter = 0; iter < cellList.size(); ++iter)
	{
		cellList[iter].releaseMemory();
	}
}
