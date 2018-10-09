#include<string>
#include<unordered_map>
using namespace std;
#define N 5005


class Cells {
	unsigned short(*cell)[N];
	unordered_map<string, int> cellToNum;
	unordered_map<string, int> geneToNum;
	unordered_map<int, string> numToCell;
	unordered_map<int, string> numToGene;
public:
	Cells() {
		cell = new unsigned short[N][N];
	}
	~Cells() {
		delete cell;
	}
	void readFile();
	void findCell(string cellName);
	void findGene(string geneName);
	void findCellAndGene(string cellName, string geneName);
};
