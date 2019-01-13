#ifndef GENECOUNT
#define GENECOUNT
#include<string>
using namespace std;
class geneCount {
	string geneName;
	long long index;
	int count=0;
public:
	void setGeneName(string geneName);
	void setCount(int count);
	void setIndex(long long index);

	string getGeneName(); 
	int getCount();
	long long getIndex();
	void countAdd(int num);
	bool operator < (const geneCount &x)const;

};
#endif // !GENECOUNT

