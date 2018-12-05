#ifndef GENEINFO
#define GENEINFO
class geneInfo {
	int data;
	long long indices;
	unsigned short rank;
public:
	geneInfo();
	geneInfo(int data, long long indices);
	~geneInfo();
	
	int getData();
	long long getIndices();
	unsigned short getRank();

	void setData(int data);
	void setIndices(long long indices);
	void setRank(unsigned short rank);
};
#endif // !GENEINFO
