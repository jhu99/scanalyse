#ifndef GENEINFO
#define GENEINFO
class geneInfo {
	int data;
	long long indices;
	double rank;
public:
	geneInfo();
	geneInfo(int data, long long indices);
	
	int getData();
	long long getIndices();
	double getRank();

	void setData(int data);
	void setIndices(long long indices);
	void setRank(double rank);
};
#endif // !GENEINFO
