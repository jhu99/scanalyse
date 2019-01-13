#ifndef GENEVARIATION
#define GENEVARIATION
#include <string>
using namespace std;
class geneVariation {
	string geneName;
	long long index;
	double Variation=0.0;
public:
	void set_geneName(string geneName);
	void set_Variation(double Variation);
	void set_Index(long long index);

	string get_geneName();
	double get_Variation();
	long long get_index();

	void VarAdd(int num, double avg,int n);
	bool operator<(const geneVariation &a)const;
};
#endif // GENEVARIATION

