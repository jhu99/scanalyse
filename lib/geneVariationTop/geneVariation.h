#ifndef GENEVARIATION
#define GENEVARIATION
#include <string>
using namespace std;
class geneVariation {
	string geneName;
	double Variation=0.0;
public:
	void set_geneName(string geneName);
	void set_Variation(double Variation);

	string get_geneName();
	double get_Variation();
	void VarAdd(int num, double avg,int n);
	bool operator<(const geneVariation &a)const;
};
#endif // GENEVARIATION

