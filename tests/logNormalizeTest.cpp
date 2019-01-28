#include"logNormalize/logNormalize.h"
#include"SparseMatrix/SparseMatrix.h"
int main(int argc, const char ** argv)
{
	SparseMatrix sm;
	sm.read_10x_h5(argv[1]);
	logNormalize log(sm);
	log.logNormalizeData(5);
	cin.get();
	cin.get();
	return 0;
}