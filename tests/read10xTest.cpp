#include"argparser/argparser.h"
#include"SparseMatrix/SparseMatrix.h"
using namespace std;
int main(int argc, const char ** argv)
{
	SparseMatrix sm;
	string file_type, read_path;
	ArgParser a;
	a.refOption("file_type",
		"choose h5 or mtx",
		file_type,
		"h5", true);
	a.refOption("read_path",
		"path of the file you want to read(mtx_file only need folder path + /).",
		read_path,
		"", true);
	a.run(argc, argv);
	if (file_type.compare("h5") == 0)
	{
		sm.read_10x_h5(read_path);
	}
	else if (file_type.compare("mtx") == 0)
	{
		sm.read_10x_mtx(read_path);
	}
	sm.deleteSparseMatrix("original");
	cin.get();
	cin.get();
}