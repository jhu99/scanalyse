#include"Filtration/Filtration.h"
#include"argparser/argparser.h"
int main(int argc, const char ** argv)
{
	ArgParser a;
	int gene_top_numVariable,min_cell_sizeVariable,method_typeVariable;
	a.refOption("gene_top_num", 
				"gene_top_num  represent the number of top_genes you want to obtain.", 
				gene_top_numVariable, 
				1000,
				true);
	a.refOption("min_cell_size", 
				"min_cell_size represent the minimize gene count a cell need to have", 
				min_cell_sizeVariable, 
				1, 
				true);
	a.refOption("method_type", 
				"method_type :method 1 is provided by JieLiu and 2 is provided by TianweiLiu.", 
				method_typeVariable, 
				1, 
				true);

	string read_pathVariable,write_pathVariable;
	a.refOption("read_path", 
				"read path represent the path of h5 file you want to read.", 
				read_pathVariable, 
				"", 
				true);
	a.refOption("write_path", 
				"write_path represent the path of h5 file you want to write.", 
				write_pathVariable, 
				"",
				true);

	a.run(argc,argv);
	cout << gene_top_numVariable << endl; 
 	SparseMatrix sm;
	sm.readHDF5File(read_pathVariable, "original");
	Filtration f(sm);
	f.filtGeneAndCell(min_cell_sizeVariable, gene_top_numVariable, method_typeVariable);
	f.printFiltResult();
	f.writeFiltH5File(write_pathVariable);
	sm.deleteSparseMatrix("original");
	cin.get();
	cin.get();

}
