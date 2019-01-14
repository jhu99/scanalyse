#include"Filtration/Filtration.h"
#include"argparser/argparser.h"
int main(int argc, const char ** argv)
{
	ArgParser a;
	const std::string gene_top_num_parameter = "gene_top_num";
	int gene_top_numVariable;
	int dft_gene_top_num = 500;
	const std::string help_gene_top_num = "gene_top_num  represent the number of top_genes you want to obtain.";
	a.refOption(gene_top_num_parameter, help_gene_top_num, gene_top_numVariable, dft_gene_top_num, true);
	
	const std::string min_cell_size_parameter = "min_cell_size";
	int min_cell_sizeVariable;
	int dft_min_cell_size = 1;
	const std::string help_min_cell_size = "min_cell_size represent the minimize gene count a cell need to have";
	a.refOption(min_cell_size_parameter, help_min_cell_size, min_cell_sizeVariable, dft_min_cell_size, true);

	const std::string method_type_parameter = "method_type";
	int method_typeVariable;
	int dft_method_type = 1;
	const std::string help_method_type = "method_type :method 1 is provided by JieLiu and 2 is provided by TianweiLiu.";
	a.refOption(method_type_parameter, help_method_type, method_typeVariable, dft_method_type, true);

	const std::string read_path_parameter = "read_path";
	string read_pathVariable;
	string dft_read_path = "";
	const std::string help_read_path = "read path represent the path of h5 file you want to read.";
	a.refOption(read_path_parameter, help_read_path, read_pathVariable, dft_read_path, true);

	const std::string write_path_parameter = "write_path";
	string write_pathVariable;
	string dft_write_path = "";
	const std::string help_write_path = "write_path represent the path of h5 file you want to write.";
	a.refOption(write_path_parameter, help_write_path, write_pathVariable, dft_write_path, true);

	a.run(argc,argv);
	cout << gene_top_numVariable << endl; 
 	SparseMatrix sm;
	string type = "original";
	sm.readHDF5File(read_pathVariable, type);
	Filtration f(sm);
	f.filtGeneAndCell(min_cell_sizeVariable, gene_top_numVariable, method_typeVariable);
	f.printFiltResult();
	f.writeFiltH5File(write_pathVariable);
	string deleteType= "original";
	sm.deleteSparseMatrix(deleteType);
	cin.get();
	cin.get();

}