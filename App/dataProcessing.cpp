#include"SparseMatrix/SparseMatrix.h"
#include"argparser/argparser.h"
#include"Filtration/Filtration.h"
#include"rank/rankNormalize.h"
#include"logNormalize/logNormalize.h"

struct Option
{
	string file_type;
	string read_path;
	string write_path;
	string normalize_type;
	int min_cell_size;
	int gene_top_num;
	int gene_filt_type;
	int thread_count;
};

int main(int argc, const char **argv)
{
	SparseMatrix sm;
	ArgParser a;
	Option option;
	a.refOption("file_type",
		"choose h5 or mtx",
		option.file_type,
		"h5", true);
	a.refOption("read_path",
		"path of the file you want to read(mtx_file only need folder path + /).",
		option.read_path,
		"", true);
	a.refOption("write_path",
		"path of the HDF5 file you want to write after filtration.",
		option.write_path,
		"", true);
	a.refOption("normalize_type",
		"two normalize method were provided: ",
		option.normalize_type,
		"", true);
	a.refOption("min_cell_size",
		"min_cell_size represent the minimize gene count a cell need to have",
		option.min_cell_size,
		1, true);
	a.refOption("gene_top_num",
		"gene_top_num  represent the number of top_genes you want to obtain.",
		option.gene_top_num,
		500, true);
	a.refOption("gene_filt_type",
		"gene_filt_type has two types: type 1 is provided by JieLiu and 2 is provided by TianweiLiu.",
		option.gene_filt_type,
		1, true);
	a.refOption("thread_count",
		"thread_count: the count of threads you want to create while get geneTop and logNormalize.",
		option.thread_count,
		1, true);
	a.run(argc, argv);
	cout << "start read-----------------" << endl;
	if (option.file_type.compare("h5") == 0)
	{
		sm.read_10x_h5(option.read_path);
	}
	else if (option.file_type.compare("mtx") == 0)
	{
		sm.read_10x_mtx(option.read_path);
	}
	cout << "end read-----------------" << endl;
	cout << "start filtration-----------------" << endl;
	Filtration f(sm);
	SparseMatrix filt_sm = f.filtGeneAndCell(option.min_cell_size, option.gene_top_num, option.gene_filt_type);
	f.printFiltResult();
	cout << "end filtration-----------------" << endl; 
	
	cout << "start normalize--------------" << endl;
	if (option.normalize_type == "rank")
	{	
		rankNormalize rn(filt_sm);
		rn.ranks(option.thread_count);
		rn.print();
		cout << "end normalize--------------" << endl;
		filt_sm.set_rank(rn.getRank());
		cout << "start write to h5 file-----------------" << endl;
		filt_sm.write_norm_data(option.write_path, 500, "s", "rank");
		cout << "end write to h5 file-----------------" << endl;
	}
	else if (option.normalize_type=="log")
	{
		logNormalize log(filt_sm);
		log.logNormalizeData(option.thread_count);
		cout << "end normalize--------------" << endl;
		filt_sm.set_log_data(log.get_log_data());
		cout << "start write to h5 file-----------------" << endl;
		filt_sm.write_norm_data(option.write_path, 500, "s", "log");
		cout << "end write to h5 file-----------------" << endl;
	}

	
	sm.deleteSparseMatrix("original");
	cin.get();
	cin.get();
}