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
	a.setName("appTest","This is a description of appTest.");
	a.setVersion("0.0.1");
	a.refOption("file_type",
		"Choose h5 or mtx. Default is h5.",
		option.file_type,
		"h5", false);
	a.refOption("read_path",
		"Path of the file you want to read (mtx_file only need folder path + /).",
		option.read_path,
		"", true);
	a.refOption("write_path",
		"Path of the HDF5 file you want to write after filtration.",
		option.write_path,
		"", true);
	a.refOption("normalize_type",
		"Two normalize methods are provided: log (log normalization) and rank (rank normalization). Default is log.",
		option.normalize_type,
		"log", false);
	a.refOption("min_cell_size",
		"It filter these cell which have less read counts than min_cell_size. Default is 1.",
		option.min_cell_size,
		1, false);
	a.refOption("gene_top_num",
		"It will select the top_num genes ranked by read counts or variance. Default is 1000.",
		option.gene_top_num,
		1000, false);
	a.refOption("gene_filt_type",
		"The parameter gene_filt_type has two options: 1 and 2. Type 1 is ranked by variance and type 2 is ranked by gene expression.",
		option.gene_filt_type,
		1, false);
	a.refOption("thread_num",
		"The parameter thread_num is the number of threads you want to use for the parallelization.",
		option.thread_count,
		1, true);
	if(!a.run(argc, argv))return 0;
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
	cout << "end write to h5 file-----------------" << endl;
	cout << "start normalize--------------" << endl;
	if (option.normalize_type == "rank")
	{	
		rankNormalize rn(filt_sm);
		rn.ranks(option.thread_count);
	
		cout << "end normalize--------------" << endl;
		filt_sm.set_rank(rn.getRank());
		cout << "start write to h5 file-----------------" << endl;
		filt_sm.write_norm_data(option.write_path, "rank", 500, "s");
		cout << "end write to h5 file-----------------" << endl;
	}
	else if (option.normalize_type=="log")
	{
		logNormalize log(filt_sm);
		log.logNormalizeData(option.thread_count);
		cout << "end normalize--------------" << endl;
		filt_sm.set_log_data(log.get_log_data());
		cout << "start write to h5 file-----------------" << endl;
		filt_sm.write_norm_data(option.write_path, "log", 500, "s");
		cout << "end write to h5 file-----------------" << endl;
	}
	sm.deleteSparseMatrix("original");
	SparseMatrix sm2;
	cin.get();
	cin.get();
}
