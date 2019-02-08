import train
import load_model
import argparse
import time

t0=time.clock()
parser = argparse.ArgumentParser(description = "This tool is used to analyse single cell data.")
parser.add_argument('--task','-t', required=True,
					choices=['train','prediction'],
					help='The task will be performed. It has two options, \'train\' or \'prediction\'. Default is \'train\'.')
parser.add_argument('--input_file','-i', required = True, 
					help='The path of an input filename.')
parser.add_argument('--gene_file','-g', required= True,
					help='The path of a gene filename.', default='input_tensor.csv')
parser.add_argument('--weight_file','-w', 
					help='The path of a weight filename.')
parser.add_argument('--output_path','-o', required = True,
					help='The path of the output directory.')
parser.add_argument('--latent_size','-l', type=int, nargs='+', default=[64,32,64],
					help='The size of latent layers.')
parser.add_argument('--batch_size','-b', type=int, default=128,
					help='The size of each batch.')
parser.add_argument('--algorithm','-a', default='rmsprop', 
					choices=['adam','rmsprop'],
					help='The optimization algorithms. It has two options, \'adam\' or \'rmsprop\'. Default is \'rmsprop\'.')
parser.add_argument('--format','-f', required= True, 
					choices=['10x_h5','10x_mtx'],
					help='The format of input files. It has two options, \'10x_h5\' or \'10x_mtx\'.')
args = parser.parse_args()

#Argv 1: Pathway of the input file
if args.task =='train':
	train.train_model(input_file=args.input_file, 
					weight_file=args.weight_file, 
					output_path=args.output_path,
					gene_file=args.gene_file,
					hidden_size=args.latent_size,
					batch_size=args.batch_size, 
					optimizer=args.algorithm,
					format_type=args.format)
#Argv 2: Pathway of the load_weight
elif args.task =='prediction':
	load_model.load_weight(input_file=args.input_file,
						weight_file=args.weight_file,
						gene_file=args.gene_file,
						hidden_size=args.latent_size,
						output_path=args.output_path,
						format_type=args.format)
else:
	print('Please use \'train\' or \'prediction\' to specific --task.')
	
t1=time.clock()
print('Time elapsed: ', t1-t0)
