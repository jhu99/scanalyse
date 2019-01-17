import train
import sys
import load_model

FUN = sys.argv[1]
INPUT_FILE=sys.argv[2]
OUTPUT_FILE=sys.argv[3]
FMT=sys.argv[4]

#Argv 1: Pathway of the input file
if FUN=='train':
	train.train_model(INPUT_FILE,OUTPUT_FILE,FMT)
#Argv 2: Pathway of the load_weight
elif FUN=='prediction':
	load_model.load_weight(sys.argv[1],sys.argv[2],sys.argv[3])
else:
	print "please use \'train\' or \'prediction\'"
