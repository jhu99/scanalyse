import train
import sys
import load_model
#Argv 1: Pathway of the input file
#train.train_model(sys.argv[1],sys.argv[2],sys.argv[3])

#Argv 2: Pathway of the load_weight
load_model.load_weight(sys.argv[1],sys.argv[2],sys.argv[3])
