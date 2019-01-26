#!/bin/bash
#PROGRAM = python3.6 ./main.py
#293t_in = ../data/cf.10xgenomics.com/samples/cell-exp/1.1.0/293t/filtered_matrices_mex/hg19/

#$PROGRAM train ../data/ica_bone_marrow_h5.h5 ../result/ica_bone_marrow 10x_h5
#$PROGRAM train ../data/ica_cord_blood_h5.h5 ../result/ica_cord_blood 10x_h5
#$PROGRAM train ../data/cf.10xgenomics.com/samples/cell-exp/1.1.0/fresh_68k_pbmc_donor_a/filtered_matrices_mex/hg19/ ../result/fresh_68k_pbmc_donor_a 10x_mtx
#$PROGRAM train ../data/cf.10xgenomics.com/samples/cell-exp/1.1.0/293t/filtered_matrices_mex/hg19/ ../result/1.1.0/293t 10x_mtx
##python3.6 -m pdb ./main.py train ../data/ica_bone_marrow_h5.h5 ../result/ica_bone_marrow 10x_h5
##python3 -m pdb ./main.py train ../data/ica_cord_blood_h5.h5 ../result/ica_cord_blood 10x_h5
##python3 -m pdb ./main.py train ../data/cf.10xgenomics.com/samples/cell-exp/1.1.0/fresh_68k_pbmc_donor_a/filtered_matrices_mex/hg19/ ../result/1.1.0/fresh_68k_pbmc_donor_a 10x_mtx
#filepath = '%(output_dir) has %(weight_file)' % vars()

##run on 293t
#wget -r http://cf.10xgenomics.com/samples/cell-exp/1.1.0/293t/293t_filtered_gene_bc_matrices.tar.gz
#python3.6 -m pdb ./main.py -t train -i ../data/cf.10xgenomics.com/samples/cell-exp/1.1.0/293t/filtered_matrices_mex/hg19/ -o ../result/1.1.0/293t -b 32 -l 64 32 64 -a adam -w weight_adam_64_32_64_g1_c1_b32.h5 -f 10x_mtx
#python3 ./main.py -t train -i ../data/cf.10xgenomics.com/samples/cell-exp/1.1.0/293t/filtered_matrices_mex/hg19/ -o ../result/1.1.0/293t -b 32 -l 64 32 64 -a adam -w weight_adam_64_32_64_g1_c1_b32.h5 -f 10x_mtx

##run on ica_all
#python3 ./main.py -t train -i xxx -o xxx -w xxx -f 10x_mtx -a adam
#mkdir ../result/ica_all run on P710
a="rmsprop"
b=1024
for l in "64 32 64" "512 64 512" "1024 128 1024"
do
	px=''
	for w in $l
	do
		px+='_'$w
	done
	weight_file="weight_${a}${px}_g1_c1_b$b.h5"
	python3.6 ./main.py -t train -i ../data/ica_all.h5 -o ../result/ica_all -b $b -l $l -a $a -w ${weight_file} -f 10x_h5 > ../result/ica_all/$a${px}_b$b.log
done

#python3.6 ./main.py -t train -i ../data/ica_all.h5 -o ../result/ica_all -b 32 -l 64 32 64 -a adam -w weight_adam_64_32_64_g1_c1_b32.h5 -f 10x_h5
#python3.6 ./main.py -t train -i ../data/ica_all.h5 -o ../result/ica_all -b 128 -l 64 32 64 -a adam -w weight_adam_64_32_64_g1_c1_b32.h5 -f 10x_h5
#python3.6 ./main.py -t train -i ../data/ica_all.h5 -o ../result/ica_all -b 512 -l 64 32 64 -a adam -w weight_adam_64_32_64_g1_c1_b32.h5 -f 10x_h5
#python3.6 ./main.py -t train -i ../data/ica_all.h5 -o ../result/ica_all -b 1024 -l 64 32 64 -a adam -w weight_adam_64_32_64_g1_c1_b32.h5 -f 10x_h5
#python3.6 ./main.py -t train -i ../data/ica_all.h5 -o ../result/ica_all -b 32 -l 64 32 64 -a rmsprop -w weight_adam_64_32_64_g1_c1_b32.h5 -f 10x_h5
#python3.6 ./main.py -t train -i ../data/ica_all.h5 -o ../result/ica_all -b 128 -l 64 32 64 -a rmsprop -w weight_adam_64_32_64_g1_c1_b32.h5 -f 10x_h5
#python3.6 ./main.py -t train -i ../data/ica_all.h5 -o ../result/ica_all -b 512 -l 64 32 64 -a rmsprop -w weight_adam_64_32_64_g1_c1_b32.h5 -f 10x_h5
#python3.6 ./main.py -t train -i ../data/ica_all.h5 -o ../result/ica_all -b 1024 -l 64 32 64 -a rmsprop -w weight_adam_64_32_64_g1_c1_b32.h5 -f 10x_h5
#python3.6 ./main.py -t train -i ../data/ica_all.h5 -o ../result/ica_all -b 32 -l 512 64 512 -a adam -w weight_adam_64_32_64_g1_c1_b32.h5 -f 10x_h5
#python3.6 ./main.py -t train -i ../data/ica_all.h5 -o ../result/ica_all -b 128 -l 512 64 512 -a adam -w weight_adam_64_32_64_g1_c1_b32.h5 -f 10x_h5
#python3.6 ./main.py -t train -i ../data/ica_all.h5 -o ../result/ica_all -b 512 -l 512 64 512 -a adam -w weight_adam_64_32_64_g1_c1_b32.h5 -f 10x_h5
#python3.6 ./main.py -t train -i ../data/ica_all.h5 -o ../result/ica_all -b 1024 -l 512 64 512 -a adam -w weight_adam_64_32_64_g1_c1_b32.h5 -f 10x_h5
#python3.6 ./main.py -t train -i ../data/ica_all.h5 -o ../result/ica_all -b 32 -l 512 64 512 -a rmsprop -w weight_adam_64_32_64_g1_c1_b32.h5 -f 10x_h5
#python3.6 ./main.py -t train -i ../data/ica_all.h5 -o ../result/ica_all -b 128 -l 512 64 512 -a rmsprop -w weight_adam_64_32_64_g1_c1_b32.h5 -f 10x_h5
#python3.6 ./main.py -t train -i ../data/ica_all.h5 -o ../result/ica_all -b 512 -l 512 64 512 -a rmsprop -w weight_adam_64_32_64_g1_c1_b32.h5 -f 10x_h5
#python3.6 ./main.py -t train -i ../data/ica_all.h5 -o ../result/ica_all -b 1024 -l 512 64 512 -a rmsprop -w weight_adam_64_32_64_g1_c1_b32.h5 -f 10x_h5
#python3.6 ./main.py -t train -i ../data/ica_all.h5 -o ../result/ica_all -b 32 -l 1024 128 1024 -a adam -w weight_adam_64_32_64_g1_c1_b32.h5 -f 10x_h5
#python3.6 ./main.py -t train -i ../data/ica_all.h5 -o ../result/ica_all -b 128 -l 1024 128 1024 -a adam -w weight_adam_64_32_64_g1_c1_b32.h5 -f 10x_h5
#python3.6 ./main.py -t train -i ../data/ica_all.h5 -o ../result/ica_all -b 512 -l 1024 128 1024 -a adam -w weight_adam_64_32_64_g1_c1_b32.h5 -f 10x_h5
#python3.6 ./main.py -t train -i ../data/ica_all.h5 -o ../result/ica_all -b 1024 -l 1024 128 1024 -a adam -w weight_adam_64_32_64_g1_c1_b32.h5 -f 10x_h5
#python3.6 ./main.py -t train -i ../data/ica_all.h5 -o ../result/ica_all -b 32 -l 1024 128 1024 -a rmsprop -w weight_adam_64_32_64_g1_c1_b32.h5 -f 10x_h5
#python3.6 ./main.py -t train -i ../data/ica_all.h5 -o ../result/ica_all -b 128 -l 1024 128 1024 -a rmsprop -w weight_adam_64_32_64_g1_c1_b32.h5 -f 10x_h5
#python3.6 ./main.py -t train -i ../data/ica_all.h5 -o ../result/ica_all -b 512 -l 1024 128 1024 -a rmsprop -w weight_adam_64_32_64_g1_c1_b32.h5 -f 10x_h5
#python3.6 ./main.py -t train -i ../data/ica_all.h5 -o ../result/ica_all -b 1024 -l 1024 128 1024 -a rmsprop -w weight_adam_64_32_64_g1_c1_b32.h5 -f 10x_h5
