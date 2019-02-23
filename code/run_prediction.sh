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
#python3 -m pdb ./main.py -t train -i ../data/cf.10xgenomics.com/samples/cell-exp/1.1.0/293t/filtered_matrices_mex/hg19/ -o ../result/1.1.0/293t -b 128 -l 128 64 32 64 128 -a adam -w weight_adam_64_32_64_g1_c1_b128.h5 -f 10x_mtx
#
##run on ica_all
#python3 ./main.py -t train -i xxx -o xxx -w xxx -f 10x_mtx -a adam
#mkdir ../result/ica_all run on P710
for a in "adam" "rmsprop"
do
	for b in {128,512}
	do
		for l in "128 64 32 64 128"
		do
			px=''
			for w in $l
			do
				px+='_'$w
			done
			weight_file="weight_${a}${px}_g1_c1_b$b.h5"
			gene_file="gene_${a}${px}_g1_c1_b$b.csv"
			python3.6 ./main.py -t train -i ../data/ica_all.h5 -o ../result/ica_all -g ${gene_file} -b $b -l $l -a $a -w ${weight_file} -f 10x_h5 > ../result/ica_all/$a${px}_b$b.log
		done
	done
done
# Generate 10 masking data 

#mkdir ../data/geneFilterResult/maskingData
#mkdir ../data/geneFilterResult/maskingData/h5
#mkdir ../data/geneFilterResult/maskingData/csv
#mkdir ../data/geneFilterResult/maskingData/log

# for i in {1:100}
# do
	# echo "../bin/maskingData ../data/geneFilterResult/293t_filtered_gene_bc_matrices_mex.h5 ../data/geneFilterResult/maskingData/h5/293t_filtered_gene_bc_matrices_mex_${i}.h5"
# done
# python3.6 -m pdb ./main.py -t prediction -i ../data/geneFilterResult/maskingData/h5/293t_filtered_gene_bc_matrices_mex_masking_2_1.h5 -g ../result/ica_all/input_gene_g1.csv -l 64 32 64 -w ../result/ica_all/weight_adam_64_32_64_g1_c1_b32.h5 --filtered --mode denoise -f 10x_h5 -o ./result/t
#for i in $(seq 1 10)
#do
	##wget https://s3.amazonaws.com/saver-vince-out/output/24${i}-job-output.zip
	##mv 24${i}-job-output.zip 293t_filtered_gene_bc_matrices_mex_masking_adam_64_32_64_2_${i}.zip
	#mkdir 293t_filtered_gene_bc_matrices_mex_masking_adam_64_32_64_2_${i}
	#unzip -d ./293t_filtered_gene_bc_matrices_mex_masking_adam_64_32_64_2_${i} 293t_filtered_gene_bc_matrices_mex_masking_adam_64_32_64_2_${i}.zip
#done
#result="../result/zheng/maskingData/myresult"
#l="64 32 64"
#a="adam"
#px="_64_32_64"
#b=32
#weight_file="weight_${a}${px}_g1_c1_b$b.h5"
#for file in $(ls ../data/geneFilterResult/293*.h5)  
#do  
	#folder=${file%/*}
	#filename=${file##*/}
	#filename=${filename%.*}	
	#for j in 2
	#do
		#for i in $(seq 1 10)
		#do
			#h5_files="${folder}/maskingData/h5/${filename}_masking_${j}_${i}.h5"
			##echo "python3.6 ./main.py -t prediction -i ${h5_files} -o ${result}/${filename}_masking_${a}${px}_${j}_${i} -g ../result/ica_all/input_gene_g1.csv -l $l -w ../result/ica_all/${weight_file} -f 10x_h5"
			#python3.6 ./main.py -t prediction -i ${h5_files} -o ${result}/${filename}_masking_${a}${px}_${j}_${i} -g ../result/ica_all/input_gene_g1.csv -l $l -w ../result/ica_all/${weight_file} --filtered --mode full -f 10x_h5
		#done
	#done 
#done

#for file in $(ls ../data/geneFilterResult/293*.h5)  
#do  
	#folder=${file%/*}
	#filename=${file##*/}
	#filename=${filename%.*}	
	#for j in 5
	#do
		#for i in $(seq 1 10)
		#do
			#h5_files="${folder}/maskingData/h5/${filename}_masking_${j}_${i}.h5"
			#csv_files="${folder}/maskingData/csv/${filename}_masking_${j}_${i}.csv"
			#csv_files_t="${folder}/maskingData/csv/${filename}_masking_${j}_${i}_t.csv"
			#result_path="${folder}/maskingData/"
			#log_files="${folder}/maskingData/log/${filename}_masking_${j}_${i}_log.csv"
			#../bin/maskingDataTest ${file} ${j} ${h5_files} ${log_files} $i
			#../bin/write2CSVTest ${h5_files} ${csv_files} "original"
			#../bin/convert_csv_row_to_clumn ${csv_files} ${csv_files_t} &
			##python3.6 ./main.py -t prediction -i ${h5_files} -o ${result} -g ../result/ica_all/input_gene_g1.csv -l $l -w xx -f 10x_h5
		#done
		#wait
	#done 
#done

#declare -a testdata=("293t_filtered_gene_bc_matrices_mex" "b_cells_filtered_matrices_mex" "cd34_filtered_matrices_mex" "cd4_t_helper_filtered_matrices_mex" "cd56_nk_filtered_matrices_mex" "cytotoxic_t_filtered_matrices_mex" "memory_t_filtered_matrices_mex" "naive_cytotoxic_filtered_matrices_mex" "naive_t_filtered_matrices_mex" "regulatory_t_filtered_matrices_mex")

#for i in "${testdata[@]}"
#do
	##echo $i
	#for a in "adam" "rmsprop"
	#do
		##echo $a
		#for l in "64 32 64" "512 64 512"
		#do
			##echo $l
			#px=''
			#for w in $l
			#do
				#px+='_'$w
			#done
			#for b in 32 128 512 1024
			#do
				##echo $b
				#weight_file="weight_${a}${px}_g1_c1_b$b.h5"
				##echo "python3.6 ./main.py -t prediction -i ../data/zheng/${i}/hg19/ -o ../result/zheng/${i}_${a}${px}_${b} -g ../result/ica_all/input_gene_g1.csv -l $l -w ../result/ica_all/${weight_file} -f 10x_mtx &"
				#python3.6 ./main.py -t prediction -i ../data/zheng/${i}/hg19/ -o ../result/zheng/${i}_${a}${px}_${b} -g ../result/ica_all/input_gene_g1.csv -l $l -w ../result/ica_all/${weight_file} -f 10x_mtx &
			#done
		#done
	#done
	#wait
#done

#for a in "${testdata[@]}"
#do
	#echo "$a"
	#printf "%s\n" "$a"
	#echo "python3.6 ./main.py -t prediction -i ../data/zheng/${a}/hg19/ -o ../result/zheng/${a}.mean.csv -g ../result/ica_all/input_gene_g1.csv -w ../result/ica_all/weight_adam_64_32_64_g1_c1_b32.h5 -f 10x_mtx"
#done
###Test prediction model 
#python3.6 ./main.py -t prediction -i ../data/zheng/293t_filtered_gene_bc_matrices_mex/hg19/ -o ../result/zheng/test.csv -g ../result/ica_all/input_gene_g1.csv -w ../result/ica_all/weight_adam_64_32_64_g1_c1_b32.h5 -f 10x_mtx
#python3.6 -m pdb ./main.py -t prediction -i ../data/cf.10xgenomics.com/samples/cell-exp/1.1.0/293t/filtered_matrices_mex/hg19/ -o ../result/ -g ../result/ica_all/input_gene_g1.csv -w ../result/ica_all/weight_adam_64_32_64_g1_c1_b1024.h5 -f 10x_mtx
#python3.6 ./main.py -t prediction -i ../data/cf.10xgenomics.com/samples/cell-exp/1.1.0/293t/filtered_matrices_mex/hg19/ -o ../result/ -g ../result/ica_all/input_gene_g1.csv -w ../result/ica_all/weight_adam_64_32_64_g1_c1_b1024.h5 -f 10x_mtx
###Test training model
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
