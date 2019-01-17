#!/bin/bash

PROGRAM = ./main.py
python $PROGRAM ../data/ica_bone_marrow_h5.h5 ../result/ica_bone_marrow 10x_h5
python $PROGRAM ../data/ica_cord_blood_h5.h5 ../result/ica_cord_blood 10x_mtx
python $PROGRAM ../data/cf.10xgenomics.com/samples/cell-exp/1.1.0/fresh_68k_pbmc_donor_a/filtered_matrices_mex/hg19/ 
