#!/usr/bin/env Rscript
library('VIPER')
library('argparser')

p <- arg_parser("Run viper")
p <- add_argument(p, "ifile", help="an input file")
p <- add_argument(p, "outdir", help="a output file")
p <- add_argument(p, "prefix", help="a prefix of output files")
p <- add_argument(p, "prefixt", help="a prefix of output files")

argv <- parse_args(p)

inputfile=argv$ifile
outdir=argv$outdir
prefix = argv$prefix
prefixt= argv$prefixt

gene.expression <- read.csv(inputfile)

# Preliminary analysis
# testCell <- PredictCell(gene.expression, GeneNum=1000, ZeroRate=0.05)
# testGene <- PredictGene(gene.expression, GeneNum=500, ZeroRate=0.05)
# PlotR2(testcell, testGene)

system.time(res <- VIPER(gene.expression, num = 5000, percentage.cutoff = 0.5, minbool = FALSE, alpha = 0.5,  report = True, outdir = outdir, prefix = prefix))
system.time(res <- VIPER(t(gene.expression), num = 5000, percentage.cutoff = 0.88, minbool = FALSE, alpha = 0.5,  report = True, outdir = outdir, prefix = prefixt))
