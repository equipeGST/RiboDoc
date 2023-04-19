#!/usr/bin/env Rscript

# Libraries
library("optparse")

# Arguments / Options
option_list = list(
  make_option(c("-w", "--work_dir"), type="character",
              help="Path to working directory"),
  make_option(c("-n", "--name"), type="character",
              help="File name"),
  make_option(c("-f", "--feature_type"), type="character",
              help="Feature type in GFF3 annotation"),
  make_option(c("-s", "--start_or_stop"), type="character",
              help="Either 'start' or 'stop'"),
  make_option(c("-l", "--length_in_utr"), default=30, type="integer",
              help="Number of nucleotides in UTR")
)

opt = parse_args(OptionParser(option_list=option_list))
work_dir <- opt$w
name_file <- opt$n
feature_type <- opt$f
start_or_stop <- opt$s
length_in_utr <- opt$l
# args = commandArgs(trailingOnly=TRUE)
# work_dir <- args[1]
# name_file <- args[2]
# feature_type <- args[3]
# start_or_stop <- args[4]
# length_in_utr <- args[5]

########### Periodicity ###########
path <- paste0(work_dir, "RESULTS/qualitativeAnalysis/")
perio <- read.table(file = paste0(path,"periodicity/",name_file,".txt"))

jpeg(filename = paste0(path,"graphes/periodicity/",name_file,".jpeg"))
barplot(perio$V2,col = c("red","green","blue"),
        names.arg = paste0(perio$V1,"-",length_in_utr),
        cex.names = 0.75, las=3,
        xlab = paste0("Relative positions around",feature_type,start_or_stop,"s"),
        ylab = "Number of reads",
        border = F)
dev.off()
