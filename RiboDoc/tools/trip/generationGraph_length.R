# Libraries
library("optparse")

# Arguments / Options
option_list = list(
  make_option(c("-w", "--work_dir"), type="character",
              help="Path to working directory"),
  make_option(c("-s", "--sample_name"), type="character",
              help="Sample name")
)

opt = parse_args(OptionParser(option_list=option_list))
work_dir <- opt$w
sample_name <- opt$s

########### Reads Lengths Repartition ###########
path <- paste0(work_dir,"RESULTS/qualitativeAnalysis/")
readsLength <- read.table(file = paste0(path,"readsLengthRepartition/",sample_name,".readsLengthRepartition.txt"))

jpeg(filename = paste0(path,"graphes/readsLengthRepartition/",sample_name,".readsLengthRepartition.jpeg"))
  barplot(readsLength$V2,
          names.arg = readsLength$V1,
          xlab = "Read lengths",
          ylab = "Number of reads")
dev.off()
