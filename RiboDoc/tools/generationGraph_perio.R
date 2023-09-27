# Libraries
library("optparse")
##########################
####### PARAMETERS #######
##########################
option_list = list(
  make_option(c("-w", "--work_dir"), type="character",
              help="Path to working directory"),
  make_option(c("-n", "--name"), type="character",
              help="File name"),
  make_option(c("-f", "--feature_type"), type="character",
              help="Feature type in GFF3 annotation"),
  make_option(c("-s", "--start_or_stop"), type="character",
              help="Either 'start' or 'stop'"),
  make_option(c("-l", "--length_before"), default=30, type="integer",
              help="Number of nucleotides in UTR")
)

opt = parse_args(OptionParser(option_list=option_list))
work_dir <- opt$w
name_file <- opt$n
feature_type <- opt$f
start_or_stop <- opt$s
length_before <- opt$l
##########################

###########################
####### Periodicity #######
###########################
path <- paste0(work_dir, "RESULTS/qualitativeAnalysis/")
perio <- read.table(file = paste0(path,"periodicity/",name_file,".txt"))

# Define the colors, depending on the length of the region before start or stop
if(start_or_stop == "start"){
  if(length_before %% 3 == 0) {
    phase_color <- c('blue','red','green')
  } else if(length_before %% 3 == 1) {
    phase_color <- c('green','blue','red')
  } else {
    phase_color <- c('red','green','blue')
  }
} else {
  if(length_before %% 3 == 0) {
    phase_color <- c('green','blue','red')
  } else if(length_before %% 3 == 1) {
    phase_color <- c('red','green','blue')
  } else {
    phase_color <- c('blue','red','green')
  }
}

jpeg(filename = paste0(path,"graphes/periodicity/",name_file,".jpeg"))
  barplot(perio$V2,
          col = phase_color,
          names.arg = perio$V1-length_before,
          cex.names = 0.75, las=3,
          xlab = paste0("Relative positions around",feature_type,start_or_stop,"s"),
          ylab = "Number of reads",
          border = F)
dev.off()
###########################
