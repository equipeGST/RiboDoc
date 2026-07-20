#########################
####### LIBRARIES #######
#########################

library("optparse")


##########################
####### PARAMETERS #######
##########################

option_list = list(
  make_option(c("-w", "--work_dir"), type="character",
              help="Path to working directory"),
  make_option(c("-c", "--config_file"), type="character",
              help="Path to config file"),
  make_option(c("-f", "--feature_type"), type="character",
              help="Path to designed R functions"),
  make_option(c("-t", "--type"), type="character",
              help="Feature type: riboseq or rnaseq")
)
opt = parse_args(OptionParser(option_list=option_list))

# Path to working directory
local_path <- opt$w

# Path to config file
config_file <- opt$c

# Feature type
type <- opt$type

# Path to designed R functions
functions_path <- paste0(opt$f,"DESeq2_analysis_functions.R")

# Read the config file
params <- scan(file = config_file,
               what = "character",
               sep = ":"
)


#########################
####### FUNCTIONS #######
#########################

# Load designed R functions
source(functions_path)


#########################
####### EXECUTION #######
#########################

# Path to count tables folder
paths_list <- DESeq2_folder_paths(local_path, type)

# Make counts rounded for future differential analyses
expData = read.table(paths_list$pathway_matrix, header = T, row.names = 1, check.names = FALSE, sep="\t")

# Create count matrices by gene of by transcript
expData_transcript <- transcript_or_gene(expData)
