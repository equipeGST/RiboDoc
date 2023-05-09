# Path to project directory
local_path <- "/data/"

# Read the config file
params <- scan(file = paste0(local_path, "config.yaml"),
               what = "character",
               sep = ":"
)

# Load designed R functions
source(paste0("/RiboDoc/RiboDoc/tools/DESeq2_analysis_functions.R"))

# Path to count tables folder
paths_list <- DESeq2_folder_paths(local_path)

# Make counts rounded for future differential analyses
expData = round(read.table(paths_list$pathway_matrix, header = T, row.names = 1, check.names = FALSE, sep="\t"))

# Load names corresponding to transcript IDs
names_list = read.table(paths_list$pathway_names, header = T, row.names = 1, check.names = FALSE, sep="\t")

# Create count matrices by gene of by transcript
expData_transcript <- transcript_or_gene(data = expData, list_of_names = names_list)
