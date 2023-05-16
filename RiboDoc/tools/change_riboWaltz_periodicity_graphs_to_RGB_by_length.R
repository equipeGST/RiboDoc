# Load args : path to working directory
args <- commandArgs(trailingOnly=TRUE)
local_path <- args[1]

# Load parameters from configuration file
params <- scan(file = paste0(local_path, "config.yaml"),
               what = "character",
               sep = ":"
)

qualitative_analysis <- gsub(" ", "", params[which(params=="qualitative_analysis")+1], fixed = TRUE)

readsLength_min <- as.integer(gsub(" ", "", params[which(params=="readsLength_min")+1], fixed = TRUE))
readsLength_max <- as.integer(gsub(" ", "", params[which(params=="readsLength_max")+1], fixed = TRUE))

window_utr <- as.integer(gsub(" ", "", params[which(params=="window_utr")+1], fixed = TRUE))
window_cds <- as.integer(gsub(" ", "", params[which(params=="window_cds")+1], fixed = TRUE))

# Define the order of colors, depending on the phase it is in (assuming start and stop are in the same phase)
if(window_utr %% 3 == 0) {
  phase_color_start <- c('red','green','blue')
} else if(window_utr %% 3 == 1) {
  phase_color_start <- c('blue','red','green')
} else {
  phase_color_start <- c('green','blue','red')
}
if(window_cds %% 3 == 0) {
  phase_color_stop <- c('blue','red','green')
} else if(window_cds %% 3 == 1) {
  phase_color_stop <- c('green','blue','red')
} else {
  phase_color_stop <- c('red','green','blue')
}


# Add the "periodicity" folder to the "RESULTS" folder from RiboDoc
dir.create(paste0(local_path, "RESULTS/periodicity_-", window_utr, "+", window_cds, "/"), showWarnings = F)

samples <- list.dirs(paste0(local_path, "RESULTS/riboWaltz.", readsLength_min, "-", readsLength_max, "/"), full.names = F, recursive = F)

# For each sample
for(sample in samples) {
  print(sample)
  dir.create(paste0(local_path, "RESULTS/periodicity_-", window_utr, "+", window_cds, "/", sample, "/"),showWarnings = F)

  # For each length in the desired window (in config file)
  for(specific_length in readsLength_min:readsLength_max) {
    
    # Load data
    pathway_metaprofile_table_specific <- paste0(local_path, "RESULTS/riboWaltz.", readsLength_min, "-", readsLength_max, "/", sample,"/results_by_length/metaprofiles_-", window_utr,"+", window_cds, "/metaprofile_psite_length", specific_length, "_-", window_utr,"+", window_cds, ".csv")
    perio_specific <- read.table(pathway_metaprofile_table_specific, header = TRUE, sep = "\t")
    
    # Select relative positions from start or stop
    perio_start_specific <- perio_specific[which(perio_specific$reg == "Distance from start (nt)"),]
    perio_stop_specific <- perio_specific[which(perio_specific$reg == "Distance from stop (nt)"),]
    
    # Filter to fit the window of interest
    perio_start_filtered_specific <- perio_start_specific[which(perio_start_specific$distance > -window_utr & perio_start_specific$distance < window_cds),]
    perio_stop_filtered_specific <- perio_stop_specific[which(perio_stop_specific$distance > -window_cds & perio_stop_specific$distance < window_utr),]
    
    # Save plots as tiff file
    tiff(paste0(local_path, "RESULTS/periodicity_-", window_utr, "+", window_cds, "/", sample, "/length", specific_length, "_start.tiff"))
      barplot(perio_start_filtered_specific[,gsub("-","_",sample)],
              col = phase_color_start,
              names.arg = perio_start_filtered_specific$distance,
              cex.names = 0.75, las=3,
              xlab = 'Relative positions around start codon', ylab = 'Number of reads',
              border = F)
    dev.off()
    
    tiff(paste0(local_path, "RESULTS/periodicity_-", window_utr, "+", window_cds, "/", sample, "/length", specific_length, "_stop.tiff"))
      barplot(perio_stop_filtered_specific[,gsub("-","_",sample)],
              col = phase_color_stop,
              names.arg = perio_stop_filtered_specific$distance,
              cex.names = 0.75, las=3,
              xlab = 'Relative positions around stop codon', ylab = 'Number of reads',
              border = F)
    dev.off()
  }
}
