#########################
####### LIBRARIES #######
#########################

library("optparse")
library("stringr")
library("data.table")
library("riboWaltz")


##########################
####### PARAMETERS #######
##########################

option_list = list(
  make_option(c("-w", "--work_dir"), type="character",
              help="Path to working directory"),
  make_option(c("-f", "--feature_type"), type="character",
              help="Path to designed R functions"),
  make_option(c("-g", "--gtf"), type="character",
              help="Path and name of the input GTF file")
)
opt = parse_args(OptionParser(option_list=option_list))

# Path to wotking directory
local_path <- opt$w

# Path to remastered psite() function from riboWaltz
function_psite <- paste0(opt$f,"ribowaltz_psite_with_NA_control.R")

#  Name of the gtf file
gtf_file <- opt$g

# Read the config file
params <- scan(file = paste0(local_path, "config.yaml"),
               what = "character",
               sep = ":"
)

readsLength_min <- gsub(" ", "", params[which(params=="readsLength_min")+1], fixed = TRUE)
readsLength_max <- gsub(" ", "", params[which(params=="readsLength_max")+1], fixed = TRUE)

ribowaltz_folder <- paste0(local_path, "RESULTS/riboWaltz.", readsLength_min, "-", readsLength_max, "/")
dir.create(ribowaltz_folder)


#########################
####### EXECUTION #######
#########################

# Creates annotation table by transcript names
annotation_db <- riboWaltz::create_annotation(gtf_file)
annotation_db_transcript_with_cds0l <- data.table(annotation_db)
annotation_db_transcript <- annotation_db_transcript_with_cds0l[annotation_db_transcript_with_cds0l$l_cds > 0,]
# Free unused memory
rm(list=c("annotation_db","annotation_db_transcript_with_cds0l"))
gc()


# Bam files to be computed
bam_folder <- paste0(local_path,"RESULTS/BAM_transcriptome.",readsLength_min,"-",readsLength_max,"/")

bam_list <- list.files(bam_folder, pattern = "\\.bam$")

samples <- str_replace(bam_list, ".[0-9]{1,3}-[0-9]{1,3}.bam", "")
names(samples) <- str_remove(bam_list, ".bam")
samples


# List of reads coordinates
reads_list <- try(riboWaltz::bamtolist(bamfolder = bam_folder, annotation = annotation_db_transcript, name_samples = samples))
reads_list <- riboWaltz::bamtolist(bamfolder = bam_folder, annotation = annotation_db_transcript, name_samples = samples)

# Changing dashes ("-") in names to underscores ("_") as ggplot uses parse which throws an error with digits followed by dashes
samples_renamed <- gsub("-","_", samples)
names(samples_renamed) <- gsub("-","_", names(samples_renamed))
# names(samples) <- gsub("([[:digit:]])_([[:digit:]]+)$","\\1-\\2", names(samples), perl = T)
names(reads_list) <- gsub("-","_", names(reads_list))


# Reads lengths distribution
length_dist <- rlength_distr(reads_list, sample = samples_renamed)
for(i in 1:length(samples_renamed))
{
  dir.create(paste0(ribowaltz_folder,samples[i],"/"), showWarnings = F)
  col_selec <- c(1,(i*2),(i*2)+1)
  write.table(length_dist$dt[,..col_selec],
              paste0(ribowaltz_folder,samples[i],"/reads_distribution_",samples[i],".csv"),
              quote = F, row.names = F, sep ="\t")
  tiff(file=paste0(ribowaltz_folder,samples[i],"/reads_distribution_",samples[i],".tiff"))
  plot(length_dist[[paste0("plot_",samples_renamed[i])]])
  dev.off()
}
rm(length_dist)
gc()


# p-site offset calculation
source(function_psite)
psite_offset <- psite_ribowaltz(reads_list,
                                flanking = 6,
                                start = TRUE,
                                extremity = "auto",
                                plot = TRUE,
                                plot_dir = ribowaltz_folder,
                                plot_format = "tiff",
                                cl = 100,
                                txt = TRUE,
                                txt_file = paste0(ribowaltz_folder, "best_offset.csv")
)

reads_psite_list <- riboWaltz::psite_info(reads_list, psite_offset)
write.table(psite_offset, paste0(ribowaltz_folder, "psite_offset.csv"), quote = F, row.names = F, sep ="\t")
rm(list=c("psite_offset","reads_list"))
gc()


# Proportion of region covered by P-sites
psite_region <- region_psite(reads_psite_list, annotation_db_transcript, sample = samples_renamed)

tiff(file=paste0(ribowaltz_folder, "region_psite.tiff"))
psite_region[["plot"]]
dev.off()

rm(psite_region)
gc()


# Phasing by read length
frames_stratified <- frame_psite_length(reads_psite_list, sample = samples_renamed, region = "all")

tiff(file=paste0(ribowaltz_folder, "frame_psite_length.tiff"))
frames_stratified[["plot"]]
dev.off()

rm(frames_stratified)
gc()


# Global phasing
frames <- frame_psite(reads_psite_list, sample = samples_renamed, region = "all")

tiff(file=paste0(ribowaltz_folder, "frame_psite.tiff"))
frames[["plot"]]
dev.off()

rm(frames)
gc()


# Global metaprofiles
window_utr <- as.integer(gsub(" ", "", params[which(params=="window_utr")+1], fixed = TRUE))
window_cds <- as.integer(gsub(" ", "", params[which(params=="window_cds")+1], fixed = TRUE))
metaprofile <- metaprofile_psite(reads_psite_list,
                                 annotation_db_transcript,
                                 sample = samples_renamed,
                                 utr5l = window_utr, utr3l = window_utr,
                                 cdsl = window_cds,
                                 plot_title = "sample.transcript")

readsLength_min <- as.integer(gsub(" ", "", params[which(params=="readsLength_min")+1], fixed = TRUE))
readsLength_max <- as.integer(gsub(" ", "", params[which(params=="readsLength_max")+1], fixed = TRUE))
for(i in 1:length(samples_renamed))
{
  dir.create(paste0(ribowaltz_folder, samples[i],"/results_by_length/"), showWarnings = F)
  tiff(file=paste0(ribowaltz_folder,samples[i],"/metaprofile_psite_-", window_utr, "+", window_cds, ".tiff"))
  plot(metaprofile[[paste0("plot_",samples_renamed[i])]])
  dev.off()
  
  # Metaprofiles by length
  for(len in readsLength_min:readsLength_max) {
    dir.create(paste0(ribowaltz_folder, samples[i],"/results_by_length/metaprofiles_-", window_utr, "+", window_cds, "/"), showWarnings = F)
    dir.create(paste0(ribowaltz_folder, samples[i],"/results_by_length/reads_psite/"), showWarnings = F)
    reads_psite_list_specific_length <- setNames(list(reads_psite_list[[samples_renamed[i]]][length==len]),samples_renamed[i])
    
    metaprofile_specific <- metaprofile_psite(reads_psite_list_specific_length,
                                              annotation_db_transcript,
                                              sample = samples_renamed[i],
                                              utr5l = window_utr, utr3l = window_utr,
                                              cdsl = window_cds,
                                              plot_title = "sample.transcript")
    
    tiff(file=paste0(ribowaltz_folder, samples[i],"/results_by_length/metaprofiles_-", window_utr, "+", window_cds, "/metaprofile_psite_length", len, "_-", window_utr, "+", window_cds, ".tiff"))
    plot(metaprofile_specific[[paste0("plot_",samples_renamed[i])]])
    dev.off()
    
    write.table(metaprofile_specific$dt,
                paste0(ribowaltz_folder, samples[i], "/results_by_length/metaprofiles_-", window_utr, "+", window_cds, "/metaprofile_psite_length", len, "_-", window_utr, "+", window_cds, ".csv"),
                quote = F, row.names = F, sep ="\t")

    write.table(reads_psite_list_specific_length[[samples_renamed[i]]],
                paste0(ribowaltz_folder, samples[i],"/results_by_length/reads_psite/reads_psite_list_specific_length_", len, ".csv"),
                quote = F, row.names = F, sep ="\t")
  }
}