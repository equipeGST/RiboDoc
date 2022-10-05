library(riboWaltz)
library(stringr)
library(data.table)


local_path <- "/data/"
dir.create(paste0(local_path, "RESULTS/riboWaltz/"))

params <- scan(file = paste0(local_path, "config.yaml"),
               what = "character",
               sep = ":"
)
args <- commandArgs(trailingOnly=TRUE)

gtf_file <- args[1]


# Creates annotation table by transcript names
annotation_db <- riboWaltz::create_annotation(gtf_file)
annotation_db_transcript_with_cds0l <- data.table(annotation_db)
annotation_db_transcript <- annotation_db_transcript_with_cds0l[annotation_db_transcript_with_cds0l$l_cds > 0,]
# Free unused memory
rm(list=c("annotation_db","annotation_db_transcript_with_cds0l"))
gc()


# Bam files to be computed
bam_folder <- paste0(local_path, "RESULTS/BAM_transcriptome/")

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
  dir.create(paste0(local_path, "RESULTS/riboWaltz/",samples[i],"/"), showWarnings = F)
  tiff(file=paste0(local_path, "RESULTS/riboWaltz/",samples[i],"/reads_distribution_",samples[i],".tiff"))
  plot(length_dist[[paste0("plot_",samples_renamed[i])]])
  dev.off()
}
write.table(length_dist$dt, paste0(local_path, "RESULTS/riboWaltz/length_dist.csv"), quote = F, row.names = F, sep ="\t")
rm(length_dist)
gc()


# p-site offset calculation
source(paste0("/RiboDoc/RiboDoc/tools/ribowaltz_psite_with_NA_control.R"))
psite_offset <- psite_ribowaltz(reads_list,
                                flanking = 6,
                                start = TRUE,
                                extremity = "auto",
                                plot = TRUE,
                                plot_dir = paste0(local_path, "RESULTS/riboWaltz/"),
                                plot_format = "tiff",
                                cl = 100,
                                txt = TRUE,
                                txt_file = paste0(local_path, "RESULTS/riboWaltz/best_offset.txt")
)

reads_psite_list <- riboWaltz::psite_info(reads_list, psite_offset)
write.table(psite_offset, paste0(local_path, "RESULTS/riboWaltz/psite_offset.csv"), quote = F, row.names = F, sep ="\t")
for(i in 1:length(samples_renamed))
{
  write.table(reads_list[[i]], paste0(local_path, "RESULTS/riboWaltz/",samples[i],"/reads_list_",samples[i],".csv"), quote = F, row.names = F, sep ="\t")
}
rm(list=c("psite_offset","reads_list"))
gc()


# Codon coverage
codon_coverage <- codon_coverage(reads_psite_list, annotation_db_transcript, psite = TRUE)
write.table(codon_coverage, paste0(local_path, "RESULTS/riboWaltz/codon_coverage.csv"), quote = F, row.names = F, sep ="\t")
rm(codon_coverage)
gc()


# CDS coverage
cds_coverage <- cds_coverage(reads_psite_list, annotation_db_transcript)
write.table(cds_coverage, paste0(local_path, "RESULTS/riboWaltz/cds_coverage.csv"), quote = F, row.names = F, sep ="\t")
rm(cds_coverage)
gc()


# Proportion of region covered by P-sites
psite_region <- region_psite(reads_psite_list, annotation_db_transcript, sample = samples_renamed)

tiff(file=paste0(local_path, "RESULTS/riboWaltz/region_psite.tiff"))
psite_region[["plot"]]
dev.off()

write.table(psite_region$dt, paste0(local_path, "RESULTS/riboWaltz/psite_region.csv"), quote = F, sep ="\t")
rm(psite_region)
gc()


# Phasing by read length
frames_stratified <- frame_psite_length(reads_psite_list, sample = samples_renamed, region = "all")

tiff(file=paste0(local_path, "RESULTS/riboWaltz/frame_psite_length.tiff"))
frames_stratified[["plot"]]
dev.off()

write.table(frames_stratified$dt, paste0(local_path, "RESULTS/riboWaltz/frames_stratified.csv"), quote = F, row.names = F, sep ="\t")
rm(frames_stratified)
gc()


# Global phasing
frames <- frame_psite(reads_psite_list, sample = samples_renamed, region = "all")

tiff(file=paste0(local_path, "RESULTS/riboWaltz/frame_psite.tiff"))
frames[["plot"]]
dev.off()

write.table(frames$dt, paste0(local_path, "RESULTS/riboWaltz/frames.csv"), quote = F, row.names = F, sep ="\t")
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
  dir.create(paste0(local_path, "RESULTS/riboWaltz/", samples[i],"/results_by_length/"), showWarnings = F)
  tiff(file=paste0(local_path, "RESULTS/riboWaltz/",samples[i],"/metaprofile_psite_-", window_utr, "+", window_cds, ".tiff"))
  plot(metaprofile[[paste0("plot_",samples_renamed[i])]])
  dev.off()
  
  # Metaprofiles by length
  for(len in readsLength_min:readsLength_max) {
    dir.create(paste0(local_path, "RESULTS/riboWaltz/", samples[i],"/results_by_length/metaprofiles_-", window_utr, "+", window_cds, "/"), showWarnings = F)
    dir.create(paste0(local_path, "RESULTS/riboWaltz/", samples[i],"/results_by_length/reads_psite/"), showWarnings = F)
    reads_psite_list_specific_length <- setNames(list(reads_psite_list[[samples_renamed[i]]][length==len]),samples_renamed[i])
    
    metaprofile_specific <- metaprofile_psite(reads_psite_list_specific_length,
                                              annotation_db_transcript,
                                              sample = samples_renamed[i],
                                              utr5l = window_utr, utr3l = window_utr,
                                              cdsl = window_cds,
                                              plot_title = "sample.transcript")
    
    tiff(file=paste0(local_path, "RESULTS/riboWaltz/", samples[i],"/results_by_length/metaprofiles_-", window_utr, "+", window_cds, "/metaprofile_psite_length", len, "_-", window_utr, "+", window_cds, ".tiff"))
    plot(metaprofile_specific[[paste0("plot_",samples_renamed[i])]])
    dev.off()
    
    write.table(metaprofile_specific$dt,
                paste0(local_path, "RESULTS/riboWaltz/", samples[i], "/results_by_length/metaprofiles_-", window_utr, "+", window_cds, "/metaprofile_psite_length", len, "_-", window_utr, "+", window_cds, ".csv"),
                quote = F, row.names = F, sep ="\t")

    write.table(reads_psite_list_specific_length[[samples_renamed[i]]],
                paste0(local_path, "RESULTS/riboWaltz/", samples[i],"/results_by_length/reads_psite/reads_psite_list_specific_length_", len, ".csv"),
                quote = F, row.names = F, sep ="\t")
  }
}
# write.table(metaprofile$dt, paste0(local_path, "RESULTS/riboWaltz/metaprofile.csv"), quote = F, row.names = F, sep ="\t")
# rm(list=c("metaprofile", "metaprofile_specific","reads_psite_list_specific_length"))
# gc()
# 
# 
# for(i in 1:length(samples))
# {
#   write.table(reads_psite_list[[i]], paste0(local_path, "RESULTS/riboWaltz/",samples[i],"/reads_psite_list_",samples[i],".csv"), quote = F, row.names = F, sep ="\t")
# }
# rm(reads_psite_list)
# gc()
