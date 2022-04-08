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


# Bam files to be computed
bam_folder <- paste0(local_path, "RESULTS/BAM_transcriptome/")

bam_list <- list.files(bam_folder)
bam_list <- bam_list[str_detect(bam_list, ".bam$")]

samples <- str_replace(bam_list, ".[0-9]{1,3}-[0-9]{1,3}.bam", "")
names(samples) <- str_remove(bam_list, ".bam")
samples

reads_list <- try(riboWaltz::bamtolist(bamfolder = bam_folder, annotation = annotation_db_transcript, name_samples = samples))
reads_list <- riboWaltz::bamtolist(bamfolder = bam_folder, annotation = annotation_db_transcript, name_samples = samples)


# p-site calculation
source(paste0("/RiboDoc/RiboDoc/tools/ribowaltz_psite_with_NA_control.R"))
psite_offset <- psite_ribowaltz(reads_list,
                                 flanking = 6,
                                 start = TRUE,
                                 extremity = "auto",
                                 plot = TRUE,
                                 plot_dir = paste0(local_path, "RESULTS/riboWaltz/"),
                                 plot_format = "png",
                                 cl = 100,
                                 txt = TRUE,
                                 txt_file = paste0(local_path, "RESULTS/riboWaltz/best_offset.txt")
)

reads_psite_list <- riboWaltz::psite_info(reads_list,
                                          psite_offset
)


codon_coverage <- codon_coverage(reads_psite_list, annotation_db_transcript, psite = TRUE)
cds_coverage <- cds_coverage(reads_psite_list, annotation_db_transcript)


length_dist <- rlength_distr(reads_list, sample = samples)

for(i in 1:length(samples))
{
  dir.create(paste0(local_path, "RESULTS/riboWaltz/",samples[i],"/"), showWarnings = F)
  pdf(file=paste0(local_path, "RESULTS/riboWaltz/",samples[i],"/reads_distribution_",samples[i],".pdf"))
  plot(length_dist[[paste0("plot_",samples[i])]])
  dev.off()
}

psite_region <- region_psite(reads_psite_list, annotation_db_transcript, sample = samples)

pdf(file=paste0(local_path, "RESULTS/riboWaltz/region_psite.pdf"))
psite_region[["plot"]]
dev.off()


frames_stratified <- frame_psite_length(reads_psite_list, sample = samples, region = "all")

pdf(file=paste0(local_path, "RESULTS/riboWaltz/frame_psite_length.pdf"))
frames_stratified[["plot"]]
dev.off()

frames <- frame_psite(reads_psite_list, sample = samples, region = "all")

pdf(file=paste0(local_path, "RESULTS/riboWaltz/frame_psite.pdf"))
frames[["plot"]]
dev.off()

window_utr <- as.integer(gsub(" ", "", params[which(params=="window_utr")+1], fixed = TRUE))
window_cds <- as.integer(gsub(" ", "", params[which(params=="window_cds")+1], fixed = TRUE))
metaprofile <- metaprofile_psite(reads_psite_list,
                                 annotation_db_transcript,
                                 sample = samples,
                                 utr5l = window_utr, utr3l = window_utr,
                                 cdsl = window_cds,
                                 plot_title = "sample.transcript")

for(i in 1:length(samples))
{
  pdf(file=paste0(local_path, "RESULTS/riboWaltz/",samples[i],"/metaprofile_psite_",samples[i],".pdf"))
  plot(metaprofile[[paste0("plot_",samples[i])]])
  dev.off()
}

write.table(psite_offset, paste0(local_path, "RESULTS/riboWaltz/psite_offset.csv"), quote = F, row.names = F, sep ="\t")
for(i in 1:length(samples))
{
  write.table(reads_list[[i]], paste0(local_path, "RESULTS/riboWaltz/",samples[i],"/reads_list_",samples[i],".csv"), quote = F, row.names = F, sep ="\t")
  write.table(reads_psite_list[[i]], paste0(local_path, "RESULTS/riboWaltz/",samples[i],"/reads_psite_list_",samples[i],".csv"), quote = F, row.names = F, sep ="\t")
}
write.table(psite_region$dt, paste0(local_path, "RESULTS/riboWaltz/psite_region.csv"), quote = F, sep ="\t")
write.table(metaprofile$dt, paste0(local_path, "RESULTS/riboWaltz/metaprofile.csv"), quote = F, row.names = F, sep ="\t")
write.table(codon_coverage, paste0(local_path, "RESULTS/riboWaltz/codon_coverage.csv"), quote = F, row.names = F, sep ="\t")
write.table(cds_coverage, paste0(local_path, "RESULTS/riboWaltz/cds_coverage.csv"), quote = F, row.names = F, sep ="\t")
write.table(length_dist$dt, paste0(local_path, "RESULTS/riboWaltz/length_dist.csv"), quote = F, row.names = F, sep ="\t")
write.table(frames_stratified$dt, paste0(local_path, "RESULTS/riboWaltz/frames_stratified.csv"), quote = F, row.names = F, sep ="\t")
write.table(frames$dt, paste0(local_path, "RESULTS/riboWaltz/frames.csv"), quote = F, row.names = F, sep ="\t")

