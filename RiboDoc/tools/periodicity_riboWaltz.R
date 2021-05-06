library(riboWaltz)
library(stringr)
library(data.table)

dir.create("/data/RESULTS/riboWaltz/")

system(paste0("cat /proc/",Sys.getpid(),"/status | grep VmSize"))

params <- scan(file = "/data/config.yaml",
               what = "character",
               sep = ":"
)

gff_file <- gsub(" ", "", params[which(params=="gff")+1], fixed = TRUE)

# Creates annotation table by transcript names
annotation_db <- riboWaltz::create_annotation(paste0("/data/database/exons_",gff_file,".gtf"))
annotation_db_transcript <- data.table(annotation_db)


# Bam files to be computed
bam_folder <- "/data/RESULTS/BAM_transcriptome/"

bam_list <- list.files(bam_folder)
bam_list <- bam_list[seq(1, length(bam_list), 3)]

samples <- str_replace(bam_list, ".[0-9]{1,3}-[0-9]{1,3}.bam", "")
names(samples) <- str_remove(bam_list, ".bam")
samples

reads_list <- try(riboWaltz::bamtolist(bamfolder = bam_folder, annotation = annotation_db_transcript, name_samples = samples))
reads_list <- riboWaltz::bamtolist(bamfolder = bam_folder, annotation = annotation_db_transcript, name_samples = samples)


# p-site calculation
psite_offset <- riboWaltz::psite(reads_list,
                                 flanking = 0,
                                 start = TRUE,
                                 extremity = "auto",
                                 plot = TRUE,
                                 plot_dir = "/data/RESULTS/riboWaltz/",
                                 plot_format = "png",
                                 cl = 100,
                                 log_file = TRUE,
                                 log_file_dir = "/data/RESULTS/riboWaltz/"
                                 )

reads_psite_list <- riboWaltz::psite_info(reads_list,
                                          psite_offset
                                          )


codon_coverage <- codon_coverage(reads_psite_list, annotation_db_transcript, psite = TRUE)
cds_coverage <- cds_coverage(reads_psite_list, annotation_db_transcript)


length_dist <- rlength_distr(reads_list, sample = samples)
length_dist[[paste0("plot_",samples[1])]]
length_dist[[paste0("plot_",samples[3])]]

for(i in 1:length(samples))
{
  dir.create(paste0("/data/RESULTS/riboWaltz/",samples[i],"/"))
  pdf(file=paste0("/data/RESULTS/riboWaltz/",samples[i],"/reads_distribution_",samples[i],".pdf"))
  plot(length_dist[[paste0("plot_",samples[i])]])
  dev.off()
}

psite_region <- region_psite(reads_psite_list, annotation_db_transcript, sample = samples)
psite_region[["plot"]]

pdf(file="/data/RESULTS/riboWaltz/region_psite.pdf")
psite_region[["plot"]]
dev.off()


frames_stratified <- frame_psite_length(reads_psite_list, sample = samples, region = "all")
frames_stratified[["plot"]]

pdf(file="/data/RESULTS/riboWaltz/frame_psite_length.pdf")
frames_stratified[["plot"]]
dev.off()

frames <- frame_psite(reads_psite_list, sample = samples, region = "all")
frames[["plot"]]

pdf(file="/data/RESULTS/riboWaltz/frame_psite.pdf")
frames[["plot"]]
dev.off()

metaprofile <- metaprofile_psite(reads_psite_list, annotation_db_transcript, sample = samples,
                                  plot_title = "sample.transcript")
metaprofile[["plot_Mut.1"]]
metaprofile[["plot_WT.1"]]

for(i in 1:length(samples))
{
  pdf(file=paste0("/data/RESULTS/riboWaltz/",samples[i],"/metaprofile_psite_",samples[i],".pdf"))
  plot(metaprofile[[paste0("plot_",samples[i])]])
  dev.off()
}

write.csv(psite_offset, "/data/RESULTS/riboWaltz/psite_offset.csv", quote = F, row.names = F)
for(i in 1:length(samples))
{
  write.csv(reads_list[[i]], paste0("/data/RESULTS/riboWaltz/",samples[i],"/reads_list_",samples[i],".csv"), quote = F, row.names = F)
  write.csv(reads_psite_list[[i]], paste0("/data/RESULTS/riboWaltz/",samples[i],"/reads_psite_list_",samples[i],".csv"), quote = F, row.names = F)
}
write.csv(psite_region$dt, "/data/RESULTS/riboWaltz/psite_region.csv", quote = F)
write.csv(metaprofile$dt, "/data/RESULTS/riboWaltz/metaprofile.csv", quote = F, row.names = F)
write.csv(codon_coverage, "/data/RESULTS/riboWaltz/codon_coverage.csv", quote = F, row.names = F)
write.csv(cds_coverage, "/data/RESULTS/riboWaltz/cds_coverage.csv", quote = F, row.names = F)
write.csv(length_dist$dt, "/data/RESULTS/riboWaltz/length_dist.csv", quote = F, row.names = F)
write.csv(frames_stratified$dt, "/data/RESULTS/riboWaltz/frames_stratified.csv", quote = F, row.names = F)
write.csv(frames$dt, "/data/RESULTS/riboWaltz/frames.csv", quote = F, row.names = F)

