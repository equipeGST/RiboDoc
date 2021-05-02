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
#transcript_id <- str_replace(annotation_db$transcript, "(ENST)", "transcript:\\1")
#annotation_db_transcript$transcript <- transcript_id


# Bam files to be computed
bam_folder <- "/data/RESULTS/BAM_transcriptome/"

bam_list <- list.files(bam_folder)
bam_list <- bam_list[seq(1, length(bam_list), 3)]

samples <- str_replace(bam_list, ".[0-9]{1,3}-[0-9]{1,3}.bam", "")
names(samples) <- str_remove(bam_list, ".bam")
samples

#bam_names <- str_remove(bam_list, ".bam")
#bams <- c("Treated1","Treated2","Treated3","Control1","Control2","Control3")
##bams <- c(rep("Treated",3),rep("Control",3))
#names(bams) <- bam_names

#test <- c("WT.1","Mut.1","WT.3")
#test <- c("Control.1","Treated.1")
#names(test) <- c(str_remove(bam_list[4], ".bam"), str_remove(bam_list[1], ".bam"))

reads_list <- riboWaltz::bamtolist(bamfolder = bam_folder, annotation = annotation_db_transcript, name_samples = samples)
reads_list <- riboWaltz::bamtolist(bamfolder = bam_folder, annotation = annotation_db_transcript, name_samples = samples)
#reads_list_test <- riboWaltz::bamtolist(bamfolder = bam_folder, name_samples = test,  annotation = annotation_db_transcript)


# p-site calculation
psite_offset <- riboWaltz::psite(reads_list,
                                 flanking = 0,
#                                 flanking = 6,
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
                                          #site = c("psite","esite","asite"),
                                          #fastapath = "database/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa",
                                          #fasta_genome = FALSE
                                          #gtfpath = "database/exons_Saccharomyces_cerevisiae.R64-1-1.103_no_description.gtf"
                                          )


codon_coverage <- codon_coverage(reads_psite_list, annotation_db_transcript, psite = TRUE)
cds_coverage <- cds_coverage(reads_psite_list, annotation_db_transcript)


length_dist <- rlength_distr(reads_list, sample = samples)
                            #, cl = 99)
#length_dist_test <- rlength_distr(reads_list_test, sample = test)
length_dist[[paste0("plot_",samples[1])]]
length_dist[[paste0("plot_",samples[3])]]

for(i in 1:length(samples))
{
  dir.create(paste0("/data/RESULTS/riboWaltz/",samples[i],"/"))
  pdf(file=paste0("/data/RESULTS/riboWaltz/",samples[i],"/reads_distribution_",samples[i],".pdf"))
  plot(length_dist[[paste0("plot_",samples[i])]])
  dev.off()
}

#example_ends_heatmap <- rends_heat(reads_list, annotation_db_transcript, sample = samples)
##example_ends_heatmap <- rends_heat(reads_list_test, annotation_db_transcript, sample = test)
#                                   #cl = 85, utr5l = 25, cdsl = 40, utr3l = 25)
#example_ends_heatmap[["plot"]]


psite_region <- region_psite(reads_psite_list, annotation_db_transcript, sample = samples)
#psite_region <- region_psite(reads_psite_list_test, annotation_db_transcript, sample = test)
psite_region[["plot"]]

pdf(file="/data/RESULTS/riboWaltz/region_psite.pdf")
psite_region[["plot"]]
dev.off()


frames_stratified <- frame_psite_length(reads_psite_list, sample = samples, region = "all")
#frames_stratified <- frame_psite_length(reads_psite_list_test, sample = test, region = "all")
                                        #cl = 90)
frames_stratified[["plot"]]

pdf(file="/data/RESULTS/riboWaltz/frame_psite_length.pdf")
frames_stratified[["plot"]]
dev.off()

frames <- frame_psite(reads_psite_list, sample = samples, region = "all")
#frames <- frame_psite(reads_psite_list_test, sample = test, region = "all")
frames[["plot"]]

pdf(file="/data/RESULTS/riboWaltz/frame_psite.pdf")
frames[["plot"]]
dev.off()

metaprofile <- metaprofile_psite(reads_psite_list, annotation_db_transcript, sample = samples,
                                  plot_title = "sample.transcript")
                                  #utr5l = 20, cdsl = 40, utr3l = 20)
#metaprofile <- metaprofile_psite(reads_psite_list_test, annotation_db_transcript, sample = test,
#                                  plot_title = "sample.transcript")
metaprofile[["plot_Mut.1"]]
metaprofile[["plot_WT.1"]]

for(i in 1:length(samples))
{
  pdf(file=paste0("/data/RESULTS/riboWaltz/",samples[i],"/metaprofile_psite_",samples[i],".pdf"))
  plot(metaprofile[[paste0("plot_",samples[i])]])
  dev.off()
}

#codon_usage <- codon_usage_psite(reads_psite_list, annotation_db_transcript, sample = test,
#                                          fastapath = "database/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa",
#                                          fasta_genome = FALSE,
#                                          frequency_normalization = FALSE)
#                                          #frequency_normalization = TRUE)
#codon_usage[["plot"]]


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

