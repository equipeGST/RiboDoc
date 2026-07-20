#!/usr/bin/env Rscript

# Calculation of optimal P-site offsets, diagnostic analysis and visual inspection of ribosome profiling data

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(riboWaltz))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))

# =========
# Parse parameters
# =========

option_list = list(
    make_option(c("-D", "--output_dir"), type="character", help="Output Directory"),
    make_option(c("-r", "--reference_condition"), type="character", help="Reference condition"),
    make_option(c("-t", "--threads"), type="integer", help="Number of CPUs"),
    make_option(c("-b", "--bam"), type="character", help="Path to BAM alignment (transcripts) folder"),
    make_option(c("-g", "--gtf"), type="character", help="Path to the input GTF annotation file"),
    make_option(c("-f", "--fasta"), type="character", help="Path to the input FASTA annotation file"),
    make_option(c("-m", "--min_len"), type="integer", help="Minimum read length for P-site identification"),
    make_option(c("-M", "--max_len"), type="integer", help="Maximum read length for P-site identification"),
    make_option(c("-c", "--cds_window"), type="integer", help="Window size for CDS region", default=90),
    make_option(c("-u", "--utr_window"), type="integer", help="Window size for UTR regions", default=30)
)
opt = parse_args(OptionParser(option_list=option_list))

# =========
# Set parameters
# =========

output_dir     = opt$output_dir
gtf            = opt$gtf
bam_folder     = opt$bam
config_file    = opt$config
dir.create(output_dir, showWarnings = FALSE)

palette <- c("#be95c4", "#5fa8d3", "#a7c957", "#ffbd00", "#e63946", "#e76f51")

refCond         <- opt$reference_condition
window_utr      <- opt$utr_window
window_cds      <- opt$cds_window
readsLength_min <- opt$min_len
readsLength_max <- opt$max_len

# =========
# Annotation
# =========

annotation_db <- riboWaltz::create_annotation(gtf)
annotation_db_transcript_with_cds0l <- data.table(annotation_db)
annotation_db_transcript <- annotation_db_transcript_with_cds0l[annotation_db_transcript_with_cds0l$l_cds > 0,]

rm(list=c("annotation_db","annotation_db_transcript_with_cds0l"))
gc()

# =========
# Load BAM files
# =========

bam_list <- list.files(bam_folder, pattern = "_riboseq.*.*.transcripts.bam$")

samples        <- str_replace(bam_list, ".bam", "")
names(samples) <- str_remove(bam_list, ".bam")
print(samples)

# Single load (no double bamtolist call)
reads_list <- riboWaltz::bamtolist(bamfolder = bam_folder, annotation = annotation_db_transcript, name_samples = samples)

# Sanitize names (replace "-" with "_")
samples_renamed        <- gsub("-","_", samples)
names(samples_renamed) <- gsub("-","_", names(samples_renamed))
names(reads_list)      <- gsub("-","_", names(reads_list))

# sample_names : character vector of all individual sample names (length = n samples)
sample_names <- names(reads_list)

# input_samples : named list grouping replicates by condition (length = n conditions)
# e.g. list(WT = c("WT_riboseq.1.transcripts", ...), mutant = c(...))
conditions    <- unique(str_extract(sample_names, "^[^_]+"))
input_samples <- setNames(
    lapply(conditions, function(cond) grep(cond, sample_names, value = TRUE)),
    conditions)

print("Individual sample names:")
print(sample_names)
print("Samples grouped by condition:")
print(input_samples)

# =========
# P-site offset calculation
# =========

psite_offset <- psite(reads_list,
    flanking    = 6,
    start       = TRUE,
    extremity   = "auto",
    plot        = TRUE,
    plot_dir    = output_dir,
    plot_format = "tiff",
    cl          = 100,
    txt         = TRUE,
    txt_file    = paste0(output_dir, "best_offset.tsv"))

reads_psite_list <- psite_info(reads_list, psite_offset)
write.table(psite_offset, paste0(output_dir, "psite_offset.tsv"), quote = F, row.names = F, sep ="\t")
rm(psite_offset)
gc()

# =========
# Read length distribution — averaged per condition (1 plot per condition)
# =========

length_dist_all_conds <- rlength_distr(
    reads_list,
    sample      = input_samples,
    multisamples = "average",
    plot_style  = "dodge",
    cl          = 99,
    colour      = palette)
ggsave(filename = paste0(output_dir, "/read_length_distribution.tiff"),
       plot = length_dist_all_conds[["plot"]], device = "tiff", width = 12, height = 8)
rm(length_dist_all_conds)
gc()

# Read length distribution — one plot and TSV per individual sample
# sample_names used here so that $dt has 2 columns per sample and plots are named "plot_<samplename>"

for(i in seq_along(sample_names))
{
    sname <- sample_names[i]
    print(sname)
    length_dist <- rlength_distr(
        reads_list,
        sample = sname)
    dir.create(paste0(output_dir, sname, "/"), showWarnings = FALSE)
    write.table(length_dist$count_dt,
                paste0(output_dir, sname, "/reads_distribution_", sname, ".tsv"),
                quote = FALSE, row.names = FALSE, sep = "\t")
    ggsave(filename = paste0(output_dir, sname, "/read_length_distribution_", sname, ".tiff"),
           plot    = length_dist[[paste0("plot_", sname)]],
           device  = "tiff", width = 12, height = 8)
}
rm(length_dist)
gc()

# =========
# Read ends heatmap (averaged per condition)
# =========

# ends_heatmap <- rends_heat(
#     reads_list,
#     annotation_db_transcript,
#     sample       = input_samples,
#     cl           = 85,
#     multisamples = "average",
#     utr5l        = 25, cdsl = 40, utr3l = 25,
#     colour       = "#333f50")
# ggsave(filename = paste0(output_dir, "/heatmap_ends.tiff"),
#        plot = ends_heatmap[["plot"]], device = "tiff", width = 12, height = 8)
# rm(ends_heatmap, reads_list)
# gc()

# =========
# Codon and CDS coverage
# =========

codon_coverage_table <- codon_coverage(reads_psite_list, annotation_db_transcript, psite = FALSE)
write.table(codon_coverage_table, file = paste0(output_dir, "codon_coverage.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)

cds_coverage_table <- cds_coverage(reads_psite_list, annotation_db_transcript)
write.table(cds_coverage_table, file = paste0(output_dir, "cds_coverage.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE)

# =========
# P-sites per region (averaged per condition)
# =========

psite_per_region_stack <- riboWaltz::region_psite(
    reads_psite_list,
    annotation_db_transcript,
    sample       = input_samples,
    multisamples = "average",
    plot_style   = "stack",
    cl           = 85,
    colour       = c("#333f50", "gray70", "#39827c"))
ggsave(filename = paste0(output_dir, "/psite_region_stack.tiff"),
       plot = psite_per_region_stack[["plot"]], device = "tiff", width = 12, height = 8)
# Defensive access: try $dt directly if $region$dt does not exist
if (!is.null(psite_per_region_stack$dt)) {
    write.table(psite_per_region_stack$dt,
                paste0(output_dir, "region_psite.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
} else if (!is.null(psite_per_region_stack$region$dt)) {
    write.table(psite_per_region_stack$region$dt,
                paste0(output_dir, "region_psite.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
}
rm(psite_per_region_stack)
gc()

psite_per_region_dodge <- riboWaltz::region_psite(
    reads_psite_list,
    annotation_db_transcript,
    sample       = input_samples,
    multisamples = "average",
    plot_style   = "dodge",
    cl           = 85,
    colour       = c("#333f50", "gray70", "#39827c"))
ggsave(filename = paste0(output_dir, "/psite_region_dodge.tiff"),
       plot = psite_per_region_dodge[["plot"]], device = "tiff", width = 12, height = 8)
rm(psite_per_region_dodge)
gc()

# =========
# Frame / P-site length heatmap
# Condition-level averaged plot + one TSV and TIFF per individual sample
# =========

heatmap_all_conds <- frame_psite_length(
    reads_psite_list,
    annotation_db_transcript,
    sample       = input_samples,
    multisamples = "average",
    plot_style   = "facet",
    region       = "all",
    cl           = 85,
    colour       = "#bc3e46")
ggsave(filename = paste0(output_dir, "/frame_psite_length.tiff"),
       plot = heatmap_all_conds[["plot"]], device = "tiff", width = 10, height = 12)
write.table(heatmap_all_conds$dt, paste0(output_dir, "frame_psite_length.tsv"),
            quote = FALSE, row.names = FALSE, sep = "\t")
rm(heatmap_all_conds)
gc()

# Per-sample frame/length heatmap (called individually so plots and $dt are indexed by sample)
for(sname in sample_names)
{
    print(paste("Processing sample:", sname))
    heatmap_sample <- frame_psite_length(
        reads_psite_list,
        annotation_db_transcript,
        sample     = sname,
        plot_style = "facet",
        region     = "all",
        cl         = 85,
        colour     = "#bc3e46")
    
    # Vérifier si $dt existe
    if (!is.null(heatmap_sample$count_dt)) {
        write.table(heatmap_sample$count_dt,
                    paste0(output_dir, sname, "/heatmap_psite_length_", sname, ".tsv"),
                    quote = FALSE, row.names = FALSE, sep = "\t")
    }
    
    # Vérifier si plot existe
    if (!is.null(heatmap_sample$plot) || !is.null(heatmap_sample[[paste0("plot_", sname)]])) {
        plot_obj <- heatmap_sample$plot %||% heatmap_sample[[paste0("plot_", sname)]]
        ggsave(filename = paste0(output_dir, sname, "/heatmap_psite_length_", sname, ".tiff"),
               plot = plot_obj, device = "tiff", width = 12, height = 8)
    }
    
    rm(heatmap_sample)
    gc()
}

# =========
# Global phasing (averaged per condition)
# =========

frames <- frame_psite(
    reads_psite_list,
    annotation_db_transcript,
    sample       = input_samples,
    multisamples = "average",
    plot_style   = "facet",
    region       = "all",
    colour       = c("#333f50", "#39827c"))
ggsave(filename = paste0(output_dir, "/frame_psite.tiff"),
       plot = frames[["plot"]], device = "tiff", width = 12, height = 8)
write.table(frames$count_dt, paste0(output_dir, "frame_psite.tsv"),
            quote = FALSE, row.names = FALSE, sep = "\t")
rm(frames)
gc()

# =========
# Metaprofile and metaheatmap of P-site density (averaged per condition)
# =========

metaprofile_all_conds <- metaprofile_psite(
    reads_psite_list,
    annotation_db_transcript,
    sample       = input_samples,
    multisamples = "average",
    plot_style   = "overlap",
    utr5l        = window_utr, cdsl = window_cds, utr3l = window_utr,
    colour       = c("#333f50", "#39827c"))
ggsave(filename = paste0(output_dir, "/metaprofile_psite.tiff"),
       plot = metaprofile_all_conds[["plot"]], device = "tiff", width = 12, height = 8)
rm(metaprofile_all_conds)
gc()

metaheatmap_all_conds <- metaheatmap_psite(
    reads_psite_list,
    annotation_db_transcript,
    sample       = input_samples,
    multisamples = "average",
    utr5l        = window_utr, cdsl = window_cds, utr3l = window_utr,
    colour       = "#333f50")
ggsave(filename = paste0(output_dir, "/heatmap_psite.tiff"),
       plot = metaheatmap_all_conds[["plot"]], device = "tiff", width = 12, height = 8)
rm(metaheatmap_all_conds)
gc()

# =========
# Per-sample metaprofile by read length
# =========

for(sname in sample_names)
{
    dir.create(paste0(output_dir, sname, "/results_by_length/metaprofiles_-", window_utr, "+", window_cds, "/"),
               showWarnings = FALSE, recursive = TRUE)
    dir.create(paste0(output_dir, sname, "/results_by_length/reads_psite/"),
               showWarnings = FALSE, recursive = TRUE)

    # Condition-level metaheatmap saved per sample folder
    metaheatmap_sample <- metaheatmap_psite(
        reads_psite_list,
        annotation_db_transcript,
        sample = sname,
        utr5l  = window_utr, cdsl = window_cds, utr3l = window_utr,
        colour = "#333f50")
    ggsave(filename = paste0(output_dir, sname, "/metaprofile_psite_-", window_utr, "+", window_cds, ".tiff"),
           plot = metaheatmap_sample[["plot"]], device = "tiff", width = 12, height = 8)
    rm(metaheatmap_sample)
    gc()

    for(len in readsLength_min:readsLength_max)
    {
        reads_psite_len <- setNames(
            list(reads_psite_list[[sname]][length == len]),
            sname)

        metaprofile_specific <- metaprofile_psite(
            reads_psite_len,
            annotation_db_transcript,
            sample = sname,
            utr5l  = window_utr, utr3l = window_utr, cdsl = window_cds)
        plot_name <- paste0("plot_", sname)

        ggsave(filename = paste0(output_dir, sname, "/results_by_length/metaprofiles_-", window_utr, "+", window_cds,
                                 "/metaprofile_psite_length", len, "_-", window_utr, "+", window_cds, ".tiff"),
               plot = metaprofile_specific[[plot_name]], device = "tiff", width = 12, height = 8)
        
        dt <- metaprofile_specific$count_dt[, c(3,2,4)]
        colnames(dt) <- c("distance", "reg", sname)
        write.table(dt,
                    paste0(output_dir, sname, "/results_by_length/metaprofiles_-", window_utr, "+", window_cds,
                           "/metaprofile_psite_length", len, "_-", window_utr, "+", window_cds, ".tsv"),
                    quote = FALSE, row.names = FALSE, sep = "\t")

        write.table(reads_psite_len[[sname]],
                    paste0(output_dir, sname, "/results_by_length/reads_psite/reads_psite_list_specific_length_", len, ".tsv"),
                    quote = FALSE, row.names = FALSE, sep = "\t")
    }
}

# =========
# Codon usage (averaged per condition)
# =========

# cu_barplot_all_conds <- codon_usage_psite(
#     reads_psite_list,
#     annotation_db_transcript,
#     sample                 = input_samples,
#     multisamples           = "average",
#     plot_style             = "facet",
#     fastapath              = opt$fasta,
#     fasta_genome           = FALSE,
#     frequency_normalization = FALSE)
# ggsave(filename = paste0(output_dir, "/codon_usage_psite.tiff"),
#        plot = cu_barplot_all_conds[["plot"]], device = "tiff", width = 12, height = 8)
# rm(cu_barplot_all_conds)
# gc()