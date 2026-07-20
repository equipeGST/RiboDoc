#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(anota2seq))

option_list = list(
    make_option(c("-f", "--folder"), type="character", help="Folder containing the input files"),
    make_option(c("-o", "--output_dir"), type="character", help="Output directory for the results"),
    make_option(c("-r", "--region"), type="character", help="Region used for the analysis, e.g. 'CDS' or '5UTR' or '3UTR'")
    #make_option(c("-R", "--ref"), type="character", help="Reference condition")
    # make_option(c("-s", "--sample_file"), type="character", help="Number of CPUs"),
    # make_option(c("-t", "--sample_treatment_col"), type="character", help="Path to BAM alignment file of transcripts"),
    # make_option(c("-r", "--reference_level"), type="character", help="Path to the input GTF annotation file"),
    # make_option(c("-f", "--fasta"), type="character", help="Path to the input FASTA annotation file"),
    # make_option(c("-l", "--length_range"), type="character", help="Range of read lengths for P-site identification, formatted as two integers separated by a colon (e.g., '26:34'). If unspecified, all lengths are considered"),
    # make_option(c("-p", "--periodicity_threshold"), type="integrer", help="Periodicity threshold for P-site identification to filter out read lengths below. Must be a value between 10 and 100. If null, no periodicity filtering is done"),
    # make_option(c("-s", "--start"), type="logical", help="use the translation initiation site as reference codon for calculating offsets. If 'FALSE', the second last codon is used instead"),
    # make_option(c("-e", "--extremity"), type="character", help="Specifies if the offset correction step should be based on 5' extremities ('5end') or 3'extremities ('3end') or 'auto' i.e. the optimal extremity is automatically selected"),
    # make_option(c("-n", "--start_nts"), type="integrer", help="Number of nt from start codon used to exclude P-sites near initiating ribosome when calculating CDS window P-site counts"),
    # make_option(c("-N", "--stop_nts"), type="integrer", help="Number of nucleotides from stop codon used to exclude P-sites near terminating ribosome when calculating CDS window P-site counts"),
    # make_option(c("-F", "--frequency_normalization"), type="logical", help="For codon usage index calculation"),
)
opt = parse_args(OptionParser(option_list=option_list))

folder      = opt$folder
output_dir  = opt$output_dir
region      = opt$region
#refCond     = opt$ref

files <- list.files(path = folder, pattern=".tsv", full.names=T)
riboseq_df = data.frame()
rnaseq_df = data.frame()
for (file in files) {
    if (endsWith(file, paste0(region, ".tsv"))){
        if (grepl("riboseq", file)){
            df <- read.table(file = file, header=TRUE, sep = "\t")
            # colnames(df)[3:ncol(df)] <- paste0(colnames(df)[3:ncol(df)], "_", sample)
            if (nrow(riboseq_df) != 0){
                riboseq_df <- merge(riboseq_df, df, by=c("id", "name"))
            }else {
                riboseq_df <- df
            }}
        if (grepl("rnaseq", file)){
            df <- read.table(file = file, header=TRUE, sep = "\t")
            # colnames(df)[3:ncol(df)] <- paste0(colnames(df)[3:ncol(df)], "_", sample)
            if (nrow(rnaseq_df) != 0){
                rnaseq_df <- merge(rnaseq_df, df, by=c("id", "name"))
            }else {
                rnaseq_df <- df
            }}
    }
}

rnaseq_df   <- tibble::column_to_rownames(rnaseq_df, var = "id")
riboseq_df  <- tibble::column_to_rownames(riboseq_df, var = "id")

rnaseq_df   <- subset(rnaseq_df, select=-c(name))
riboseq_df  <- subset(riboseq_df, select=-c(name))

print("Files:")
print(colnames(riboseq_df))
print("Using the following design for contrast:")
vec <- sapply(strsplit(colnames(riboseq_df), "_"), "[[", 1)
print(vec)

conditions <- unique(vec)

# Toutes les paires de conditions
pairs <- combn(conditions, 2, simplify = FALSE)

for (pair in pairs) {
    cond_a <- pair[1]
    cond_b <- pair[2]
    pair_label <- paste0(cond_a, "_vs_", cond_b)
    
    message("\n========== Contrast: ", pair_label, " ==========")
    
    # Sélectionner uniquement les colonnes des deux conditions
    sel_cols <- vec %in% c(cond_a, cond_b)
    ribo_sub <- riboseq_df[, sel_cols, drop = FALSE]
    rna_sub  <- rnaseq_df[,  sel_cols, drop = FALSE]
    vec_sub  <- vec[sel_cols]
    
    # Créer le répertoire de sortie propre à la paire
    pair_outdir <- file.path(output_dir, pair_label)
    dir.create(pair_outdir, recursive = TRUE, showWarnings = FALSE)
    setwd(pair_outdir)

    rep_by_cond <- table(vec_sub)
    message("Replicates per condition: ")
    print(rep_by_cond)

    if (all(rep_by_cond >= 3)) {
        message("Importing data into anota2seq...")
    } else {
        message("anota2seq requires 3 replicate experiments per group if there are 2 conditions")
        next
    } 
    ads <- anota2seqDataSetFromMatrix(
        dataP    = ribo_sub,
        dataT    = rna_sub,
        phenoVec = vec_sub,
        dataType = "RNAseq",
        normalize = TRUE)
    
    message("Running anota2seq...")
    if (all(rep_by_cond >= 3)) {
        ads <- anota2seqRun(ads)
    } else {
      ads <- anota2seqRun(ads, onlyGroup = TRUE)
    }
    
    message("Plotting results...")
    anota2seqPlotPvalues(ads, selContrast = 1, plotToFile = TRUE)
    anota2seqPlotFC(ads,      selContrast = 1, plotToFile = TRUE)
    
    message("Saving tables...")
    for (analysis in c("buffering", "translation", "mRNA abundance", "total mRNA")) {
        fname <- paste0("ANOTA2SEQ_", gsub(" ", "_", analysis), ".tsv")
        tbl <- anota2seqGetOutput(
            object      = ads,
            output      = "full",
            selContrast = 1,
            analysis    = analysis,
            getRVM      = TRUE)
        write.table(tbl, fname, sep = "\t", quote = FALSE, row.names = FALSE)
    }
    
    message("Done: ", pair_label)
}

message("\nAll pairwise contrasts completed.")