# Functions used for matrices creation and differential analysis

# Scan config file and determine paths
DESeq2_folder_paths <- function(path) {
  DESeq2_folder <- paste0(
    path,
    "RESULTS/DESeq2_",
    gsub(" ", "", params[which(params == "gff_cds_feature") +
                           1], fixed = TRUE),
    ".",
    gsub(" ", "", params[which(params == "readsLength_min") +
                           1], fixed = TRUE),
    "-",
    gsub(" ", "", params[which(params == "readsLength_max") +
                           1], fixed = TRUE),
    "/"
  )
  DESeq2_gene <- paste0(DESeq2_folder, "DESeq2_by_gene/")
  DESeq2_transcript <-
    paste0(DESeq2_folder, "DESeq2_by_transcript/")
  path_list <- list(
    DESeq2_folder = DESeq2_folder,
    DESeq2_gene = DESeq2_gene,
    DESeq2_transcript = DESeq2_transcript,
    pathway_matrix = paste0(DESeq2_folder, "count_matrix_by_transcript_IDs.csv"),
    pathway_names = paste0(DESeq2_folder, "names_correspondence_list.csv")
  )
}

# Merge counts by gene or by transcript
transcript_or_gene <- function(data) {
  expData_named_gene <- data[, -(dim(data)[2] - 1)]
  expData_counts_by_gene <-
    aggregate(expData_named_gene[, -dim(expData_named_gene)[2]], expData_named_gene["Gene_name"], sum)
  expData_counts_by_gene[, 2:(dim(expData_counts_by_gene)[2])] <-
    round(expData_counts_by_gene[, 2:(dim(expData_counts_by_gene)[2])])
  
  write.table(
    x = data.frame(expData_counts_by_gene, row.names = 1),
    file = paste0(paths_list$DESeq2_gene, "count_matrix_by_gene.csv"),
    sep = "\t",
    col.names = NA,
    row.names = T,
    quote = F
  )
  
  expData_named_transcript <- data[, -(dim(data)[2])]
  expData_counts_by_transcript <-
    aggregate(expData_named_transcript[, -dim(expData_named_transcript)[2]], expData_named_transcript["Transcript_name"], sum)
  expData_counts_by_transcript[, 2:(dim(expData_counts_by_transcript)[2])] <-
    round(expData_counts_by_transcript[, 2:(dim(expData_counts_by_transcript)[2])])
  
  write.table(
    x = data.frame(expData_counts_by_transcript, row.names = 1),
    file = paste0(
      paths_list$DESeq2_transcript,
      "count_matrix_by_transcript.csv"
    ),
    sep = "\t",
    col.names = NA,
    row.names = T,
    quote = F
  )
  
  # new_expData <- data.frame(expData_counts_by_transcript, row.names = 1)
  
  return(data.frame(expData_counts_by_transcript, row.names = 1))
}


# Determination of the number of samples
sort_columns_by_sample <- function(data, reference_condition) {
  expData_Sorted <- data[, sort(colnames(data))]
  Names_Col <- colnames(expData_Sorted)
  Cond <-
    str_replace(Names_Col, regex("[.]{0,1}[[:digit:]]{1,}$"), "")
  nb_col_condA <- table(factor(Cond))[[1]]
  nb_col_condB <- table(factor(Cond))[[2]]
  
  
  if (Cond[1] != reference_condition) {
    Ord <- c((1:nb_col_condB) + nb_col_condA , (1:nb_col_condA))
    expData_Sorted <- expData_Sorted[, Ord]
    
    Cond <- Cond[Ord]
    nb_col_condA <- nb_col_condA + nb_col_condB
    nb_col_condB <- nb_col_condA - nb_col_condB
    nb_col_condA <- nb_col_condA - nb_col_condB
    Names_Col <- colnames(expData_Sorted)
  }
  sorted_data_list <- list(
    expData_Sorted = expData_Sorted,
    Names_Col = Names_Col,
    Cond = Cond,
    nb_col_condA = nb_col_condA,
    nb_col_condB = nb_col_condB,
    sample_test = Cond[nb_col_condA + 1]
  )
  return(sorted_data_list)
}

# Library sizes
library_size_plot <- function(data_list) {
  barplot(
    colSums(data_list$expData_Sorted / 1000000),
    ylab = "Total read number (million)",
    main = "Library sizes",
    col = c(
      rep("yellow", data_list$nb_col_condA),
      rep("blue", data_list$nb_col_condB)
    ),
    names = colnames(data_list$expData_Sorted),
    cex.names = 0.6,
    las = 2
  )
}

# Pairwise Scatter plot
panel.cor <- function(x, y) {
  r <- round(cor(10 ^ x, 10 ^ y, method = "spearman"), digits = 4)
  txt <- paste0("R = ", r)
  text(2, 2, txt, cex = 1)
}

upper.panel <- function(x, y) {
  points((x),
         (y),
         pch = ".",
         cex = 1,
         col = rgb(0, 0, 0, 0.15))
}

# PCA
make_PCA <- function(data, axis = 2) {
  res.pca <- PCA(t(data) , ncp = 3, graph = FALSE)
  if (dim(res.pca$eig)[1] > axis - 1) {
    plot(
      res.pca,
      choix = "ind",
      autoLab = "yes",
      axes = c(1, axis),
      width = 3,
      height = 3
    )
  }
}

# DESeq2 matrix creation
deseq_object <- function(data_list) {
  conds <-
    factor(c(
      rep("CondA", data_list$nb_col_condA),
      rep("CondB", data_list$nb_col_condB)
    ))
  colData <- data.frame(condition = conds)
  ddsObjet <-
    DESeqDataSetFromMatrix(countData = data_list$expData_Sorted,
                           colData   = colData,
                           formula( ~ condition))
  ddsObjet <- estimateSizeFactors(ddsObjet)
  Size_Factors <- sizeFactors(ddsObjet)
  # Normalised data
  # normCountData <- counts(ddsObjet, normalized = TRUE)
  return(ddsObjet)
}

# Library Size Normalization
library_size_barplot <- function(dds_object, data_list) {
  normalized_counts <- counts(dds_object, normalized = TRUE)
  barplot(
    colSums(normalized_counts / 1000000),
    ylab = "Total read number (million)",
    main = "Library sizes (after normalization)",
    col = c(
      rep("yellow", data_list$nb_col_condA),
      rep("blue", data_list$nb_col_condB)
    ),
    names = colnames(data_list$expData_Sorted),
    las = 2,
    cex.names = 0.6
  )
}

# Library Size Boxplot
library_size_boxplot <- function(dds_object, data_list) {
  normalized_counts <- counts(dds_object, normalized = TRUE)
  row_sub_1 = apply(data_list$expData_Sorted, 1, function(row)
    all(row != 0))
  row_sub_2 = apply(normalized_counts, 1, function(row)
    all(row != 0))
  par(mfrow = c(1, 2))
  
  boxplot((data_list$expData_Sorted[row_sub_1, ]),
          main = "Library sizes",
          log = "y",
          col = c(
            rep("yellow", data_list$nb_col_condA),
            rep("blue", data_list$nb_col_condB)
          ),
          names = colnames(data_list$expData_Sorted),
          cex.names = 0.6,
          las = 2
          ,
          ylim = c(1, max(data_list$expData_Sorted))
  )
  boxplot((normalized_counts[row_sub_2, ]),
          ylab = "Total read number",
          main = "Library sizes (after normalization)",
          log = "y",
          col = c(
            rep("yellow", data_list$nb_col_condA),
            rep("blue", data_list$nb_col_condB)
          ),
          names = colnames(data_list$expData_Sorted),
          cex.names = 0.6,
          las = 2,
          yaxt = "n",
          ylim = c(1, max(data_list$expData_Sorted))
  )
  
  par(mfrow = c(1, 1))
}

# Dispersion
dispersion_estimation <- function(ddsObjet) {
  ddsEstim = DESeq(ddsObjet)
  resDESeq = results(ddsEstim, contrast = c("condition", "CondB", "CondA"))
  
  mcols(resDESeq, use.names = TRUE)
  dds = estimateDispersions(ddsObjet)
  plotDispEsts(dds)
  
  return(resDESeq)
}

# Determine means
find_means <- function(dds_object, data_list) {
  normalized_counts <- counts(dds_object, normalized = TRUE)
  # Control mean
  if (data_list$nb_col_condA > 1) {
    CT_Mean <-
      data.frame(apply(normalized_counts[, 1:(data_list$nb_col_condA)], 1, mean))
  } else {
    CT_Mean <- data.frame(normalized_counts[, 1])
  }
  # Tested condition mean
  if (data_list$nb_col_condB > 1) {
    Mut_Mean <-
      data.frame(apply(normalized_counts[, (data_list$nb_col_condA + 1):dim(normalized_counts)[2]], 1, mean))
  } else {
    Mut_Mean <-
      data.frame(normalized_counts[, (data_list$nb_col_condA + 1)])
  }
  
  return(list(CT_Mean = CT_Mean, Mut_Mean = Mut_Mean))
}

# logFC frequency
frequency_histogram <- function(deseq_results) {
  hist(
    deseq_results$log2FoldChange,
    nclass = 200,
    xlab = "logFC values",
    main = "Frequency / LogFC"
  )
  abline(v = 2, col = "red")
  abline(v = -2, col = "red")
}

# p-value
pval_distribution <- function(deseq_results) {
  hist(deseq_results[, "pvalue"],
       nclass = 100,
       xlab = "p-values",
       main = "Histogram of p-values (DESeq2)")
}

# adjusted p-value
padj_distribution <- function(deseq_results) {
  hist(deseq_results[, "padj"],
       nclass = 100,
       xlab = "padj",
       main = "Histogram adjusted p-values")
}

ma_plot <- function(deseq_results, data_list) { # MA-plot function
  data = data.frame(deseq_results)
  fdr = Var_padj
  fc = Var_log2FC
  
  detection_call = rep(1, nrow(data))
  
  data$baseMean <- log2(data$baseMean +1)
  genenames <- rownames(data)
  
  sig <- rep(3, nrow(data))
  sig[which(data$padj <= fdr & data$log2FoldChange < 0 & abs(data$log2FoldChange) >= log2(fc) & detection_call ==1)] = 2
  sig[which(data$padj <= fdr & data$log2FoldChange > 0 & abs(data$log2FoldChange) >= log2(fc) & detection_call ==1)] = 1
  data <- data.frame(name = genenames, mean = data$baseMean, lfc = data$log2FoldChange, padj = data$padj, sig = sig)
  
  . <- NULL
  data$sig <- as.factor(data$sig)
  .lev <- levels(data$sig) %>% as.numeric()
  new.levels <- c(paste0("Up: ", sum(sig == 1)), paste0("Down: ", sum(sig == 2)), "NS") %>% .[.lev]
  
  data$sig <- factor(data$sig, labels = new.levels)
  
  data <- data[order(data$padj), ]
  
  complete_data <- stats::na.omit(data)
  labs_data <- subset(complete_data, padj <= fdr & name!="" & abs(lfc) >= log2(fc))
  labs_data <- utils::head(labs_data, 5)
  
  # Plot
  mean <- lfc <- sig <- name <- padj <-  NULL
  
  ma_plot <- ggplot(data, aes(x = mean, y = lfc)) +
    geom_point(aes(alpha = sig, colour=sig, fill=sig, size=sig), shape=21) +
    scale_colour_manual(values = c("red", "red", "grey")) +
    scale_size_manual(values= c(1, 1, 0.5), guide = guide_legend(override.aes = list(size = 5))) +
    scale_alpha_manual(values = c(1, 1, 0.5)) +
    scale_fill_manual(values = c("red", "red", "grey")) +
    ggrepel::geom_label_repel(data = labs_data, mapping = aes(label = name),
                              box.padding = unit(0.35, "lines"),
                              point.padding = unit(0.3, "lines"),
                              force = 1, seed = 42, fontface = "plain",
                              size = 12/3, color = "black",
                              max.overlaps = Inf,
                              show.lgend = FALSE) + 
    scale_x_continuous(breaks=seq(0, max(data$mean), 2)) +
    labs(x = "Mean of normalized counts",
         y = "log2FC", 
         title = paste0("MA-plot - ",data_list$Cond[data_list$nb_col_condB+1]," vs ", data_list$Cond[1]), 
         color = "Expression change")+ # to remove legend title use color = ""
    geom_hline(yintercept = c(0, -log2(fc), log2(fc)), linetype = c(1, 2, 2),
               color = c("black", "black", "black")) + 
    guides(
      color = FALSE,  # Remove the legend for color
      alpha = FALSE,  # Remove the legend for alpha
      size = FALSE
    ) + 
    theme_classic() +
    theme(legend.position="bottom",
          plot.title = element_text(face = "bold", hjust=0.5, size=20),
          panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"),
          legend.title = element_blank(),
          legend.text=element_text(size=12))
  
  return(ma_plot)
}

# Volcano plot global
volcano_plot <- function(deseq_results, data_list) { # Volcano plot function
  deseq_results_modified <- data.frame(deseq_results) # Create a new dataframe with the results of DESeq2
  deseq_results_modified$diffexpressed <- ifelse( # Add a new column with the expression change
    abs(deseq_results_modified$log2FoldChange) > Var_log2FC & deseq_results_modified$padj < Var_padj,
    "Significant", # If the gene is differentially expressed
    "NS") # If the gene is not differentially expressed
  
  limit_x <- max(abs(ceiling(range(deseq_results_modified$log2FoldChange[!is.na(deseq_results_modified$log2FoldChange)])))) # Set the limit of the x-axis
  limit_y <- -log10(range(deseq_results_modified$padj[!is.na(deseq_results_modified$padj)])) # Set the limit of the y-axis
  
  top_significant_values <- deseq_results_modified %>% filter(padj %in% sort(deseq_results_modified$padj)[1:5]) # Select the 5 most significant values
  
  up_genes <- nrow(subset(deseq_results_modified, log2FoldChange > Var_log2FC & padj < Var_padj)) # Count the number of up-regulated genes
  down_genes <- nrow(subset(deseq_results_modified, log2FoldChange < Var_log2FC & padj < Var_padj)) # Count the number of down-regulated genes
  
  volcano_plot_global <- ggplot(deseq_results_modified, aes(x = log2FoldChange, y = -log10(padj))) + # Create the volcano plot object
    geom_point(aes(colour = diffexpressed, alpha = diffexpressed, size = diffexpressed), shape = 16) + # Add the values
    geom_point(data = top_significant_values, aes(colour = "Significant"), shape = 21, size = 3, fill = "red", colour = "black") + # Add the 5 most significant values
    geom_hline(yintercept = -log10(Var_padj), linetype = "dashed") + # Horizontal line at p-value = 0.05
    geom_vline(xintercept = c(-Var_log2FC, Var_log2FC), linetype = "dashed") + # Vertical line at log2FC = 1
    geom_label_repel(
      data = top_significant_values, # Add labels for the 5 most significant values
      aes(x = log2FoldChange, y = -log10(padj), label = row.names(top_significant_values)),
      force = 2,
      nudge_y = 1
    ) +
    scale_colour_manual(values = c("Significant" = "red", "NS" = "grey"), name =NULL) + 
    scale_alpha_manual(values = c("Significant" = 1, "NS" = 0.5), name =NULL) +
    scale_size_manual(values = c("Significant" = 2, "NS" = 1), name =NULL) +
    scale_fill_manual(values = c("Significant" = "red", "NS"="grey"), name="Expression change") + # If you want the fill color for significant points
    scale_x_continuous(breaks = c(seq(-limit_x, limit_x, 5)), limits = c(-limit_x, limit_x)) + # Set the limits of the x-axis 
    guides(
      color = guide_legend(override.aes = list(size = 5)) # Remove the legend for the alpha aesthetic
    ) +
    labs(
      title = paste0("Volcano plot - ", data_list$Cond[data_list$nb_col_condB + 1], " vs ", data_list$Cond[1]),
      x = "log2(fold change)",
      y = "-log10(adjusted P-value)",
      colour = "Expression change"
    ) +
    geom_text(x = Inf, y = Inf, label = paste("Up:", up_genes, "\n Down:", down_genes), hjust = 1, vjust = 1.5, size = 4, color = "black") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", hjust = 0.5, size = 20),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5),    
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    )
  
  return(volcano_plot_global)
}

# Zoomed Volcano plot function
volcano_plot_zoom <- function(deseq_results, data_list) {
  limite <- Var_padj/(10^(-log10(Var_padj)*2)) # Set the limit of the y-axis
  
  zoomed_in_volcano_plot <- volcano_plot(deseq_results, data_list) + # Create the zoomed volcano plot object
    scale_y_continuous(limits = c(0, -log10(limite))) + # Set the limits of the y-axis
    geom_point(data = data.frame(subset(data.frame(deseq_results), padj < limite)), aes(y=-log10(limite)), colour="red", alpha = 0.5, shape = 17, size = 3) +
    labs(title = paste0("Volcano plot Zoom - ",data_list$Cond[data_list$nb_col_condB+1]," vs ",data_list$Cond[1]))
  
  return(zoomed_in_volcano_plot)
}

# Final tables creation
tables_creation <-
  function(dds_object,
           deseq_results,
           data_list,
           means,
           reference_condition,
           gene_transcript = "gene") {
  
  normalized_counts <- counts(dds_object, normalized = TRUE)
  Bruts_Norm <- cbind(data_list$expData_Sorted, normalized_counts)
  Names_Col <- colnames(data_list$expData_Sorted)
  
  Table_Complete <-
    data.frame(cbind(
      row.names(Bruts_Norm) ,
      Bruts_Norm,
      means$CT_Mean,
      means$Mut_Mean ,
      deseq_results
    ))
  
  colnames(Table_Complete) =  c(
    "ID" ,
    Names_Col,
    paste("norm", data_list$Names_Col, sep = "_"),
    paste0(reference_condition, "_Mean"),
    paste0(data_list$sample_test, "_Mean"),
    colnames(deseq_results)
  )
  
  names_list = read.table(
    paths_list$pathway_names,
    header = T,
    row.names = 1,
    check.names = FALSE,
    sep = "\t"
  )
  if (gene_transcript == "gene" & ncol(names_list) > 0)
  {
    allGenes <- Table_Complete
  } else {
    allGenes <- merge(Table_Complete, names_list, by = 1)
  }
  
  inducedGenes = allGenes[which((allGenes[, "log2FoldChange"] > Var_log2FC) &
                                  (allGenes[, "padj"] < Var_padj)), ]
  dim(inducedGenes)
  
  repressedGenes = allGenes[which((allGenes[, "log2FoldChange"] < -Var_log2FC) &
                                    (allGenes[, "padj"] < Var_padj)), ]
  dim(repressedGenes)
  
  # Sort tables by padj
  inducedGenes <- inducedGenes[(order(inducedGenes[, "padj"])), ]
  repressedGenes <-
    repressedGenes[(order(repressedGenes[, "padj"])), ]
  
  return(
    list(
      allGenes = allGenes,
      inducedGenes = inducedGenes,
      repressedGenes = repressedGenes
    )
  )
}

# Write tables
write_DE_tables <- function(tables_list, gene_transcript = "gene") {
  if (gene_transcript == "gene") {
    path_to_tables <- paths_list$DESeq2_gene
  } else {
    path_to_tables <- paths_list$DESeq2_transcript
  }
  write.table(
    tables_list$allGenes,
    file = paste0(path_to_tables, gene_transcript, "_complete.csv"),
    quote = F,
    sep = "\t",
    row.names = F
  )
  
  write.table(
    tables_list$inducedGenes,
    file = paste0(path_to_tables, gene_transcript, "_up.csv"),
    quote = F,
    sep = "\t",
    row.names = F
  )
  
  write.table(
    tables_list$repressedGenes,
    file = paste0(path_to_tables, gene_transcript, "_down.csv"),
    quote = F,
    sep = "\t",
    row.names = F
  )
}
