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
  Names_Col <- colnames(data)
  conditions <- gsub("\\..*","",Names_Col)
  
  nb_col_refCond <- table(factor(conditions[grepl(reference_condition, conditions)]))
  nb_col_conditions <- table(factor(conditions[!grepl(reference_condition, conditions)]))
  # Create a logical vector where TRUE corresponds to the reference condition
  is_ref_condition <- conditions == reference_condition
  
  # Order the columns so that reference condition columns come first
  ordered_cols <- c(
    Names_Col[grepl(reference_condition, conditions)],
    Names_Col[!grepl(reference_condition, conditions)])
  
  # Sort the data frame by the ordered columns
  expData_Sorted <- data[, match(ordered_cols, colnames(data))]
  
  sorted_data_list <- list(
    expData_Sorted = expData_Sorted,
    Names_Col = ordered_cols,
    conditions = conditions,
    nb_col_refCond = nb_col_refCond,
    nb_col_conditions = nb_col_conditions
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
      rep("yellow", sum(data_list$nb_col_refCond)),
      rep("blue", sum(data_list$nb_col_conditions))
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
make_PCA <- function(data, axis = 2, conds) {  
  transposed_data <- t(data)
  res.pca <- PCA(transposed_data, ncp = 3, graph = FALSE)
  
  if (dim(res.pca$eig)[1] > axis - 1) {    
    return(pca_plot(res.pca, axis, conds))
    }
}

pca_plot <- function(res.pca, axis, conds) {
  
  symmetric_limits <- function (x) {
    c(-max(abs(x)),  max(abs(x)))
  }

  percent_var_explained <- (res.pca$eig^2 / sum(res.pca$eig^2))*100
  x_axis = "Dim.1"
  percent_var_x = res.pca$eig[1,2]
  percent_var_y = res.pca$eig[axis,2]
  colours <- c("#be95c4", "#5fa8d3", "#a7c957", "#ffbd00", "#e63946", "#e76f51")
  
  y_axis = paste0("Dim.", axis)
  data = as.data.frame(res.pca$ind$coord)
  conditions = c()
  for (i in c(1:length(rownames(data)))) {
      for (cond in conds){
          if (grepl(cond, rownames(data)[i])){
              conditions <- c(conditions, cond)
              }
          }
      }
  data <- cbind(data, conditions)
  pca_plot <- ggplot(data, aes_string(x = x_axis, y = y_axis)) +
      geom_point(aes(fill=as.factor(data$conditions)), shape = 21, alpha = 0.7, size = 6, stroke=0.75) +
      theme_bw() +
      scale_fill_manual(aes(order=as.factor(data$conditions)), values=colours) +
      xlab(paste(x_axis, " (",round(percent_var_x, digit=2),"%)", sep=""))+
      ylab(paste(y_axis, " (",round(percent_var_y, digit=2),"%)", sep=""))+
      ggtitle("PCA graph of individuals")+
      theme(
        axis.text=element_text(size=16),
        legend.title= element_blank(),
        axis.title=element_text(size=16),
        legend.text = element_text(size =13),
        plot.title = element_text(size=18, face="bold", hjust = 0.5))+
        scale_x_continuous(limits = symmetric_limits) +
        scale_y_continuous(limits = symmetric_limits) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey")

  pca_plot <- pca_plot +
    ggrepel::geom_text_repel(
      aes(label = rownames(res.pca$ind$coord)),
      segment.alpha = 0.8,
      box.padding = 1,
      nudge_x = 1,
      force = 1,
      size=4,
      max.overlaps = Inf)

  return(pca_plot)
}

# DESeq2 matrix creation
deseq_object <- function(data_list) {
  conds <-
    factor(c(
      rep("CondA", sum(data_list$nb_col_refCond)),
      rep("CondB", sum(data_list$nb_col_conditions))
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
library_size_barplot <- function(dds_object, data_list, group, refCond) {
  normalized_counts <- counts(dds_object, normalized = TRUE)
  barplot(
    colSums(normalized_counts / 1000000),
    ylab = "Total read number (million)",
    main = "Library sizes (normalized)",
    sub = paste0(refCond, " vs ", group),
    col = c(
      rep("yellow", sum(data_list$nb_col_refCond)),
      rep("blue", sum(data_list$nb_col_conditions))
    ),
    names = colnames(data_list$expData_Sorted),
    las = 2,
    cex.names = 0.6
  )
}

# Library Size Boxplot
library_size_boxplot <- function(dds_object, data_list, group, refCond) {
  normalized_counts <- counts(dds_object, normalized = TRUE)
  row_sub_1 = apply(data_list$expData_Sorted, 1, function(row)
    all(row != 0))
  row_sub_2 = apply(normalized_counts, 1, function(row)
    all(row != 0))
  par(mfrow = c(1, 2), cex.axis=0.75, cex.lab=0.75)
  boxplot((data_list$expData_Sorted[row_sub_1, ]),
          main = "Library sizes",
          sub = paste0(refCond, " vs ", group),
          log = "y",
          col = c(
            rep("yellow", sum(data_list$nb_col_refCond)),
            rep("blue", sum(data_list$nb_col_conditions))
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
          sub = paste0(refCond, " vs ", group),
          log = "y",
          col = c(
            rep("yellow", sum(data_list$nb_col_refCond)),
            rep("blue", sum(data_list$nb_col_conditions))
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
dispersion_estimation <- function(ddsObjet, group, refCond) {
  ddsEstim = DESeq(ddsObjet)
  resDESeq = results(ddsEstim, contrast = c("condition", "CondB", "CondA"))
  
  mcols(resDESeq, use.names = TRUE)
  dds = estimateDispersions(ddsObjet)
  plotDispEsts(dds, main = paste0(refCond, " vs ", group), legend = "FALSE")
  
  return(resDESeq)
}

# Determine means
find_means <- function(dds_object, data_list) {
  normalized_counts <- counts(dds_object, normalized = TRUE)
  # Control mean
  if (sum(data_list$nb_col_refCond )> 1) {
    CT_Mean <-
      data.frame(apply(normalized_counts[, 1:(sum(data_list$nb_col_refCond))], 1, mean))
  } else {
    CT_Mean <- data.frame(normalized_counts[, 1])
  }
  # Tested condition mean
  if (sum(data_list$nb_col_conditions) > 1) {
    Mut_Mean <-
      data.frame(apply(normalized_counts[, (sum(data_list$nb_col_refCond) + 1):dim(normalized_counts)[2]], 1, mean))
  } else {
    Mut_Mean <-
      data.frame(normalized_counts[, (sum(data_list$nb_col_refCond) + 1)])
  }
  names(CT_Mean)[1] <- "CT_Mean"
  names(Mut_Mean)[1] <- "Mut_Mean"
  
  return(list(CT_Mean = CT_Mean, Mut_Mean = Mut_Mean))
}

# logFC frequency
frequency_histogram <- function(deseq_results, group, refCond) {
  hist(
    deseq_results$log2FoldChange,
    nclass = 200,
    xlab = "logFC values",
    main = "Frequency / LogFC",
    sub = paste0(refCond, " vs ", group)
  )
  abline(v = 2, col = "red")
  abline(v = -2, col = "red")
}

# p-value
pval_distribution <- function(deseq_results, group, refCond) {
  hist(deseq_results[, "pvalue"],
       nclass = 100,
       xlab = "p-values",
       main = "Histogram of p-values (DESeq2)",
       sub = paste0(refCond, " vs ", group))
}

# adjusted p-value
padj_distribution <- function(deseq_results, group, refCond) {
  hist(deseq_results[, "padj"],
       nclass = 100,
       xlab = "padj",
       main = "Histogram adjusted p-values",
       sub = paste0(refCond, " vs ", group))
}

ma_plot <- function(deseq_results, data_list, group, refCond) {
  data = data.frame(deseq_results)
  # thresholds
  fdr = Var_padj
  fc = Var_log2FC
  
  data$baseMean <- log2(data$baseMean +1) # log2 transformation of baseMean for plotting
  
  genenames <- rownames(data)
  
  # Significance levels (3 : NS, 2 : Down, 1 : Up : 2)
  sig <- rep(3, nrow(data))
  sig[which(data$padj <= fdr & data$log2FoldChange < fc & abs(data$log2FoldChange) >= log2(fc))] = 2
  sig[which(data$padj <= fdr & data$log2FoldChange > fc & abs(data$log2FoldChange) >= log2(fc))] = 1
  data <- data.frame(name = genenames, mean = data$baseMean, lfc = data$log2FoldChange, padj = data$padj, sig = sig)
  . <- NULL
  data$sig <- as.factor(data$sig)
  .lev <- levels(data$sig) %>% as.numeric()
  new.levels <- c("Up", "Down", "NS") %>% .[.lev]
  data$sig <- factor(data$sig, labels = new.levels)
  
  data <- data[order(data$padj), ]
  
  complete_data <- stats::na.omit(data)
  filtered_data <- subset(complete_data, padj <= fdr & name!="" & abs(lfc) >= log2(fc))

  # Sort by padj and log2fc in ascending order and select top 5
  best_up <- utils::head(filtered_data[order(filtered_data$lfc, filtered_data$padj, decreasing=c(TRUE, FALSE)), ], 3)
  # Sort by padj in ascending order and log2fc in descending order and select top 5
  best_down <- utils::head(filtered_data[order(filtered_data$lfc, filtered_data$padj,  decreasing=c(FALSE, FALSE)), ], 3)

  # Combine the rows
  labs_data <- rbind(best_up, best_down)
  mean <- lfc <- sig <- name <- padj <-  NULL

  # Plot creation
  ma_plot <- ggplot(data, aes(x = mean, y = lfc)) +
    geom_point(aes(alpha = sig, colour = sig, fill = sig, size = sig), shape=21) +
    # color of the border of points
    scale_colour_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) + 
    # size of points
    scale_size_manual(values= c("Up" = 1, "Down" = 1, "NS" = 0.5)) +
    # transparency of points
    scale_alpha_manual(values = c("Up" = 1, "Down" = 1, "NS" = 0.5)) +
    # fill color of points
    scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
    # adding labels for most significant element
    ggrepel::geom_label_repel(
      data = labs_data, aes(label = name),
      box.padding = unit(0.35, "lines"),
      segment.alpha = 0.8,
      fill = alpha(c("white"),0.4),
      point.padding = unit(0.3, "lines"),
      nudge_y = ifelse(labs_data$lfc > 0, 2, -2),
      min.segment.length = unit(0, 'lines'),
      force = 0.8, fontface = "plain",
      size = 2.5, color = "black",
      max.overlaps = Inf,
      show.legend = FALSE) + 
    # changing axis scale
    scale_x_continuous(breaks=seq(0, max(data$mean), 2)) +
    # adding titles 
    labs(
      # axis titles
      x = "Mean of normalized counts",
      y = "log2FC",
      # main title
      title = "MA-plot",
      subtitle = paste0(refCond," vs ", group),
      # legend title
      color = "Expression change") +
    # adding horizontal line(s) for log2FC threshold
    geom_hline(yintercept = c(0, -log2(fc), log2(fc)), linetype = c(1, 2, 2), color = c("black", "black", "black")) +
    # legend settings
    guides(
      color = FALSE,  # Remove the legend for color
      alpha = FALSE,  # Remove the legend for alpha
      size = FALSE,   # Remove the legend for size
      fill = guide_legend(override.aes = list(size=6)) # Incresing point size in legend
      ) + 
    # Remove secondary axis
    theme_classic() +
    theme(
      # graphical options for legend
      legend.position="bottom",
      legend.title = element_blank(),
      legend.text=element_text(size=14),
      plot.title = element_text(face = "bold", hjust=0.5, size=15),
      plot.subtitle = element_text(hjust=0.5),
      panel.border = element_rect(colour = "black", fill = NA, size= 0.5),    
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.text=element_text(size=12),
      axis.title=element_text(size=14,face="bold"))
  
  return(ma_plot)
}


# Volcano plot global
volcano_plot <- function(deseq_results, data_list, group, refCond) { # Volcano plot function
    
  data <- data.frame(deseq_results)
  fdr <- Var_padj
  fc <- Var_log2FC

  data$diffexpressed <- 
    ifelse( # Add a new column with the expression change
        abs(data$log2FoldChange) > fc & data$padj <= fdr, # If the gene is differentially expressed
            ifelse(data$log2FoldChange > fc, 
                "Up",
                "Down"), # If the gene is up-regulated or down-regulated
        "NS") # If the gene is not differentially expressed
        
  data$diffexpressed[is.na(data$diffexpressed)] <- "NS" # Replace the NA values with "NS"
  
  limit_x <- max(abs(ceiling(range(data$log2FoldChange[!is.na(data$log2FoldChange)])))) # Set the limit of the x-axis
  limit_y <- -log10(range(data$padj[!is.na(data$padj)])) # Set the limit of the y-axis

  up <- subset(data, diffexpressed == "Up") # Select the up-regulated genes
  down <- subset(data, diffexpressed == "Down") # Select the down-regulated genes

  up <- up[order(up$lfc, decreasing=TRUE), ] # Order the up-regulated genes by padj
  down <- down[order(down$lfc, decreasing=FALSE), ] # Order the down-regulated genes by padj

  best_up <- utils::head(up[order(up$padj, decreasing=FALSE), ], 3)
  best_down <- utils::head(down[order(down$padj, decreasing=FALSE), ], 3)
  top_significant_values <- rbind(best_up, best_down)
  top_significant_values$colour <- ifelse(top_significant_values$diffexpressed == "Up", "red", "blue")

  up_genes <- nrow(subset(data, log2FoldChange > fc & padj < fdr)) # Count the number of up-regulated genes
  down_genes <- nrow(subset(data, log2FoldChange < fc & padj < fdr)) # Count the number of down-regulated genes

  volcano_plot_global <- ggplot(data, aes(x = log2FoldChange, y = -log10(padj))) + # Create the volcano plot object
    geom_point(aes(colour = diffexpressed, alpha = diffexpressed, size = diffexpressed), shape = 16) +
    scale_colour_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey"), name = NULL) + 
    scale_alpha_manual(values = c("Up" = 1, "Down" = 1, "NS" = 0.5), name = NULL) +
    scale_size_manual(values = c("Up" = 2, "Down" = 2, "NS" = 1), name = NULL) +
    scale_fill_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey"), name = "Expression change") +
    # Add the most significant values
    geom_point(data = top_significant_values, shape = 21, size = 3, fill=top_significant_values$colour, colour = "black") +
    # Horizontal line for padj corresponding to padj threshold
    geom_hline(yintercept = -log10(fdr), linetype = "dashed") +
    # Vertical line for log2(fold change) corresponding to log2FC threshold
    geom_vline(xintercept = c(-fc, fc), linetype = "dashed") +
    # Add labels for the most significant values
    ggrepel::geom_label_repel(
        data = top_significant_values,
        aes(x = log2FoldChange, y = -log10(padj)), 
        size = 2.5,
        segment.alpha = 0.8,
        fill = alpha(c("white"),0.4),
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines"),
        fontface = "plain",
        label = row.names(top_significant_values),
        force = 0.8,
        min.segment.length = unit(0, 'lines'),
        max.overlaps = Inf,
        show.legend = FALSE
    ) +
    # Set the limits of the x-axis 
    scale_x_continuous(breaks = c(seq(-limit_x, limit_x, 5)), limits = c(-limit_x, limit_x)) +
    # Remove the legend for the color aesthetic
    guides(color = guide_legend(override.aes = list(size = 5))) +
    labs(
        title = "Volcano plot",
        subtitle = paste0(refCond," vs ", group),
        x = "log2(fold change)",
        y = "-log10(adjusted P-value)",
        colour = "Expression change"
    ) +
    geom_text(x = Inf, y = Inf, label = paste("Up:", up_genes, "\n Down:", down_genes),
        hjust = 1, vjust = 1.5, size = 4, color = "black") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 14),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 15),
      plot.subtitle = element_text(hjust=0.5),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5),    
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold")
    )

  return(volcano_plot_global)
}

# Zoomed Volcano plot function
volcano_plot_zoom <- function(deseq_results, data_list, group, refCond) {
  limite <- Var_padj/(10^(-log10(Var_padj)*2)) # Set the limit of the y-axis
  
  zoomed_in_volcano_plot <- volcano_plot(deseq_results, data_list, group, refCond) + # Create the zoomed volcano plot object
    scale_y_continuous(limits = c(0, -log10(limite))) + # Set the limits of the y-axis
    geom_point(data = data.frame(subset(data.frame(deseq_results), padj < limite)), aes(y=-log10(limite)), colour="red", alpha = 0.5, shape = 17, size = 3) +
    labs(title = "Volcano plot Zoom",
      subtitle = paste0(refCond," vs ", group))
  
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

main_function <- function(data, refCond, groups) {
  
  expData_sorted_gene_list <- c()
  for (i in c(1:length(groups))) {
    name <- groups[i]
    expData_gene <- data %>% select(contains(c(refCond, name)))
    expData_sorted_gene_list[[i]]  <- sort_columns_by_sample(expData_gene, refCond)
  }
  return(expData_sorted_gene_list)
}

deseq_list <- function(expData_sorted_gene_list) {
  deseq_object_list <- c()
  for (i in c(1:length(expData_sorted_gene_list))){
    dds <- deseq_object(expData_sorted_gene_list[[i]])
    deseq_object_list[[i]] <- dds
  }
  return(deseq_object_list)
}

dispersion_object_list <- function(dds_object_gene_list, groups, refCond) {
  dispersion_list <- c()
  for (i in c(1:length(dds_object_gene_list))){
    dispersion_list[[i]] <-  dispersion_estimation(dds_object_gene_list[[i]], groups[i], refCond)
  }
  return(dispersion_list)
}

mean_list <- function(dds_object_gene_list, expData_sorted_gene_list) {
  means_list <- c()
  for (i in c(1:length(expData_sorted_gene_list))){
    means <- find_means(dds_object_gene_list[[i]], expData_sorted_gene_list[[i]])
    means_list[[i]] <- means
  }
  return(means_list)
}

# Write tables
write_DE_tables <- function(tables_list, gene_transcript = "gene", refCond, group) {
  if (gene_transcript == "gene") {
    path_to_tables <- paths_list$DESeq2_gene
  } else {
    path_to_tables <- paths_list$DESeq2_transcript
  }
  write.table(
    tables_list$allGenes,
    file = paste0(path_to_tables, gene_transcript, "_complete_", refCond, "-", group, ".csv"),
    quote = F,
    sep = "\t",
    row.names = F
  )
  
  write.table(
    tables_list$inducedGenes,
    file = paste0(path_to_tables, gene_transcript, "_up_", refCond, "-", group, ".csv"),
    quote = F,
    sep = "\t",
    row.names = F
  )
  
  write.table(
    tables_list$repressedGenes,
    file = paste0(path_to_tables, gene_transcript, "_down_", refCond, "-", group, ".csv"),
    quote = F,
    sep = "\t",
    row.names = F
  )
}
