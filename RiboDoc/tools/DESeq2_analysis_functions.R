

# Scan config file and determine paths
DESeq2_folder_paths <- function(path) {
  DESeq2_folder <- paste0(path, "RESULTS/DESeq2.",
                          gsub(" ", "", params[which(params=="readsLength_min")+1], fixed = TRUE), "-",
                          gsub(" ", "", params[which(params=="readsLength_max")+1], fixed = TRUE), "/")
  DESeq2_gene <- paste0(DESeq2_folder, "DESeq2_by_gene/")
  DESeq2_transcript <- paste0(DESeq2_folder, "DESeq2_by_transcript/")
  path_list <- list(DESeq2_folder = DESeq2_folder,
                    DESeq2_gene = DESeq2_gene,
                    DESeq2_transcript = DESeq2_transcript,
                    pathway_matrix = paste0(DESeq2_folder,"count_matrix_by_transcript_unnamed.txt"),
                    pathway_names = paste0(DESeq2_folder,"names_correspondence_list.txt"))
}

# Merge by gene
transcript_or_gene <- function(data, list_of_names, gene_transcript="gene") {
  expData_named <- data.frame(merge(data, list_of_names, by = "row.names"), row.names = 1)
  if(gene_transcript == "gene") {
    expData_named_gene <- expData_named[,-(dim(expData_named)[2]-1)]
    expData_counts_by_gene <- aggregate(expData_named_gene[,-dim(expData_named_gene)[2]], expData_named_gene["Gene_name"], sum)
    new_expData <- data.frame(expData_counts_by_gene, row.names = 1)
    write.table(new_expData, paste0(paths_list$DESeq2_gene,"count_matrix_by_gene.txt"), sep = "\t", col.names = NA, row.names = T, quote = F)
  } else {
    expData_named_transcript <- expData_named[,-(dim(expData_named)[2])]
    expData_counts_by_transcript <- aggregate(expData_named_transcript[,-dim(expData_named_transcript)[2]], expData_named_transcript["Transcript_name"], sum)
    new_expData <- data.frame(expData_counts_by_transcript, row.names = 1)
    write.table(new_expData, paste0(paths_list$DESeq2_transcript,"count_matrix_by_transcript.txt"), sep = "\t", col.names = NA, row.names = T, quote = F)
  }
  return(new_expData)
}

# Determination of the number of samples
sort_columns_by_sample <- function(data,reference_condition) {
  expData_Sorted <- data[,sort(colnames(data))]
  Names_Col <- colnames(expData_Sorted)
  Cond <- str_replace(Names_Col,regex("[.]{0,1}[[:digit:]]{1,}$"),"")
  nb_col_condA <- table(factor(Cond))[[1]]
  nb_col_condB <- table(factor(Cond))[[2]]
  sample_test <- Cond[nb_col_condA+1]
  
  
  if (Cond[1] != reference_condition) {
    Ord <- c( (1:nb_col_condB)+nb_col_condA , (1:nb_col_condA) )
    expData_Sorted <- expData_Sorted[, Ord]
    
    Cond <- Cond[Ord]
    nb_col_condA <- nb_col_condA + nb_col_condB
    nb_col_condB <- nb_col_condA - nb_col_condB
    nb_col_condA <- nb_col_condA - nb_col_condB
    Names_Col <- colnames(expData_Sorted)
    sample_test <- Cond[1]
  }
  sorted_data_list <- list(expData_Sorted = expData_Sorted,
                           Names_Col = Names_Col,
                           Cond = Cond,
                           nb_col_condA = nb_col_condA,nb_col_condB = nb_col_condB,
                           sample_test = sample_test)
  return(sorted_data_list)
}

# Library sizes
library_size_plot <- function(data_list) {
  barplot(colSums(data_list$expData_Sorted/1000000), ylab = "Total read number (million)",
          main = "Library sizes", col = c(rep("yellow", data_list$nb_col_condA), rep("blue", data_list$nb_col_condB)),
          names = colnames(data_list$expData_Sorted),
          cex.names = 0.6,
          las = 2
  )
}

# Pairwise Scatter plot
panel.cor <- function(x, y){
  r <- round(cor(10^x, 10^y, method="spearman"), digits=4 )
  txt <- paste0("R = ", r)
  text(2, 2, txt, cex = 1)
}

upper.panel<-function(x, y){
  points((x),(y), pch = ".", cex =1, col = rgb(0, 0, 0, 0.15))
}

# PCA
make_PCA <- function(data) {
  par(mfrow = c(2,1))
  res.pca <- PCA(t(data) , ncp = 3, graph = FALSE)
  
  plot(res.pca, choix ="ind", autoLab = "yes", axes = c(1,2), width=3, height=3)  
  
  plot(res.pca, choix ="ind", autoLab = "yes", axes = c(1,3))
  par(mfrow = c(1,1))
}

# DESeq2 matrix creation
deseq_object <- function(data_list) {
  conds <- factor(c(rep("CondA", data_list$nb_col_condA), rep("CondB", data_list$nb_col_condB)))
  colData <- data.frame(condition = conds)
  ddsObjet <- DESeqDataSetFromMatrix(countData = data_list$expData_Sorted,
                                     colData   = colData, formula(~ condition))
  ddsObjet <- estimateSizeFactors(ddsObjet)
  Size_Factors <- sizeFactors(ddsObjet)
  # Normalised data
  # normCountData <- counts(ddsObjet, normalized = TRUE)
  return(ddsObjet)
}

# Library Size Normalization
library_size_barplot <- function(dds_object, data_list) {
  normalized_couts <- counts(dds_object, normalized = TRUE)
  barplot(colSums(normalized_couts/1000000), ylab = "Total read number (million)",
          main = "Library sizes (after normalization)",
          col = c(rep("yellow", data_list$nb_col_condA), rep("blue", data_list$nb_col_condB)),
          names = colnames(data_list$expData_Sorted),
          las = 2,
          cex.names = 0.6)
}

# Library Size Boxplot
library_size_boxplot <- function(dds_object, data_list) {
  normalized_couts <- counts(dds_object, normalized = TRUE)
  row_sub_1 = apply(data_list$expData_Sorted, 1, function(row) all(row != 0 ))
  row_sub_2 = apply(normalized_couts, 1, function(row) all(row != 0 ))
  par(mfrow = c(1,2))
  
  boxplot((data_list$expData_Sorted[row_sub_1,]),
          main = "Library sizes",
          log = "y",
          col = c(rep("yellow", data_list$nb_col_condA), rep("blue", data_list$nb_col_condB)),
          names = colnames(data_list$expData_Sorted),
          cex.names = 0.6,
          las =2
          ,ylim = c(1,max(data_list$expData_Sorted) )
  )
  boxplot((normalized_couts[row_sub_2,]), ylab = "Total read number",
          main = "Library sizes (after normalization)",
          log = "y",
          col = c(rep("yellow", data_list$nb_col_condA), rep("blue", data_list$nb_col_condB)),
          names = colnames(data_list$expData_Sorted),
          cex.names = 0.6,
          las =2,
          yaxt = "n",
          ylim = c(1,max(data_list$expData_Sorted) )
  )
  
  par(mfrow = c(1,1))
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
find_means <- function(dds_object,data_list) {
  normalized_couts <- counts(dds_object, normalized = TRUE)
  # Control mean
  if(data_list$nb_col_condA > 1) {
    CT_Mean <- data.frame(apply(normalized_couts[,1:(data_list$nb_col_condA)],1,mean))
  } else {
    CT_Mean <- data.frame(normalized_couts[,1])
  }
  # Tested condition mean
  if(data_list$nb_col_condB > 1) {
    Mut_Mean <- data.frame(apply(normalized_couts[,(data_list$nb_col_condA+1):dim(normalized_couts)[2]],1,mean))
  } else {
    Mut_Mean <- data.frame(normalized_couts[,(data_list$nb_col_condA+1)])
  }
  
  return(list(CT_Mean = CT_Mean, Mut_Mean = Mut_Mean))
}

# logFC frequency
frequency_histogram <- function(deseq_results) {
  hist(deseq_results$log2FoldChange,
       nclass = 200,
       xlab = "logFC values",
       main ="Frequency / LogFC" )
  abline(v = 2, col = "red")
  abline(v = -2, col = "red")
}

# p-value
pval_distribution <- function(deseq_results) {
  hist(deseq_results[,"pvalue"], nclass = 100, xlab = "p-values",
       main = "Histogram of p-values (DESeq2)")
}

# adjusted p-value
padj_distribution <- function(deseq_results) {
  hist(deseq_results[,"padj"], nclass = 100, xlab = "padj",
       main = "Histogram adjusted p-values")
}

# MA-plot
ma_plot <- function(deseq_results) {
  deseq_results_modified <- deseq_results
  deseq_results_modified$padj[which(abs(deseq_results_modified$log2FoldChange) <= abs(Var_log2FC))] <- 1
  plotMA(deseq_results_modified, alpha = Var_padj, ylab = "log2FC", colSig = rgb(1,0,0,0.5))
}

# Volcano plot global
volcano_plot <- function(deseq_results, data_list) {
  myColors <- ifelse((deseq_results$padj < Var_padj & deseq_results$log2FoldChange > Var_log2FC), rgb(1, 0, 0, 0.25) ,
                     ifelse((deseq_results$padj < Var_padj & deseq_results$log2FoldChange < -Var_log2FC), rgb(1, 0, 0, 0.25) ,
                            rgb(0, 0, 0, 0.25) ) )
  
  plot(deseq_results[, "log2FoldChange"], -log10(deseq_results[, "padj"]),
       pch = 20, cex = 1,
       col=myColors,
       xlab = "log2 FC", ylab = "-log10(padj)",
       main = paste0("Volcano plot - ",data_list$Cond[data_list$nb_col_condB+1]," vs ", data_list$Cond[1])
  )
}

# ZOOM Volcano plot
volcano_plot_zoom <- function(deseq_results, data_list) {
  limite <- Var_padj/(10^(-log10(Var_padj)*2))
  
  
  myColors <- ifelse((deseq_results$padj < Var_padj & deseq_results$log2FoldChange > Var_log2FC), rgb(1, 0, 0, 0.25) ,
                     ifelse((deseq_results$padj < Var_padj & deseq_results$log2FoldChange < -Var_log2FC), rgb(1, 0, 0, 0.25) ,
                            rgb(0, 0, 0, 0.25) ) )
  
  Valeurs_Limite_Y <- ifelse((deseq_results$padj < limite), limite ,
                             deseq_results$padj) 
  
  Forme_Pixel_Y <- ifelse((deseq_results$padj < limite), 17 ,
                          20) 
  
  
  plot(deseq_results[, "log2FoldChange"], -log10(Valeurs_Limite_Y),
       pch = Forme_Pixel_Y,
       col=myColors,
       xlab = "log2 FC", ylab = "-log10(padj)",
       ylim = c(0,-log10(limite)),
       main = paste0("Volcano plot Zoom - ",data_list$Cond[data_list$nb_col_condB+1]," vs ",data_list$Cond[1])
  )
}

# Final tables creation
tables_creation <- function(dds_object, deseq_results, data_list, means, reference_condition, gene_transcript = "gene") {
  normalized_couts <- counts(dds_object, normalized = TRUE)
  Bruts_Norm <- cbind(data_list$expData_Sorted, normalized_couts)
  Names_Col <- colnames(data_list$expData_Sorted)
  
  Table_Complete <- data.frame(cbind(row.names(Bruts_Norm) , Bruts_Norm, means$CT_Mean, means$Mut_Mean , deseq_results) )
  
  
  colnames(Table_Complete) =  c("ID" ,
                                Names_Col,
                                paste("norm", data_list$Names_Col, sep = "_"),
                                paste0(reference_condition,"_Mean"),
                                paste0(data_list$sample_test, "_Mean"),
                                colnames(deseq_results))
  
  if(gene_transcript == "gene" & ncol(names_list) > 0)
  {
    allGenes <- Table_Complete
  } else {
    allGenes <- merge(Table_Complete, names_list, by = 1)
  }
  
  inducedGenes = allGenes[which((allGenes[, "log2FoldChange"] > Var_log2FC) & (allGenes[, "padj"] < Var_padj)),]
  dim(inducedGenes)
  
  repressedGenes = allGenes[which((allGenes[, "log2FoldChange"] < -Var_log2FC) & (allGenes[, "padj"] < Var_padj)),]
  dim(repressedGenes)
  
  # Sort tables by padj
  inducedGenes <- inducedGenes[(order(inducedGenes[,"padj"])),]
  repressedGenes <- repressedGenes[(order(repressedGenes[,"padj"])),]
  
  return(list(allGenes = allGenes, inducedGenes = inducedGenes, repressedGenes = repressedGenes))
}

# Write tables
write_DE_tables <- function(tables_list, gene_transcript = "gene") {
  if(gene_transcript == "gene") {
    path_to_tables <- paths_list$DESeq2_gene
  } else {
    path_to_tables <- paths_list$DESeq2_transcript
  }
  write.table(tables_list$allGenes,
              file = paste0(path_to_tables,gene_transcript,"_complete.txt"),
              quote = F, sep = "\t", row.names = F)
  
  
  write.table(tables_list$inducedGenes,
              file = paste0(path_to_tables,gene_transcript,"_up.txt"),
              quote = F, sep = "\t", row.names = F)
  
  write.table(tables_list$repressedGenes,
              file = paste0(path_to_tables,gene_transcript,"_down.txt"),
              quote = F, sep = "\t", row.names = F)
}

