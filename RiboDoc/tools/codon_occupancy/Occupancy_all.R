T1<-Sys.time()

#####################
# Libraries
#####################
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggsignif))
suppressPackageStartupMessages(library(ggh4x))
suppressPackageStartupMessages(library(ggtext))

#####################
# Parameters
#####################
option_list = list(
  
  make_option(c("-O", "--Offset"),
              default = "psite_offset.csv",
              type = "character",
              help = "File containing P-site offsets. Default : 'psite_offset.csv"),

  make_option(c("-F", "--Folder"),
              default = "./",
              type = 'character',
              help = "Path to designed R functions"),
  
  make_option(c("-N", "--Name_ref"),
              default = "WT",
              type = 'character',
              help = "Name reference. Default : 'WT'"),
  
  make_option(c("-n", "--name_test"),
              default = "Mut",
              type = 'character',
              help = "Name mutant. Default : 'Mut'"),
  
  make_option(c("-r", "--reads_occurence_threshold"),
              default = 0,
              type = "integer",
              help = "Minimum number of reads taken into account. Default : 0"),
  
  make_option(c("-m", "--minimun_length"),
              default = 25,
              type = "integer",
              help = "Minimum length of reads taken into account. Default : 25"),
  
  make_option(c("-M", "--Maximum_length"),
              default = 35,
              type = "integer",
              help = "Maximum length of reads taken into account. Default : 35"),

  make_option(c("-e", "--elongation"),
              default = 50,
              type = "integer",
              help = "Length of sequence elongation for transcriptome creation. Default : 50"),
  
  make_option(c("-s", "--site"),
              default = "A",
              type = "character",
              help = "Ribosomal decoding site to analyze ('A' or 'P'). Default : 'A'"),
  
  make_option(c("-f", "--frame"),
              default = "TRUE",
              type = "logical",
              help = "Select only reads with decoding site in a reading frame. Default : 'TRUE'"),
  
  make_option(c("-p", "--pathway"),
              default = "sequenceBedCount/",
              type = 'character',
              help = "Pathway to the sequenceBedCount folder. Default : 'sequenceBedCount/'"),
  
  make_option(c("-o", "--outpathway"),
              default = "codon_occupancy/",
              type = 'character',
              help = "Pathway to the output folder. Default : 'codon_occupancy/'")
)
opt = parse_args(OptionParser(option_list=option_list))

file_offset <- opt$O

# Functions to check if given pathways are ending with a '/' when it needs to be a folder
function_folder <- opt$F
pathway_functions <- paste0(function_folder,"check_pathway.R")
triplets_pathway <- paste0(function_folder,"codons_and_aa.csv")

name_ref <- opt$N
name_test <- opt$n
print(paste0("Name of the reference samples : ", name_ref))
print(paste0("Name of the tested samples : ", name_test))

minimal_read <- opt$r
print(paste0("Minimal coverage threshold : ", minimal_read))

min_length <- opt$m
print(paste0("Minimal reads length : ", min_length))
max_length <- opt$M
print(paste0("Maximal reads length : ", max_length))

elongation <- opt$e
print(paste0("Transcripts elongation length  : ", elongation))

site <- opt$s
print(paste0("Decoding site of the ribosome : ", site))

framed <- opt$f
print(paste0("Selection of reads in reading frames only : ", framed))

pathway_file_sequenceBedCount <- opt$p
pathway_file_graphs <- opt$o


#####################
# Functions
#####################

source(file = pathway_functions)

# Make a table with the mean counts from every read length in each condition
means_in_df <- function(data, mean_names) {
  mean_df <- NULL
  for (name in 1:length(mean_names)) {
    mean_df <-
      cbind(mean_df, rowMeans(data[, grep(paste0(mean_names[name]), colnames(data))]))
    colnames(mean_df)[dim(mean_df)[2]] <- mean_names[name]
  }
  return(mean_df)
}

# Remove useless info from data to only keep sequence and coordinates of 1 site
adjust_to_offset <- function(data, offsets, site_of_interest = "A") {
    # Load corresponding offset
    offset <- offsets[offsets$length == kmer, 7]
    # To get the sequence of the A-site only
    if (site_of_interest == "A") {
      offset <- offset + 3
    }
    # Adjust data to only keep sequence and coordinates of the site of interest
    data[, 3] <- data[, 2] + offset + 2
    data[, 2] <- data[, 2] + offset
    data[, 5] <- str_sub(data[, 5], offset + 1, offset + 3)
    
    return(data)
  }

# Select reads in readings frames only or take every read in the dataset
only_reading_frame <- function(data, in_frame, elong) {
  if(in_frame == T) {
    data_framed <- data[which(data$Five_prime_pos %% 3 == elong %% 3),]
    return(data_framed)
  } else {
    return(data)
  }
}

# Normalization by the total of reads (library size)
normalization <- function(counts) {
  cat("Normalizing counts... ")
  return(counts / sum(counts))
}

# Filter by minimal number of reads at a specific position
filtering <- function(data, min_nbr = 0) {
  cat(
    paste0(
      "Subsetting data to only keep positions with at least ",
      minimal_read,
      "reads... "
    )
  )
  if (min_nbr > 0) {
    data_filtered <- subset(data, Counts >= min_nbr)
  } else {
    data_filtered <- data
  }
  return(data_filtered)
}

# Filter then normalize or normalize then filter data
filt_norm <- function(data, min_nbr, filter_first) {
  data_filt_norm <- data
  if (filter_first) {
    data_filt_norm <-
      filtering(data = data_filt_norm, min_nbr = min_nbr)
    data_filt_norm$Counts <-
      normalization(counts = data_filt_norm$Counts)
  } else {
    data_filt_norm$Counts <-
      normalization(counts = data_filt_norm$Counts)
    data_filt_norm < filtering(data = data_filt_norm,
                               min_nbr = min_nbr / sum(data_filt_norm$Counts))
  }
  return(data_filt_norm)
}

# Remove triplets containing "N"
remove_Ns <- function(data) {
  return(data[grep("N", data$Codons, invert = T),])
}

# Aggregate counts by codon
aggregate_codons <- function(data) {
  return(aggregate(Counts ~ Codons, data = data, sum))
}

# In case not all codons are present in the BED file
missing_triplets <- function(data, counts, codons_table) {
  codons <- row.names(codons_table)
  if (dim(counts)[1] < length(codons)) {
    missing_triplets <- codons[!(codons %in% data$Codons)]
    
    # Add missing codons with 0 counts in table
    complete_counts <-
      rbind(
        aggregate(Counts ~ Codons, data = data, sum),
        data.frame(Codons = missing_triplets, Counts = 0)
      )
    complete_counts <-
      complete_counts[match(row.names(codons_table), complete_counts$Codons),]
    
  } else {
    complete_counts <-
      counts[match(row.names(codons_table), counts$Codons),]
  }
  return(complete_counts)
}

# General function for normalization, filtering and counting of each codon
codon_counts <- function(data, min_nbr, codons_table, filter_first, in_frame, elong) {

    framed_data <- only_reading_frame(data = data, in_frame = in_frame, elong = elong)
    
    data_filt_norm <- filt_norm(data = framed_data, min_nbr, filter_first = filter_first)
    
    data_filt_norm <- remove_Ns(data = data_filt_norm)
    
    codon_counts <- aggregate_codons(data = data_filt_norm)
    
    complete_counts <- missing_triplets(data = data_filt_norm, counts = codon_counts, codons_table = df_triplets)
    
    return(complete_counts)
  }

 createRectAnnotation <- function(data, facet = NULL, col_y = NULL, col_x, cond = NULL, x_min = NULL, x_max = NULL, y_min = 0, xmin_adjust = 1, xmax_adjust = xmin_adjust, ymin_adjust = 0, ymax_adjust = 0.005) {
  max_y = max(data[[col_y]])
  if(is.null(x_min)) {
    if(!is.null(facet)) {
      x_min = min(which(levels(data[[col_x]]) %like% cond))
      x_max = max(which(levels(data[[col_x]]) %like% cond))
    } else {x_min = 1}
  }
  subset_data <- data[data[[col_x]] %like% cond, ]
  if(is.null(x_max)) x_max <- x_min
  if(is.null(col_y)) {
    y_max <- +Inf
  } else {
    values <- subset_data[[col_y]]
    y_max <- values[which.max( abs(values) )]
  }
  return(geom_rect(
    data = subset_data,
    aes(x = NULL,
      y = NULL,
      xmin = x_min - xmin_adjust, 
      xmax = x_max + xmax_adjust, 
      ymin = y_min - ymin_adjust,
      ymax = y_max + ymax_adjust),
    alpha = 0,
    linewidth = 1,
    linetype = "dashed",
    color = "red",
    show.legend=FALSE))
}

pvalues_calc <- function(data, col) {

    data <- data.frame(data)
    p_values <- c()
    for (element in col) {
        t <- t.test(data[col == element, grep(name_ref, names(data))], data[col == element, grep(name_test, names(data))])
        p_values <- c(p_values, t$p.value)
    }
    return(p_values)
}

map_pvalue_to_label <- function(p_value) {

  if (p_value < 0.05) {
    if (p_value < 0.01) {
      if (p_value < 0.001) {
        return("***")
      } else {
        return("**")
      }
    } else {
      return("*")
    }
  } 
  else {
    return("")
  }
}

get_pos <- function(df){

  pos = 1
  list_pos = c(pos)
  for (i in c(2:length(df))) {
    if (df[i] != df[i-1])
      {pos = 1}
    else{
      pos = pos + 1}
    list_pos <- c(list_pos, pos)}
    
  return(list_pos)
  
}

# Process data : put rownames as a column with the name "Codon"
row_as_col <- function(data, col) {
  data <- setDT(data.frame(data), keep.rownames = TRUE)[]
  colnames(data)[1] <- col
  return(data)
}

#####################
# Execution
#####################
# Store all customizations paramters of plot in a list
customPlot <- list(
    theme_classic(),
    theme(
      legend.position = "top",
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      axis.title.x = element_blank(),
      axis.title.y =  element_text(size = 16),
      axis.text.x = element_markdown(size = 12, vjust = 0.5, hjust=1, angle = 90),
      axis.text.y = element_text(size = 14),
      axis.line.x = element_blank(),
      strip.placement = 'outside',
      strip.background.x = element_blank(),
      strip.text.x = element_text(size = 11, face="bold"),
      panel.spacing.y = unit(0, "lines"),
      strip.clip = "off"
    )
  )

codon_start <- "ATG"
codons_stop <- c("TAA", "TAG", "TGA")
aa_start <- "MET"
aa_stop <- "STOP"

colors = c("red", rep("black", 19), "red")

kmers <- seq(min_length, max_length)
cat("List of read lengths in window :\n")
cat(paste0("\t",kmers))
cat("\n")

# Check if pathways end with a "/" and create output folder
pathway_file_sequenceBedCount <- pathway_handling(pathway_file_sequenceBedCount)
pathway_file_graphs <- paste0(pathway_handling(pathway_file_graphs),site,"site/")
dir.create(pathway_file_graphs, recursive = T, showWarnings = F)

# Creation of the 64 codons
# bases <- c("A","T","G","C")
# df_triplets <- expand.grid(bases,bases,bases)
# triplets <- sort(as.vector(t(paste0(df_triplets[,1], df_triplets[,2],df_triplets[,3]))))
df_triplets <- data.frame(read.table(triplets_pathway, header = T), row.names = 1)
df_triplets <- df_triplets[order(df_triplets$Family),]
triplets <- row.names(df_triplets)

# Names strains determination via offset.csv file
fo <- read.table(file_offset, header = T)
fo[,9] <- str_replace(fo[,9], "transcriptome_elongated.","")
#doublons <- which(duplicated(fo[,9]))
#names_strains <- fo[,9][-doublons]

names_strains <- levels(factor(fo[,9]))
# Sort names strains with ref in first position
names_strains_sorted <- str_sort(names_strains, decreasing = F)

# Determination number of ref strains and mutants
occurence_ref <- str_count(names_strains_sorted, paste0("^",name_ref,".[0-9]+$"))
nb_sample <- length(occurence_ref)
number_WT <- sum(occurence_ref)
print(paste0("Number of reference samples : ", number_WT))
number_mut <- length(occurence_ref) - number_WT
print(paste0("Number of tested samples : ", number_mut))

# Search if ref names is in first position and change it if it is not the case
if (occurence_ref[1] == 0) {
  names_strains_sorted <- str_sort(names_strains, decreasing = T)
  
  # Modify numeric order of strains
  new_order <- c(seq(number_WT,1),seq(nb_sample, sum(occurence_ref)+1))
  names_strains_sorted <- names_strains_sorted[new_order]
  
  occurence_ref <- str_count(names_strains_sorted, name_ref)
}


for (sample in 1: nb_sample){
  name_strain <- names_strains_sorted[sample]
  fo_select <- subset(fo, sample == name_strain)
  
  for (kmer in seq(min_length, max_length)){
    cat(paste0("Loading data from ", names_strains_sorted[sample],".",kmer,".count.sequence.bed... "))
    
    seq_bed <- read.table(paste0(pathway_file_sequenceBedCount,names_strains_sorted[sample],".",kmer,".count.sequence.bed"), header = F)
    colnames(seq_bed) <- c("ID","Five_prime_pos","Three_prime_pos","Counts","Codons","Strand")
    
    seq_bed_adjusted <- adjust_to_offset(data = seq_bed,offsets = fo_select, site_of_interest = site)

    complete_counts <-
      codon_counts(data = seq_bed_adjusted,
        min_nbr = minimal_read,
        codons_table = df_triplets,
        filter_first = TRUE,
        in_frame = framed,
        elong = elongation)
    
    assign(paste0(names_strains_sorted[sample],".",kmer), complete_counts)
    cat("Done.\n")
  }
}

# Create tables with counts of every sample for each read length
for (kmer in kmers){
    df_final <- NULL
    
    for (sample in 1 : nb_sample){
        
      # Make a table with counts for a specific length of every condition
      dfi <- get(paste0(names_strains_sorted[sample],".",kmer))
      df_final <- cbind(df_final,dfi$Counts)
      colnames(df_final)[dim(df_final)[2]] <- paste0(names_strains_sorted[sample],".",kmer)
      row.names(df_final) <- triplets
        
    }
    assign(paste0("df_final_",kmer), df_final)
}

############################
# GRAPHS CREATION
############################
dir.create(paste0(pathway_file_graphs,"graphs_by_length/"),showWarnings = F)

cat("Length :\n")
for (kmer in kmers){
  cat(paste0("\t",kmer,"\n"))

  name_file_ref  <-  names_strains_sorted[1]
  df_final <- get(paste0("df_final_",kmer))

  df_final <- row_as_col(df_final, "Codon")

  df_final <- df_final[order(df_final$Codon),]
  df_final <- df_final %>% arrange(factor(Codon, levels = codon_start))
  bottom_rows <- df_final %>% filter(Codon %in% codons_stop)
  df_final <- df_final %>% filter(!(Codon %in% codons_stop))
  df_final <- bind_rows(df_final, bottom_rows)
  df_final$Codon <- factor(df_final$Codon, levels = unique(df_final$Codon))

  triplets <- row_as_col(df_triplets, "Codon")
  triplets <- triplets[, c(1, 3)]
  triplets <- triplets[order(match(triplets$Codon, df_final$Codon)),]
  df_final <- cbind(df_final, triplets[,2])
  df_final[df_final$AA_three_letters %like% "STOP",ncol(df_final)] <- "STOP"
  
  ###########################################
  # GRAPH NORMALISED
  ###########################################

  df_occupancy <- gather(df_final, key = "Sample", value = "Occupancy", -c(Codon,AA_three_letters))
  df_occupancy <- df_occupancy %>% separate(Sample, into = c("Condition", "Replicate"), sep = "\\.", remove = FALSE)
  # Keep the original order of the dataframe for the plot
  df_occupancy$Codon <- factor(df_occupancy$Codon, levels = unique(df_occupancy$Codon))
  df_occupancy$AA_three_letters <- factor(df_occupancy$AA_three_letters, levels = unique(df_occupancy$AA_three_letters))
  df_occupancy$Condition <- factor(df_occupancy$Condition, levels = unique(df_occupancy$Condition))

  histo <- df_occupancy %>%
    mutate(Codon = paste("<span style = 'color: ",
      ifelse(Codon %in% c(codon_start, codons_stop), "red", "black"),
      ";'>",
      Codon,
      "</span>", sep = "")) %>%
    ggplot(aes(x = Codon, y = Occupancy, fill = Condition, group = interaction(Sample, AA_three_letters, Condition))) + 
      geom_bar(stat = "identity", position=position_dodge(width=0.9), width = 0.5) +
      facet_grid2(
        .~AA_three_letters,
        space = 'free_x', scales = 'free_x', switch = 'x',
        strip = strip_themed(
        text_x = elem_list_text(colour = colors),
        by_layer_x = FALSE
        )) +
    labs(title = paste0(name_ref, " vs ", name_test, " - length=", kmer)) +
    scale_y_continuous(limits = c(0,max(df_occupancy$Occupancy)), expand = c(0, 0)) +
    scale_fill_manual(values = c("black", "red"))

  annotations <- list(
    createRectAnnotation(data = df_occupancy, col_x = "AA_three_letters", col_y = "Occupancy", cond = aa_start),
    createRectAnnotation(data = df_occupancy, col_x = "AA_three_letters", col_y = "Occupancy", cond = aa_stop, x_max = 3))

  #############################################
  # GRAPH MEAN WITH T-TEST SIGNIFICANCE RESULTS
  #############################################

  df_final_moyenne <- data.frame(
    Codon = df_final$Codon, 
    ref = df_final %>% select(starts_with(name_ref)) %>% rowMeans(), 
    test = df_final %>% select(starts_with(name_test)) %>% rowMeans(),
    AA_three_letters = df_final$AA_three_letters
  )

  df_mean_per_type <- pivot_longer(df_final_moyenne, cols = c(ref, test), names_to = "Type", values_to = "Mean")
  df_mean_per_type[df_mean_per_type$Type == "ref", 3] <- name_ref
  df_mean_per_type[df_mean_per_type$Type == "test", 3] <- name_test
  df_mean_per_type$AA_three_letters <- factor(df_mean_per_type$AA_three_letters, levels = unique(df_mean_per_type$AA_three_letters))
  df_mean_per_type$Type <- factor(df_mean_per_type$Type, levels = unique(df_mean_per_type$Type))

  df_pvalues_mean = data.frame(Codon = df_final$Codon, p.value = pvalues_calc(df_final, df_final$Codon))
  df_pvalues_mean$Significance <- sapply(df_pvalues_mean$p.value, map_pvalue_to_label)
  df_pvalues_mean <- join(df_pvalues_mean, df_mean_per_type, by = "Codon")
  
  df_pvalues_mean <- df_pvalues_mean %>% group_by(Codon) %>% filter(Mean == max(Mean))

  df_pvalues_mean <- df_pvalues_mean[order(match(df_pvalues_mean$AA_three_letters, df_mean_per_type$AA_three_letters)),]
  df_pvalues_mean$BracketPosition <- get_pos(df_pvalues_mean$AA_three_letters)
  df_pvalues_mean$Significance <- factor(df_pvalues_mean$Significance, levels = unique(df_pvalues_mean$Significance))

  plot_mean <- df_mean_per_type %>%
    mutate(
      Codon = paste("<span style = 'color: ",
        ifelse(Codon %in% c(codon_start, codons_stop), "red", "black"),
        ";'>",
        Codon,
        "</span>", sep = "")) %>%
    ggplot(aes(x = Codon, y = Mean, group = interaction(Type, AA_three_letters))) + 
      geom_bar(aes(fill = Type), stat="identity", position = position_dodge(width = 0.9), width = 0.5) +
      facet_grid2(.~AA_three_letters,
        space = 'free_x', scales = 'free_x', switch = 'x',
        strip = strip_themed(
        text_x = elem_list_text(colour = colors),
        by_layer_x = FALSE)) +
      labs(title = paste0(name_ref, " vs ", name_test, " - length=", kmer)) +
      scale_fill_manual(values = c("black", "red"))

  plot_mean_signif <- plot_mean + 
    scale_y_continuous(limits = c(0,max(df_mean_per_type$Mean + 0.001)), expand = c(0, 0)) +
    geom_signif(
      data = df_pvalues_mean,
      aes(comparisons = Codon,
        annotations = Significance, 
        xmin = BracketPosition - 0.5, 
        xmax = BracketPosition + 0.5, 
        y_position = Mean + 0.0002),
      manual = TRUE,
      size = 0.3,
      tip_length = 0,
      vjust = 0.5,
      textsize = 6,
      inherit.aes = FALSE)
  
  annotations <- list(
    createRectAnnotation(data = df_mean_per_type, col_x = "AA_three_letters", col_y = "Mean", cond = aa_start),
    createRectAnnotation(data = df_mean_per_type, col_x = "AA_three_letters", col_y = "Mean", cond = aa_stop, x_max = 3))
  
  ###########################################
  # GRAPH MEAN WITH STANDARD DEVIATION
  ###########################################
  
  # Standard deviation
  df_final_sd_wt <- apply(df_final[,2:(number_WT+1), drop = F], 1, sd, na.rm=F)
  df_final_sd_mut <- apply(df_final[,(number_WT+2):(number_WT+number_mut+1), drop = F], 1, sd, na.rm=F)

  erreur_wt <- df_final_sd_wt/ sqrt(number_WT)
  erreur_mut <- df_final_sd_mut /sqrt(nb_sample-number_mut)
  erreurs_fusion <- data.frame("Codon" = df_final$Codon, "ref" = erreur_wt, "test" = erreur_mut)
  erreurs_fusion <- pivot_longer(erreurs_fusion, cols = c(ref, test), names_to = "Type", values_to = "sd")
  erreurs_fusion <- cbind(erreurs_fusion[order(match(erreurs_fusion$Codon, df_mean_per_type$Codon)),], Mean = df_mean_per_type$Mean, AA_three_letters = df_mean_per_type$AA_three_letters)
  erreurs_fusion[erreurs_fusion$Type == "ref", 2] <- name_ref
  erreurs_fusion[erreurs_fusion$Type == "test", 2] <- name_test
  erreurs_fusion$Type <- factor(erreurs_fusion$Type, levels = unique(erreurs_fusion$Type))

  erreurs_fusion <- erreurs_fusion %>%
    mutate(
      Codon = paste("<span style = 'color: ",
        ifelse(Codon %in% c(codon_start, codons_stop), "red", "black"),
        ";'>",
        Codon,
        "</span>", sep = ""))

  plot_errors <- plot_mean + 
    scale_y_continuous(expand = c(0, 0)) +
    geom_errorbar(data = erreurs_fusion, 
      aes(ymin=Mean-sd, ymax=Mean+sd), 
      width=0.5,
      position=position_dodge(0.9))
  
  ###########################################
  # GRAPH LOGARITHM
  ###########################################

  df_final_log = data.frame(Codon = df_final_moyenne$Codon, log_ratio = log2(df_final_moyenne[, 3] / df_final_moyenne[, 2]))

  df_final_log_ordonne <- df_final_log[order(df_final_log[,2]),] # Sort by log values
  df_final_log_ordonne$colors <- rep(c("red","black"), nrow(df_final_log_ordonne)/2) 
  df_final_log_ordonne$colors[abs(df_final_log_ordonne$log_ratio) == Inf] <- "darkgreen"

  df_final_log_ordonne <- df_final_log_ordonne[c(
      which( df_final_log_ordonne$log_ratio < 0 & df_final_log_ordonne$log_ratio != -Inf
      ),
      which(df_final_log_ordonne$log_ratio == -Inf),
      which(df_final_log_ordonne$log_ratio == Inf),
      which( df_final_log_ordonne$log_ratio > 0 & df_final_log_ordonne$log_ratio != Inf
      )
  ), ]

  # Change infinite values for an arbitraty one (1/5 of the absolute min or max of the date) to be able to plot them
  df_final_log_ordonne$log_ratio[df_final_log_ordonne$log_ratio == -Inf] <-
    -max(abs(df_final_log_ordonne$log_ratio[abs(df_final_log_ordonne$log_ratio) != Inf]) /5)
  df_final_log_ordonne$log_ratio[df_final_log_ordonne$log_ratio == Inf] <-
    max(abs(df_final_log_ordonne$log_ratio[df_final_log_ordonne$log_ratio != Inf]) /5)

  df_final_log_ordonne$Codon <- factor(df_final_log_ordonne$Codon, levels =df_final_log_ordonne$Codon)

  #colors = ifelse(df_final_log_ordonne$Codon %in% c(codon_start, codons_stop), "red", "black")
  plot_log <- df_final_log_ordonne %>%
    ggplot(aes(x = Codon, y = log_ratio)) +
      geom_bar(aes(fill = Codon), stat="identity", position = position_dodge(width = 0.9), width = 0.5, show.legend = FALSE) +
      scale_fill_manual(values = df_final_log_ordonne$colors) +
      labs(title = paste0(name_ref, " vs ", name_test, " - length=", kmer)) +
      geom_hline(aes(yintercept=0)) +
      geom_text(aes(label = Codon, 
        y = ifelse(log_ratio > 0, -(max(abs(log_ratio)))/13, max(abs(log_ratio))/13)),
        colour = ifelse(df_final_log_ordonne$Codon %in% c(codon_start, codons_stop), "red", "black"),
        size = 5, angle=90) +
      scale_y_continuous(limits = c(min(df_final_log_ordonne$log_ratio), max(df_final_log_ordonne$log_ratio)+0.02), expand = c(0, 0))
  
  annotations <- list(
    createRectAnnotation(data = df_final_log_ordonne, col_x = "Codon", col_y = "log_ratio", cond = codon_start, xmin_adjust = 0.5, facet=TRUE),
    createRectAnnotation(data = df_final_log_ordonne, col_x = "Codon", col_y = "log_ratio", cond = codons_stop[1], xmin_adjust = 0.5, facet=TRUE),
    createRectAnnotation(data = df_final_log_ordonne, col_x = "Codon", col_y = "log_ratio", cond = codons_stop[2], xmin_adjust = 0.5, facet=TRUE),
    createRectAnnotation(data = df_final_log_ordonne, col_x = "Codon", col_y = "log_ratio", cond = codons_stop[3], xmin_adjust = 0.5, facet=TRUE))
  
}

###########################################
# GRAPH MEAN (not length specific)
###########################################

df_final_all <- NULL
for (sample in 1 : nb_sample){
  for (kmer in kmers){
    dfi <- get(paste0(names_strains_sorted[sample],".",kmer))
    df_final_all <- cbind(df_final_all,dfi$Counts)
    colnames(df_final_all)[dim(df_final_all)[2]] <- paste0(names_strains_sorted[sample],".",kmer)
  }
}
row.names(df_final_all) <- dfi$Codons

df_mean_by_sample <- means_in_df(data = df_final_all, mean_names = names_strains_sorted)

percentages_by_sample <- data.frame(apply(df_mean_by_sample, 2, FUN = function(x) (x*100)/sum(x)),
                                          row.names = row.names(df_mean_by_sample))
write.csv(percentages_by_sample,
          file = paste0(pathway_file_graphs,name_ref,"_vs_",name_test,"_percentage_by_sample_all_reads.csv"))

###############################################
# COUNTS WITH SIGNIFICANCE LEVELS (all lengths)
###############################################
df_mean_all <- means_in_df(df_final_all, c(name_ref,name_test))

percentages_by_condition <- data.frame(apply(df_mean_all, 2, FUN = function(x) (x*100)/sum(x)), row.names = row.names(df_mean_all))

write.csv(percentages_by_condition,
  file = paste0(pathway_file_graphs,name_ref,"_vs_",name_test,"_percentage_by_condition_all_reads.csv"))

write.csv(df_mean_all,
  file = paste0(pathway_file_graphs,name_ref,"_vs_",name_test,"_normalized_counts_by_condition_by_codon.csv"),
  quote = F)

data_mean_all <- gather(row_as_col(df_mean_all, "Codon"), key = "Condition", value = "Mean", -Codon)

df_final <- data.frame(df_final)
data_mean_all <- merge(data_mean_all, df_final[, c(1, ncol(df_final))], by = "Codon")
data_mean_all <- data_mean_all[order(match(data_mean_all$AA_three_letters, df_final$AA_three_letters)), ]
data_mean_all$AA_three_letters <- factor(data_mean_all$AA_three_letters, levels = unique(data_mean_all$AA_three_letters))
data_mean_all$Condition <- factor(data_mean_all$Condition, levels = unique(data_mean_all$Condition))

data_mean_by_sample <- row_as_col(df_mean_by_sample, "Codon")
data_mean_by_sample <- merge(data_mean_by_sample, df_final[, c(1, ncol(df_final))], by = "Codon")
data_mean_by_sample <- data_mean_by_sample[order(match(data_mean_by_sample$AA_three_letters, df_final$AA_three_letters)), ]

df_pvalues_by_sample = data.frame(Codon = data_mean_by_sample$Codon, p.value = pvalues_calc(data_mean_by_sample, data_mean_by_sample$Codon), AA_three_letters = data_mean_by_sample$AA_three_letters)
df_pvalues_by_sample$Significance <- sapply(df_pvalues_by_sample$p.value, map_pvalue_to_label)
df_pvalues_by_sample <- join(df_pvalues_by_sample, data_mean_all[,c(1,3)], by = "Codon")
df_pvalues_by_sample <- df_pvalues_by_sample %>%
  group_by(Codon) %>%
  filter(Mean == max(Mean))

df_pvalues_by_sample$BracketPosition <- get_pos(data_mean_by_sample$AA_three_letters)
df_pvalues_by_sample$AA_three_letters <- factor(df_pvalues_by_sample$AA_three_letters, levels = unique(df_pvalues_by_sample$AA_three_letters))

plot_mean_all <- data_mean_all %>%
  mutate(Codon = paste("<span style = 'color: ",
    ifelse(Codon %in% c(codon_start, codons_stop), "red", "black"),
    ";'>",
    Codon,
    "</span>", sep = "")) %>%
  ggplot(aes(x = Codon, y = Mean, group= interaction(Condition, AA_three_letters))) +
  geom_bar(aes(fill = Condition), stat="identity", position = position_dodge(width = 0.9), width = 0.5) +
  labs(title = paste0(name_ref," vs ",name_test, " - all reads")) +
  facet_grid2(.~AA_three_letters,
    space = 'free_x', scales = 'free_x', switch = 'x',
    strip = strip_themed(
    text_x = elem_list_text(colour = colors),
    by_layer_x = FALSE)) +
  scale_y_continuous(limits = c(0,max(data_mean_all$Mean + 0.005)), expand = c(0, 0)) +
  scale_fill_manual(values = c("black", "red")) +
  geom_signif(data = df_pvalues_by_sample,
    aes(comparisons = Codon,
      annotations = Significance, 
      xmin = BracketPosition - 0.5, 
      xmax = BracketPosition + 0.5, 
      y_position = Mean + 0.0002),
    size = 0.3,
    tip_length = 0,
    vjust = 0.5,
    textsize = 6,
    inherit.aes = FALSE,
    manual=TRUE)
    
annotations <- list(
  createRectAnnotation(data = data_mean_all, col_x = "AA_three_letters", col_y = "Mean", cond = aa_start),
  createRectAnnotation(data = data_mean_all, col_x = "AA_three_letters", col_y = "Mean", cond = aa_stop, x_max = 3))

tiff(file = paste0(pathway_file_graphs,name_ref,"_vs_",name_test,"_mean_by_codon_all_reads.tiff"),
     width = 1000,
     height = 500
)

  suppressWarnings(print(plot_mean_all + 
    customPlot + 
    annotations))

dev.off()

###########################################
# Amino-acid occupancy
# GRAPH MEAN (all lengths)
###########################################
df_mean_by_cond_aa <- data.frame(aggregate(x = cbind(df_mean_by_sample, df_triplets)[,1:dim(df_mean_by_sample)[2]],
                                           by = df_triplets["AA_three_letters"],
                                           FUN = sum),
                                 row.names = 1,
                                 check.names = FALSE
)

df_mean_all_aa <- data.frame(aggregate(x = cbind(df_mean_all, df_triplets)[,1:dim(df_mean_all)[2]],
                                           by = df_triplets["AA_three_letters"],
                                           FUN = sum),
                             row.names = 1,
                             check.names = FALSE
)

percentages_by_aa <- data.frame(apply(df_mean_all_aa, 2, FUN = function(x) (x*100)/sum(x)),
                                       row.names = row.names(df_mean_all_aa))
write.csv(percentages_by_aa,
          file = paste0(pathway_file_graphs,name_ref,"_vs_",name_test,"_percentage_by_amino-acid_all_reads.csv"))

write.csv(df_mean_all_aa,
    file = paste0(pathway_file_graphs,name_ref,"_vs_",name_test,"_normalized_counts_by_condition_by_aa.csv"),
    quote = F)

df_mean_all_aa <- row_as_col(df_mean_all_aa, "AA_three_letters")
df_mean_all_aa <- df_mean_all_aa[order(match(df_mean_all_aa$AA_three_letters, df_final$AA_three_letters)), ]
df_mean_all_aa <- gather(df_mean_all_aa, key = "Condition", value = "Mean", -AA_three_letters)
df_mean_all_aa$AA_three_letters <- factor(df_mean_all_aa$AA_three_letters, levels = unique(df_mean_all_aa$AA_three_letters))
df_mean_all_aa$Condition <- factor(df_mean_all_aa$Condition, levels = unique(df_mean_all_aa$Condition))

df_mean_by_cond_aa <- row_as_col(df_mean_by_cond_aa, "AA_three_letters")
df_mean_by_cond_aa <- df_mean_by_cond_aa[order(match(df_mean_by_cond_aa$AA_three_letters, df_mean_all_aa$AA_three_letters )), ]

df_pvalues_by_aa = data.frame("AA_three_letters" = df_mean_by_cond_aa$AA_three_letters, p.value = pvalues_calc(df_mean_by_cond_aa, df_mean_by_cond_aa$AA_three_letters))
df_pvalues_by_aa$Significance <- sapply(df_pvalues_by_aa$p.value, map_pvalue_to_label)
df_pvalues_by_aa <- join(df_pvalues_by_aa, df_mean_all_aa[,c(1,3)], by = "AA_three_letters")
df_pvalues_by_aa <- df_pvalues_by_aa %>%
    group_by(AA_three_letters) %>%
    filter(Mean == max(Mean))

df_pvalues_by_aa$BracketPosition <- seq(1, length(df_pvalues_by_aa$AA_three_letters), 1)

annotations <- list(
  createRectAnnotation(data = df_mean_all_aa, col_x = "AA_three_letters", col_y = "Mean", cond = aa_start, xmin_adjust = 0.5, facet=TRUE),
  createRectAnnotation(data = df_mean_all_aa, col_x = "AA_three_letters", col_y = "Mean", cond = aa_stop, xmin_adjust = 0.5, facet=TRUE))

df_mean_all_aa <- df_mean_all_aa %>% 
    mutate(AA_three_letters = paste("<span style = 'color: ",
      ifelse(AA_three_letters %in% c(aa_start, paste0(aa_stop, "-", codons_stop)), "red", "black"),
      ";'>",
      AA_three_letters,
      "</span>", sep = ""))

df_mean_all_aa$AA_three_letters <- factor(df_mean_all_aa$AA_three_letters, levels = unique(df_mean_all_aa$AA_three_letters))

plot_mean_aa <- 
  df_mean_all_aa %>%
  ggplot(aes(x = AA_three_letters, y = Mean)) +
    geom_bar(aes(fill = Condition), stat="identity", position = position_dodge(width = 0.9), width = 0.7) +
    scale_fill_manual(values = c("black", "red")) + 
    labs(title = paste0(name_ref," vs ",name_test, " - all reads")) +
    scale_y_continuous(limits = c(0,max(df_mean_all_aa$Mean + 0.005)), expand = c(0, 0)) +
    geom_signif(annotations = df_pvalues_by_aa$Significance, 
      xmin = df_pvalues_by_aa$BracketPosition - 0.5, 
      xmax = df_pvalues_by_aa$BracketPosition + 0.5, 
      y_position = df_pvalues_by_aa$Mean + 0.0005,
      size = 0.3,
      tip_length = 0,
      vjust = 0.5,
      textsize = 6)

tiff(file = paste0(pathway_file_graphs,name_ref,"_vs_",name_test,"_mean_by_amino-acid_all_reads.tiff"),
     width = 1000,
     height = 500
)
  
  suppressWarnings(print(plot_mean_aa +
    customPlot + 
    annotations))

dev.off()