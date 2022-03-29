library(stringr)

local_path <- "/data/"

# For file size (entire fastq usually too big for R)
ibs <- 8192

args <- commandArgs(trailingOnly=TRUE)
fastq <- args[1]
sample <- sub(".fastq.gz", "", basename(fastq))

# Unzip and select only nucleotide sequences in fastq file
try(File <- system(paste0("gzip -cd ", fastq," 2> /dev/null | dd ibs=",ibs," count=10 | grep -P '^[ATGC]+$'"), intern = TRUE))
File <- as.data.frame(File)

# Default sequence size checked
Adapt_size <- 10

List_adapters <- NULL

# Number of lines to check (arbitraty)
lines_check=500
if(nrow(File) < lines_check) {lines_check <- nrow(File)}

for(Ligne in 1 : lines_check)
{
  #Scan line
  for(Position1 in 1 : ((nchar(File[Ligne,1]) - Adapt_size )))
  {
    # Sequence checked in whole file (found at line Ligne from position Position1 with Adapt_size as length)
    Pattern <- str_sub(File[Ligne,1], Position1, (Position1 + Adapt_size))
    # Number of this sequence occurrences in all fastq lines selected
    Occurence <- suppressWarnings(str_count(File, Pattern))
    
    # If this sequence is found in less than 70% (0.7 : arbitrary) of the fastq lines : the window is moved
    if(Occurence < (nrow(File) * 0.7)) {next}
    # Enlarging the window of the checked sequence
    for(Position2 in (Position1 + Adapt_size) : nchar(File[Ligne,1]))
    {
      # New checked sequence
      Pattern2 <- str_sub(File[Ligne,1], Position1, Position2)
      New_Occurence <- suppressWarnings(str_count(File, Pattern2))
      if(New_Occurence < (nrow(File) * 0.7))
      {
        Pattern2 <- str_sub(File[Ligne,1], Position1, Position2-1)
        break
      }
      # If the checked sequence is big enough (more than 20nt (arbitrary))
      if(nchar(Pattern2) == 20) {break}
    }
    
    # Sequence is saved as a potential adapter
    List_adapters <- rbind(List_adapters, Pattern2)
    break
  }
  df <- as.data.frame(table(List_adapters))
  
  # If more than half the checked lines (arbitrary) have been used to make patterns
  # and if the best checked sequence appears in more than 70% (0.7 : arbitrary) of the fastq lines
  # The checked sequence is kept as an adapter
  if(Ligne > (lines_check * 0.5) & (max(df$Freq)/Ligne) > 0.7)
  {
    print(Pattern2)
    write(Pattern2, paste0(local_path, "RESULTS/adapter_lists/", sample, ".txt"), append = TRUE)
    break
  }
}
