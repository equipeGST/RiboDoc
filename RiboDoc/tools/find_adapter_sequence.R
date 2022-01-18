library(stringr)

# For file size (entire fastq usually too big for R)
ibs <- 8192

fastq <- list.files("/data/fastq/")[1]

# Unzip and select only nucleotide sequences in fastq file
try(File <- system(paste0("gzip -cd ", paste0("/data/fastq/",fastq)," | dd ibs=",ibs," count=10 | grep -P '^[ATGC]+$' "), intern = TRUE))
File <- as.data.frame(File)

# Default sequence size checked 
Adapt_size <- 10

List_adapters <- NULL

# Number of lines to check (arbitraty)
lines_check=500
if(nrow(File) > lines_check) {lines_check <- nrow(File)}

for(Ligne in 1 : lines_check)
{
  #Scan line
  for(Position1 in 1 : ((nchar(File[Ligne,1]) - Adapt_size )))
  {
    # Sequence checked in whole file (found at line Ligne from position Position1 with Adapt_size as length)
    Pattern <- str_sub(File[Ligne,1], Position1, (Position1 + Adapt_size))
    # Number of this sequence occurrences in all fastq lines selected
    Occurence <- str_count(File, Pattern)
    
    # If this sequence is found in less than 80% (0.8 : arbitrary) of the fastq lines : the window is moved
    if(Occurence < (nrow(File) * 0.8)) {next}
    # Enlarging the window of the checked sequence
    for(Position2 in (Position1 + Adapt_size) : nchar(File[Ligne,1]))
    {
      # New checked sequence
      Pattern2 <- str_sub(File[Ligne,1], Position1, Position2 )
      # If the checked sequence is big enough (more than 20nt (arbitrary))
      if(nchar(Pattern2) == 20) {break}
    }
    
    # Sequence is saved as a potential adapter
    List_adapters <- rbind(List_adapters, Pattern2)
    break
  }
  df <- as.data.frame(table(List_adapters))
  
  # If more than 100 lines (arbitrary) have been used to make patterns
  # and if the best checked sequence appears in more than 80% (0.8 : arbitrary) of the fastq lines
  # The checked sequence is kept as an adapter
  if(Ligne > 100 & (max(df$Freq)/Ligne) > 0.8 )
  {
    print(Pattern2)
    write(Pattern2, "/data/RESULTS/adapter_list.txt", append = TRUE)
    break
  }
}

