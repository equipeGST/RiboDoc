# Purpose : Add a "/" at the end of a string if missing
# Input : path to a folder (chars)
# Output : path to the folder with a '/' at the end
# Launch : pathway_handling("/path/to/folder")

# Libraries
library(stringr)

# Check the last character of pathway
last_character_determine <- function(pathway_file){
  last_character <- str_sub(pathway_file, nchar(pathway_file), nchar(pathway_file))
  return(last_character)
}

# Adds a "/" at the end of the folder pathways if missing
add_slash_to_path <- function(pathway_file, last_character) {
  if (last_character != "/") {
    return(paste0(pathway_file,"/"))
  } else {
    return(pathway_file)
  }
}

# Handle pathway string to add a "/" if missing
pathway_handling <- function(pathway_file) {
  last_character <- last_character_determine(pathway_file)
  pathway_file_final <- add_slash_to_path(pathway_file,last_character)
  return(pathway_file_final)
}