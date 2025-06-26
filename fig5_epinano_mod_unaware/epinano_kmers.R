# THIS SCRIPT IS USED TO GET THE KMER INFO given the processed output file. 
# add some explanation 
suppressPackageStartupMessages({
  suppressWarnings(library(Biostrings))
  suppressWarnings(library(argparse))
  suppressWarnings(library(tidyverse))
})

parser <- ArgumentParser()

parser$add_argument("--working_dir", type = "character",
                    help = "Directory containing input tables", 
                    required = T)
parser$add_argument("--reference_file", type = "character",
                    help = "Path to the reference fasta file", 
                    default = "/users/enovoa/fpelizzari/my_references/curlcake_constructs_EcoRV_BamHI_digestion.fasta")
parser$add_argument("--output_dir", type = "character",
                    help = "Output directory", 
                    required = T)

args <- parser$parse_args()

working_dir <- args$working_dir
reference_file <- args$reference_file
output_dir <- args$output_dir 

# quick check on the final slash oh working dir and output dir 
if (!grepl("/$", working_dir)) {
  working_dir <- paste0(working_dir, "/")
}
if (!grepl("/$", output_dir)) {
  output_dir <- paste0(output_dir, "/")
}

references <- readDNAStringSet(reference_file)

# get the different models - in working dir expected one folder per model 
models <- list.files(working_dir)

# FUNCTUON DEFININTIONS 
# check if the whole kmer is present in the stat table - 
# it has statistics for the whole base otherwise return null 
check_all_coordinates <- function(coordinates, positions) {
  if (all(coordinates %in% positions)) {
    return(T)
  } 
  else {return(F)}
}
# given a position, eg. -1, return the stat in the kmer for that position
get_stat_for_position <- function(coordinate, table, the_position) {
  col_value <- coordinate[the_position]
  stat <- table %>%
    filter(pos == col_value) %>%
    pull(stat)
  return(stat)
}

# given the coordinate of a kmer, return its median
get_kmer_median <- function(coordinate, table) {
  median <- median(stat_table_red[stat_table_red$pos %in% coordinate, "stat"])
  return(median)
}

# get the different sequences of the kmers 
get_full_kmer_sequence <- function(coordinate, dna) {
  return(as.character(substr(dna, min(coordinate), max(coordinate))))
}
# and get the central base 
get_base_for_position <- function(coordinate, dna, position) {
  index <- coordinate[position]
  return(as.character(substr(dna, index, index)))
}

# get the coordinate from the references 
get_base_indices <- function(dna_set, the_base) {
  lapply(dna_set, function(dna) {
    if (the_base == "UNM") {
      return(1:nchar(as.character(dna)))
    } else {
      dna_chars <- strsplit(as.character(dna), "")[[1]]
      return(which(dna_chars == the_base) - 1)
    }
  })
}

# iterate over the different models 
for (the_model in models){
  working_dir_mod <- paste0(working_dir, the_model)
  print(working_dir_mod)
  files_list <- list.files(working_dir_mod)
  output_dir_mod <- paste0(output_dir, the_model)
  # get the modification 
  modification <- sub(".*_(.*)\\..*", "\\1", files_list)
  # get the modified base
  modified_base <- ifelse(modification == "UNM", "UNM", 
                          substr(modification, nchar(modification), 
                                 nchar(modification)))
  
  # if it's Y (pseudourudine) then put T
  modified_base <- ifelse(modified_base != "Y", modified_base, "T")
  # if it's U (urudine) then put T
  modified_base <- ifelse(modified_base != "U", modified_base, "T")
  
  df_file_mod <- data.frame(files_list, modification, 
                            modified_base)
  print(df_file_mod)
  # get the full file names, then build a df (keep everything together)
  df_file_mod$file_full_names <- paste0(working_dir, the_model, 
                                        "/", files_list)
  
  print(df_file_mod)
  
  # get the barcodes 
  df_file_mod$bc <- sub("_.*", "", df_file_mod$files_list)
  
  # create the output dir 
  dir.create(output_dir_mod, showWarnings = F, recursive = T)
  print("---------")

  for (row_index in 1:nrow(df_file_mod)) {
    print(row_index)
    stat_table <- read.table(df_file_mod[row_index, "file_full_names"], 
                             header = T)
    
    barcode <- df_file_mod[row_index, "bc"]
    mod_base <- df_file_mod[row_index, "modified_base"]
    modification <- df_file_mod[row_index, "modification"]
    
    # THIS RETURN THE OUTPUT IN AN "EPINANO - PYTHON FRIENDLY WAY"
    # everything else --> add +1 
    # ex. when getting the sequence and the central bases
    mod_indeces <- get_base_indices(references, mod_base)
    
    df_all_references <- NULL
    
    for (ref_name in names(mod_indeces)) {
      print(ref_name)
      ref <- references[[ref_name]]
      stat_table_red <- stat_table[which(stat_table$ref == ref_name), ]
      
      coords <- mod_indeces[[ref_name]]
      
      # get all the possible kmers - the coordinates 
      kmer_coords <- lapply(coords, function(coord){
        return((coord - 2):(coord + 2))
      })
      
      # add a small check - kmer_coords must be > 0 !! otherwise invalid sequence
      kmer_coords <- kmer_coords[sapply(kmer_coords, function(el) {
        return(all(el > 0))
      })]
      
      # keep only the coord that are inside the stat table fully 
      kmer_coords <- kmer_coords[sapply(kmer_coords, 
                                        check_all_coordinates, 
                                        stat_table_red$pos)]
      
      # get the sequence, the central base, the median, the central stat 
      all_central_coordinates <- sapply(kmer_coords, `[[`, 3) # after filtering
      
      # when you get things from the biostrings - add +1!!
      kmer_coords_plus_one <- lapply(kmer_coords, function(x) { x + 1})
      
      all_sequences <- sapply(kmer_coords_plus_one, 
                              get_full_kmer_sequence, ref)
      all_central_bases <- sapply(kmer_coords_plus_one, 
                                  get_base_for_position, ref, position = 3)
      
      
      
      all_medians <- sapply(kmer_coords, get_kmer_median, stat_table_red)
      # get each single stat for each kmer  
      all_minus_two_stats <- sapply(kmer_coords, get_stat_for_position, stat_table_red, 
                                    the_position = 1)
      all_minus_one_stats <- sapply(kmer_coords, get_stat_for_position, stat_table_red, 
                                    the_position = 2)
      all_central_stats <- sapply(kmer_coords, get_stat_for_position, stat_table_red, 
                                  the_position = 3)
      all_plus_one_stats <- sapply(kmer_coords, get_stat_for_position, stat_table_red, 
                                   the_position = 4)
      all_plus_two_stats <- sapply(kmer_coords, get_stat_for_position, stat_table_red, 
                                   the_position = 5)
      
      
      df_one_reference <- data.frame(reference_name = rep(ref_name, length(all_central_coordinates)), 
                                     central_coordinate = all_central_coordinates,
                                     sequence = all_sequences, 
                                     central_base = all_central_bases,
                                     median = all_medians, 
                                     minus_two_stat = all_minus_two_stats, 
                                     minus_one_stat = all_minus_one_stats, 
                                     central_stat = all_central_stats, 
                                     plus_one_stat = all_plus_one_stats, 
                                     plus_two_stat = all_plus_two_stats)
      # bc at the end the number 
      df_all_references <- rbind(df_all_references, 
                                 df_one_reference)
    }
    file_name <- paste0(barcode, "_", modification, "_kmers.tsv")
    print(file_name)
    full_file_name <- file.path(output_dir_mod, file_name)
    # remove rows with nas if any 
    df_all_references <- df_all_references %>% drop_na()
    write.table(x = df_all_references, 
                file = full_file_name, 
                quote = F, sep = "\t", col.names = T, row.names = F)
  }
  
}



