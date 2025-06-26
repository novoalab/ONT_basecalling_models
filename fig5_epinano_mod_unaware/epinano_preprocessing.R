## EPINANO PREPROCESSING: 
# give the epinano flow output directory and will return the columns that contains 
# the necessary information. 

# EXAMPLE OF CALLING: 
# Rscript epinano_preprocessing.R \
# --working_directory "/no_backup/enovoa/nextflow_outputs/curlcakes_RNA241003_Pool1_francesco/mopmod_outfolder_problematic_bc_sup/epinano_flow/" \
# --bc_names_remove "pool1_subsampled_pod5---bc_7" "pool1_subsampled_pod5---bc_8" "pool1_subsampled_pod5---bc_9" \
# --bc_mod_df_name "/nfs/users/enovoa/fpelizzari/CC_fast_hac_sup_comparison/bc_modifications_pool1.tsv" \
# --output_dir "./epinano_processed/sup" 

# ALSO IMPORTANT! the output is assumed to be constant! 
# <barcode>_<modification>.tsv! 

# library loading 
library(argparse)
library(data.table)
suppressPackageStartupMessages(suppressWarnings(library(tidyverse)))

parser <- ArgumentParser()
parser$add_argument("--working_directory", required=TRUE,
                    help="Path to the working directory")
parser$add_argument("--bc_names_remove", nargs="+", 
                    help="List of barcodes to remove", 
                    default = NULL)
parser$add_argument("--bc_mod_df_name", required=TRUE,
                    help="Path to the barcode modification file (TSV)")
parser$add_argument("--output_dir", required=TRUE,
                    help="Output directory path")
args <- parser$parse_args()



working_directory <- args$working_directory
bc_names_remove <- args$bc_names_remove
bc_mod_df_name <- args$bc_mod_df_name
output_dir <- args$output_dir


 

file_list <- list.files(path = working_directory)
file_list_f <- file_list[!grepl("\\.pdf$|minus", file_list)]
print(file_list_f)
# bc to keep and they modification
bc_mod_df <- read.table(file = bc_mod_df_name, header = T)

bc_to_keep <- paste0("bc_", bc_mod_df$BC, "\\.")
print(bc_to_keep)
file_list_f_bc <- file_list_f[grepl(paste0(bc_to_keep, collapse = "|"), 
                                    file_list_f)]

# also remove other files if needed (example: the subsampled pod5 of the problematic bcs)
if(!is.null(bc_names_remove)) {
  file_list_f_bc <- file_list_f_bc[!grepl(paste(bc_names_remove, collapse = "|"), 
                                          file_list_f_bc)]
}

sorted_file_list <- file_list_f_bc[order(as.numeric(gsub(".*bc_(\\d+).*", 
                                                         "\\1", file_list_f_bc)))]

csv_reader_get_summ_err <- function(working_directory, 
                                    csv_name, 
                                    output_dir, 
                                    modification, bc) {
  csv_to_read <- paste0(working_directory, csv_name)
  the_df <- fread(csv_to_read, sep = ",", header = T)
  the_df <- the_df %>% 
    mutate(sum_err = mis + ins + del) %>% 
    mutate(match_plus_sum_err = mat + sum_err) %>% 
    mutate(mat = mat/match_plus_sum_err) %>% 
    mutate(sum_err = sum_err/match_plus_sum_err) %>% 
    select("#Ref", "pos", "sum_err")
  colnames(the_df) <- c("ref", "pos", "stat")
  
  # create the out_dir if not already there
  dir.create(output_dir, showWarnings = F, recursive = T)
  file_name <- paste0("bc", bc, "_", modification,".tsv")
  file_name <- file.path(output_dir, file_name)
  write.table(x = the_df, file = file_name,
              quote = F, sep = "\t", row.names = F)
}

print(sorted_file_list)

for (index in 1:length(sorted_file_list)){
  csv_reader_get_summ_err(working_directory = working_directory, 
                          csv_name = sorted_file_list[index], 
                          output_dir = output_dir, 
                          modification = bc_mod_df[index, "modification"], 
                          bc = bc_mod_df[index, "BC"])
}



