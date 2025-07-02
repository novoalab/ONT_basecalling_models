# HEATMAPS - 1F 
# example of calling --> 
# Rscript heatmaps_1F.R --working_dir epinano_kmers/ 
# --file_name "epinano_heatmap.pdf" --output_dir "plots/epinano"

suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(tidyHeatmap))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(viridisLite))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

parser$add_argument("--working_dir", type = "character",
                    help = "Directory containing input tables")
parser$add_argument("--file_name", type = "character",
                    help = "File name for the plot")
parser$add_argument("--output_dir", type = "character",
                    help = "Output directory")
args <- parser$parse_args()

working_dir <- args$working_dir
output_dir <- args$output_dir
file_name <- args$file_name



folder_list <- list.files(working_dir) 
dir.create(output_dir, showWarnings = F, recursive = T)

all_full_stats_all_models <- NULL
for (model in folder_list) {
  # enter the folder of the model 
  model_dir <- paste0(working_dir, "/", model)
  model_files <- list.files(model_dir)
  all_barcodes <- sapply(strsplit(model_files, "_"), function(e) {return(e[1])})
  all_modifications <- sapply(strsplit(model_files, "_"), function(e) {return(e[2])})
  
  df_file_mod <- data.frame(file_name = model_files, 
                            full_file_name = file.path(model_dir, model_files),
                            modification = all_modifications, 
                            barcode = all_barcodes)
  # read all the different stats
  all_minus_two_stats <- lapply(df_file_mod$full_file_name, function(file_name) {
    read.table(file_name, header = T)$minus_two_stat})
  all_minus_one_stats <- lapply(df_file_mod$full_file_name, function(file_name) {
    read.table(file_name, header = T)$minus_one_stat})
  all_central_stats <- lapply(df_file_mod$full_file_name, function(file_name) {
    read.table(file_name, header = T)$central_stat})
  all_plus_one_stats <- lapply(df_file_mod$full_file_name, function(file_name) {
    read.table(file_name, header = T)$plus_one_stat})
  all_plus_two_stats <- lapply(df_file_mod$full_file_name, function(file_name) {
    read.table(file_name, header = T)$plus_two_stat})
  
  all_central_coords <- lapply(df_file_mod$full_file_name, function(file_name) {
    read.table(file_name, header = T)$central_coord})
  all_references <- lapply(df_file_mod$full_file_name, function(file_name) {
    read.table(file_name, header = T)$reference_name})
  
  df_central_stat_that_model <- data.frame(
    minus_two_stat = unlist(all_minus_two_stats), 
    minus_one_stat = unlist(all_minus_one_stats), 
    central_stat = unlist(all_central_stats),
    plus_one_stat = unlist(all_plus_one_stats), 
    plus_two_stat = unlist(all_plus_two_stats), 
    central_coord = unlist(all_central_coords),
    reference = unlist(all_references), 
    barcode = rep(df_file_mod$barcode, times = sapply(all_central_coords, length)),
    modification = rep(df_file_mod$modification, times = sapply(all_central_coords, length)), 
    model = rep(model, sum(sapply(all_central_coords, length))))
  
  # update the df with all the models 
  all_full_stats_all_models <- rbind(all_full_stats_all_models, 
                                     df_central_stat_that_model)
}




# define the colors - take from gregor 
colours <- list('model' = c("fast" = "#F8766D", 
                            "hac" =  "#00BA38", 
                            "sup" =  "#619CFF"),
                'modification' = c('m6A' = '#381a61', 'm5C' = '#f9d14a', 
                                   '5hmC' = "#ab3329", 'ac4C' = "#88a0dc",
                                   'Y' = '#e78429', 'm1Y' = '#ed968c',
                                   'm5U' = '#7c4b73', "UNM" = 'grey50'))

#############
# again from gregor 
# https://github.com/novoalab/basecalling_models/blob/main/scripts/R/in_vitro/heatmap_v2.R

# temp_df <- all_full_stats_all_models %>% 
temp_df <- all_full_stats_all_models %>% 
  group_by(model, modification) %>% 
  # for each kmer position get the median 
  summarize(minus_two_median = median(minus_two_stat, na.rm = T), 
            minus_one_median = median(minus_one_stat, na.rm = T), 
            central_median = median(central_stat, na.rm = T), 
            plus_one_median = median(plus_one_stat, na.rm = T), 
            plus_two_median = median(plus_two_stat, na.rm = T)) %>% 
  # ungroup otherwise select doesn't allow removing this column
  ungroup() 
# dont know why but it was giving problems with the others 
colnames(temp_df) <- c("model", "modification", "-2", "-1", "0", "+1", "+2")
temp_df <- temp_df %>% 
  # Finalize order for plotting
  arrange(modification, model) %>% 
  # Create new column with combined information and remove others
  mutate(mod_model = paste0(modification, '_', model), .before = model) 

# BUILD THE ANNOTATIONS
ann <- data.frame(temp_df$model)
ann$modification <- temp_df$modification
colnames(ann) <- c('model','modification')
leftAnn <- HeatmapAnnotation(df = ann,
                             which = 'row',
                             col = colours, 
                             annotation_width = unit(c(1, 4), 'cm'),
                             gap = unit(1, 'mm'),
                             show_annotation_name = FALSE,
                             annotation_legend_param = list(
                               labels_gp = grid::gpar(fontsize = 15),
                               title_gp = grid::gpar(fontsize = 17)))

###### HEATMAP NOT NORMALIZED #######
# SELECT FOR THE MATRIX ANNOTATION
temp_mx <- temp_df %>%  
  select(-model,-modification) %>%  
  # convert column to rownames and finally convert to matrix needed for corrplot
  column_to_rownames(var = "mod_model") %>% 
  data.matrix()

# Reorder rows: UNM_* first, then the rest
temp_mx <- temp_mx[c(grep("^UNM", rownames(temp_mx)),
                     grep("^(?!UNM)", rownames(temp_mx), perl = TRUE)), ]
ann$modification <- factor(ann$modification,
                           levels = c("UNM", setdiff(unique(ann$modification), "UNM")))

# BUILD THE HEATMAP 
heatmap_not_normed <- Heatmap(temp_mx, 
                              cluster_rows = F, cluster_columns = F,
                              column_names_rot = 0,
                              col = viridis(100),
                              column_names_gp = grid::gpar(fontsize = 20, fontface = "bold"),
                              show_row_names = F,
                              left_annotation = leftAnn,
                              row_split = ann$modification,
                              heatmap_legend_param = list(
                                title = expression("SumErr")), 
                              column_title = expression("SumErr"))



# AND SAVE
filename_not_normed <- paste0(output_dir, "/", "not_normalized", file_name)
tidyHeatmap::save_pdf(.heatmap = heatmap_not_normed,
                      filename = filename_not_normed, 
                      width = 8,
                      height = 6,
                      units = "in")

##### NORMALIZED MATRIX ######

### This part is ugly - ok gregor :'( but not gonna change it 
temp_mx_s <- temp_df %>%  
  # pivot_longer to be able to apply scaling on modification groups
  pivot_longer(cols = `-2`:`+2`,
               names_to = "position",
               values_to = "sumerr") %>% 
  #group and scale
  group_by(modification) %>% 
  #vector conversion required otherwise df returned
  mutate(sumerr_scaled = as.vector(scale(sumerr))) %>% 
  ungroup() %>% 
  select(-sumerr) %>% 
  #repivot to wide
  pivot_wider(names_from = "position",
              values_from = "sumerr_scaled") %>% 
  select(-model,-modification) %>%  
  # convert column to rownames and finally convert to matrix needed for corrplot
  column_to_rownames(var = "mod_model") %>% 
  data.matrix()

# Reorder rows: UNM_* first, then the rest
temp_mx_s <- temp_mx_s[c(grep("^UNM", rownames(temp_mx_s)),
                     grep("^(?!UNM)", rownames(temp_mx_s), perl = TRUE)), ]

ann$modification <- factor(ann$modification,
                           levels = c("UNM", setdiff(unique(ann$modification), "UNM")))


heatmap_normed <- Heatmap(temp_mx_s, 
                          cluster_rows = F, cluster_columns = F,
                          column_names_rot = 0,
                          col = viridis(100),
                          column_names_gp = grid::gpar(fontsize = 20, fontface = "bold"),
                          show_row_names = F,
                          left_annotation = leftAnn,
                          row_split = ann$modification,
                          heatmap_legend_param = list(
                            title = "Z-score"), 
                          column_title = expression("SumErr - Z norm. by modification"))

filename_normed <- paste0(output_dir, "/", "normalized", file_name)
tidyHeatmap::save_pdf(.heatmap = heatmap_normed,
                      filename = filename_normed, 
                      width = 8,
                      height = 6,
                      units = "in")
