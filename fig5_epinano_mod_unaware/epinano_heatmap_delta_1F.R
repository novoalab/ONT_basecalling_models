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


dfs_split_by_modification <- split(all_full_stats_all_models, 
                                   all_full_stats_all_models$modification)

# first thing: do this for the unmodified and then remove it from the list 
unmodified_df <- dfs_split_by_modification$UNM
dfs_split_by_modification$UNM <- NULL

# merge it position by position --> merge the three different replicates
unm_median_df <- unmodified_df %>% 
  group_by(central_coord, reference, model) %>% 
  summarise(minus_two_median = median(minus_two_stat), 
            minus_one_median = median(minus_one_stat),
            central_median = median(central_stat), 
            plus_one_median = median(plus_one_stat), 
            plus_two_median = median(plus_two_stat))

# consider one modification by one, then merge all of them in the same df
# do the facet boxplot
all_deltas_all_modifications <- NULL
for (modification in names(dfs_split_by_modification)) {
  # get the df corresponding to that modification 
  one_mod_df <- dfs_split_by_modification[[modification]]
  # get the median of the three replicates
  mod_median_df <- one_mod_df %>% 
    group_by(central_coord, reference, model) %>% 
    summarise(minus_two_median = median(minus_two_stat), 
              minus_one_median = median(minus_one_stat),
              central_median = median(central_stat), 
              plus_one_median = median(plus_one_stat), 
              plus_two_median = median(plus_two_stat)) 
  
  # get the delta statistic for that specific modification AND the UNM 
  # specifically do MOD - UNM!
  merged_df <- inner_join(unm_median_df, mod_median_df, 
                          by = c("reference", "central_coord", "model")) %>%
    mutate(delta_minus_two = minus_two_median.y - minus_two_median.x, 
           delta_minus_one = minus_one_median.y - minus_one_median.x, 
           delta_central = central_median.y - central_median.x, 
           delta_plus_one = plus_one_median.y - plus_one_median.x, 
           delta_plus_two = plus_two_median.y - plus_two_median.x, )%>%
    select(reference, central_coord, model, 
           delta_minus_two, delta_minus_one, delta_central, 
           delta_plus_one, delta_plus_two) 
  merged_df$modification <- rep(modification, nrow(merged_df))
  
  # let's see if we can do the plotting 
  all_deltas_all_modifications <- rbind(all_deltas_all_modifications, 
                                        merged_df)
}

# define the colors - take from gregor 
colours <- list('model' = c("fast" = "#F8766D", 
                            "hac" =  "#00BA38", 
                            "sup" =  "#619CFF"),
                'modification' = c('m6A' = '#381a61', 'm5C' = '#f9d14a', 
                                   '5hmC' = "#ab3329", 'ac4C' = "#88a0dc",
                                   'Y' = '#e78429', 'm1Y' = '#ed968c',
                                   'm5U' = '#7c4b73'))

#############
# again from gregor 
# https://github.com/novoalab/basecalling_models/blob/main/scripts/R/in_vitro/heatmap_v2.R

## BUILD THE CORRECT DF 
# temp_df <- all_full_stats_all_models %>% 
temp_df <- all_deltas_all_modifications %>% 
  group_by(model, modification) %>% 
  # for each kmer position get the median 
  summarize(minus_two_median = median(delta_minus_two, na.rm = T), 
            minus_one_median = median(delta_minus_one, na.rm = T), 
            central_median = median(delta_central, na.rm = T), 
            plus_one_median = median(delta_plus_one, na.rm = T), 
            plus_two_median = median(delta_plus_two, na.rm = T)) %>% 
  filter(modification != "UNM") %>% 
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
          title = expression(Delta ~ "SumErr")), 
        column_title = expression(Delta ~ "SumErr"))

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
        column_title = expression(Delta ~ "SumErr - Z norm. by modification"))

filename_normed <- paste0(output_dir, "/", "normalized", file_name)
tidyHeatmap::save_pdf(.heatmap = heatmap_normed,
                      filename = filename_normed, 
                      width = 8,
                      height = 6,
                      units = "in")
