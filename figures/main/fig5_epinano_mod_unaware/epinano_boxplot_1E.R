### BOXPLOT 1E 
# one boxplot per mod - in all the boxplot we have the three models but compare
# the central statistics with the UNM


## example of calling: 
# Rscript boxplot_1E.R --working_dir epinano_kmers --output_dir "./plots/epinano" --file_name epinano_1E.pdf
suppressPackageStartupMessages(library(EnvStats))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggpubr))

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

folder_list <- list.files(working_dir) # each folder should correspond to the model

all_central_stats_all_models <- NULL
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
  
  all_central_stats <- lapply(df_file_mod$full_file_name, function(file_name) {
    read.table(file_name, header = T)$central_stat})
  all_central_coords <- lapply(df_file_mod$full_file_name, function(file_name) {
    read.table(file_name, header = T)$central_coord})
  all_references <- lapply(df_file_mod$full_file_name, function(file_name) {
    read.table(file_name, header = T)$reference_name})
  
  df_central_stat_that_model <- data.frame(
    central_stat = unlist(all_central_stats),
    central_coord = unlist(all_central_coords),
    reference = unlist(all_references), 
    barcode = rep(df_file_mod$barcode, times = sapply(all_central_coords, length)),
    modification = rep(df_file_mod$modification, times = sapply(all_central_coords, length)), 
    model = rep(model, sum(sapply(all_central_coords, length))))
  
  # update the df with all the models 
  all_central_stats_all_models <- rbind(all_central_stats_all_models, 
                                        df_central_stat_that_model)
}

dfs_split_by_modification <- split(all_central_stats_all_models, 
                          all_central_stats_all_models$modification)

# first thing: do this for the unmodified and then remove it from the list 
unmodified_df <- dfs_split_by_modification$UNM
dfs_split_by_modification$UNM <- NULL
unm_median_df <- unmodified_df %>% 
  group_by(central_coord, reference, model) %>% 
  summarise(median_three_rep = median(central_stat))

# consider one modification by one, then merge all of them in the same df
# do the facet boxplot
all_deltas_all_modifications <- NULL
for (modification in names(dfs_split_by_modification)) {
  # get the df corresponding to that modification 
  one_mod_df <- dfs_split_by_modification[[modification]]
  median_df <- one_mod_df %>% 
    group_by(central_coord, reference, model) %>% 
    summarise(median_three_rep = median(central_stat))
  merged_df <- inner_join(median_df, unm_median_df, 
                          by = c("reference", "central_coord", "model")) %>%
    mutate(delta_stat = median_three_rep.x - median_three_rep.y) %>%
    select(reference, central_coord, delta_stat, model) 
  merged_df$modification <- rep(modification, nrow(merged_df))

  # let's see if we can do the plotting 
  all_deltas_all_modifications <- rbind(all_deltas_all_modifications, 
                                        merged_df)
}

# plot, code from gregor: 
# https://github.com/novoalab/basecalling_models/blob/main/scripts/R/in_vitro/per_kmer_boxplots_pos0.R

t <- theme(
  legend.title = element_blank(),
  legend.text = element_text( size = 15),
  legend.position = "right",
  axis.text = element_text(size=15),
  axis.title=element_text(size=15),
  panel.grid.major.y = element_line(color = "grey90"),
  panel.grid.minor.y = element_line(color = "grey90"),
  panel.border = element_rect(colour = "black", fill = NA),
  strip.text.x = element_text(size = 15, face = "bold")
)

# will later be used for the stat_compare_means
pvalues_comparisons <- list( c("fast", "hac"), c("hac", "sup"), c("fast", "sup") )

p <- ggplot(all_deltas_all_modifications, 
       aes(x = model, y = delta_stat, fill = model)) + 
  geom_boxplot(aes(fill = model), 
               width=0.4, color="grey20",position = position_dodge(width =0.85),
               notch = TRUE, alpha = 0.9, lwd = 0.4,
               outlier.size = 0.25) + 
  geom_violin(aes(fill = model), position=position_dodge(width =0.85),
              alpha = 0.5, show.legend = F, color = NA, scale = "width", width = 0.7) +
  scale_y_continuous() + 
  theme_minimal() +
  facet_wrap(~modification, nrow = 1) + 
  theme_pubr() + t +  labs(x="",y=expression(Delta*"SumErr")) + 
  theme(
    axis.text = element_text(size=10),
    axis.title = element_text(size=20),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.justification = "center",
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    strip.text = element_text(size = 20),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.spacing = unit(1.5, "lines")
  )  + stat_n_text(size = 4, angle = 90)
p <- p + stat_compare_means(comparisons = pvalues_comparisons, label = "p.signif")

# and save
ggsave(filename = file_name,
       plot = p,
       path = output_dir,
       width = 15,
       height = 7.5,
       units = "in", 
       create.dir = T)

