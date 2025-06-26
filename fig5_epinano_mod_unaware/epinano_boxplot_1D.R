# BOXPLOT "1D" --> plot the stat for the different kind of modifications. 
# 1. same model -> one plot 
# 2. same statistics, not further "modified" (no comp with the WT)
# one box for each modification  

# example of calling 
# Rscript boxplot_1D.R --working_dir epinano_kmers 
# --output_dir "./plots/epinano" 
# --file_name epinano_1E.pdf

library(ggplot2)
library(ggpubr)
library(argparse)

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
# go model by model
# create a df to update each time 
all_medians_all_models <- NULL
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

  all_medians <- lapply(df_file_mod$full_file_name, function(file_name) {
    read.table(file_name, header = T)$median})
  df_medians_that_model <- data.frame(
    median = unlist(all_medians),
    barcode = rep(df_file_mod$barcode, times = sapply(all_medians, length)),
    modification = rep(df_file_mod$modification, times = sapply(all_medians, length)), 
    model = rep(model, sum(sapply(all_medians, length))))
  
  # update the df with all the models 
  all_medians_all_models <- rbind(all_medians_all_models, df_medians_that_model)
}

# define the theme - taken from gregor 
# https://github.com/novoalab/basecalling_models/blob/main/scripts/R/in_vitro/error_plots_v2.R
t <- theme(
  legend.title = element_blank(),
  legend.text = element_text( size = 25),
  legend.position = "top",
  
  axis.text = element_text(size=20),
  axis.title=element_text(size=20),
  axis.text.x = element_text(size = 20,face = "bold"),
  
  panel.grid.major.y = element_line(color = "grey90"),
  
  panel.border = element_rect(colour = "black", fill = NA),
  strip.text.x = element_text(size = 15, face = "bold")
)  

# for the order of the boxes
all_medians_all_models$modification <- factor(
  all_medians_all_models$modification,
  levels = c("UNM", sort(setdiff(unique(all_medians_all_models$modification), "UNM")))
)


# and plot 
ggplot(all_medians_all_models, aes(x = modification, y = median, fill = modification)) +
  geom_boxplot(aes(fill = modification), 
               width=0.4, color="grey20",position = position_dodge(width =0.85),
               notch = TRUE, alpha = 0.9, lwd = 0.4,
               outlier.size = 0.25) +
  geom_violin(aes(fill = modification), position=position_dodge(width =0.85),
              alpha = 0.5, show.legend = F, color = NA, scale = "width", width = 0.7) +
  scale_y_continuous() +
  scale_fill_manual("modification", 
                   values = c('UNM' = 'grey50','m6A' = '#381a61', 'm5C' = '#f9d14a',
                              '5hmC' = "#ab3329", 'ac4C' = "#88a0dc",'Y' = '#e78429',
                              'm1Y' = '#ed968c','m5U' = '#7c4b73')) +
  labs(x="",y="Summed Error (absolute)") +
  theme_minimal() +
  facet_wrap(~model) + 
  theme_pubr() + t +
  theme(
    axis.text = element_text(size=6),
    axis.title = element_text(size=12),
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
  )

# and save
ggsave(filename = file_name,
       plot = last_plot(),
       path = output_dir,
       width = 10,
       height = 5,
       units = "in", 
       create.dir = T)
