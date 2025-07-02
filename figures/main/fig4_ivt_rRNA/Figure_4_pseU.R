################################################
## Plotting density curves of Stoichiometry  ##
################################################

# Task 0: Load packages

pkgs <- c("openintro","tidyverse","backports","here","skimr","dplyr", "ggplot2", "ggsci","ggforce",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr","scales")

lapply(pkgs, library, character.only = TRUE)

#######################
# Task 1: Import Data #
#######################

# Task 1.1: Import ref

ref_human <- read.csv(here("../__ref/rRNA_mods/human_all_rRNAmods.bed"), header = F, sep = "\t") %>% 
  # Create column to overlap on
  unite(col = overlap,c("V1", "V2"),  sep = "_", remove = F)

ref_mouse<- read.csv(here("../__ref/rRNA_mods/mouse_rrna_modification_annotation.bed"), header = F, sep = "\t") %>% 
  # remove S from annotation to fit .fa naming
  mutate(across(1, ~ str_replace_all(.x, "S", ""))) %>% 
  # Create column to overlap on
  unite(col = overlap,c("V1", "V2"),  sep = "_", remove = F)  
  

ref_ecoli <- read.csv(here("../__ref/rRNA_mods/rRNA_Ecoli_Mods.bed"), header = F, sep = "\t") %>% 
  # remove S from annotation to fit .fa naming
  mutate(across(1, ~ str_replace_all(.x, "S", ""))) %>% 
  # Create column to overlap on
  unite(col = overlap,c("V1", "V2"),  sep = "_", remove = F)

# Task 1.2: Import modkit output

dir_human <- setwd(here("..","__data","rRNA","modkit","pseU","U_only","human"))

dir_mouse <- setwd(here("..","__data","rRNA","modkit","pseU","U_only","mouse"))

dir_ecoli <- setwd(here("..","__data","rRNA","modkit","pseU","U_only","ecoli"))


file_list_human <- list.files(path = dir_human)

file_list_mouse <- list.files(path = dir_mouse)

file_list_ecoli <- list.files(path = dir_ecoli)

# Task 1.3: Function for reading in all files from a directory 

# Create a loop to read in every file of the directory and append it to the initialized data.frame plus add a new column that contains the name of the pool
# Comment: Could be sped up with fread and data.tables

read_dir <- function(file_list, work_dir){
  
  setwd(work_dir)
  
  dataset <- data.frame()
  
  for (i in 1:length(file_list)){
    temp_data <- read_tsv(file_list[i], col_names = F) #each file will be read in, specify which columns you need read in to avoid any errors # specifying col_types is essential to see spike_ins
    temp_data$sample <-  gsub("\\.bed", "", file_list[i])#clean the data as needed, in this case I am creating a new column that indicates which file each row of data came from
    dataset <- plyr::rbind.fill(dataset, temp_data) #for each iteration, bind the new data to the building dataset
  }
  
  rm(i)
  rm(temp_data)
  
  return(dataset)
}

human_raw <- read_dir(file_list = file_list_human, work_dir = dir_human)

mouse_raw <- read_dir(file_list = file_list_mouse, work_dir = dir_mouse)

ecoli_raw <- read_dir(file_list = file_list_ecoli, work_dir = dir_ecoli)


####################################
# Task 2: Prefiltering of raw data #
####################################


filter_table <- function(df, ref){
  # Precompute mod_adjacent positions from ref
  ref_mod_pos <- ref %>% 
    filter(!(V4 %in% c("unmodified", NA))) %>% 
    transmute(chr = as.character(V1), start = V2)
  
  ref_mod_adj <- ref_mod_pos %>%
    mutate(start_minus1 = start - 1, start_plus1 = start + 1) %>% 
    select(chr, start_minus1, start_plus1) %>% 
    pivot_longer(cols = c(start_minus1, start_plus1), names_to = "offset", values_to = "adj_start") %>% 
    select(chr, start = adj_start) %>% 
    distinct()
  
  df_filt <- df %>% 
    separate(sample, into = c("condition","replicate")) %>% 
    filter(X5 > 500) %>%
    select("chr" = X1, "start" = X2, "end" = X3, "Nval" = X5, "percent_mod" = X11, condition, replicate) %>% 
    mutate(chr = as.character(chr)) %>%
    unite(col = overlap,c("chr", "start"), sep = "_", remove = F) %>% 
    left_join(ref, by = c("overlap" = "overlap")) %>% 
    select(-V1, -V2, -V3, -overlap, -Nval, "modification" = V4) %>% 
    pivot_wider(names_from = condition, values_from = percent_mod, values_fill = 0) %>%
    
    mutate(mod_simplified = case_when(
      modification %in% c('Ψ', 'Y') ~ "Ψ",
      is.na(modification) ~ "unmodified",
      TRUE ~ "U_mods_nonΨ"
    )) %>%
    
    mutate(
      wt_threshold = if_else(wt > 10, 1, 0),
      ivt_threshold = if_else(ivt < 5, 1, 0)
    ) %>%
    
    left_join(ref_mod_adj %>% mutate(is_mod_adj = TRUE), by = c("chr", "start")) %>% 
    mutate(mod_simplified = case_when(
      mod_simplified == "unmodified" & is_mod_adj ~ "mod_adjacent",
      TRUE ~ mod_simplified
    )) %>% 
    select(-is_mod_adj)
  
  return(df_filt)
}


human_filt <- filter_table(df = human_raw, ref = ref_human)

mouse_filt <- filter_table(df = mouse_raw, ref = ref_mouse)

ecoli_filt <- filter_table(df = ecoli_raw, ref = ref_ecoli)


##########################
### Save unified table ###
##########################

human_filt$species <- "H. sapiens"
 
human_filt$species <-  "M. musculus"

ecoli_filt$species <-  "E. coli"


master_tbl <- rbind(human_filt, human_filt, ecoli_filt)

write.table(master_tbl, file= here("..","master_table_filtered_pseU.tsv"))

##########################
# Task 3: Plot WT vs IVT #
##########################


t <- theme(
  legend.title = element_blank(),
  legend.text = element_text( size = 15),
  legend.position = "bottom",
  axis.text = element_text(size=15),
  axis.title=element_text(size=15),
  panel.grid.major.x = element_line(color = "grey90"),
  panel.grid.minor.x = element_line(color = "grey90"),
  panel.grid.major.y = element_line(color = "grey90"),
  panel.grid.minor.y = element_line(color = "grey90"),
  panel.border = element_rect(colour = "black", fill = NA),
  strip.text.x = element_text(size = 15, face = "bold")
)


colors <- c("Ψ" = "#e78429","Um" = "#00429d","m3U" = "#6a5acd","m5U" = "#C71585","xp4U" = "#228B22",
            "Y" = "#e78429","D" = "#800020", "m3Y" = "#00CCCC",
            "Ψm" = "#32CD32", "Am" = "#FF6F61","Gm" = "#2e003e"
            )

# Scatterplot

ecoli_filt %>% 
filter(chr %in% c("23","16")) %>% 
ggplot(aes(x = ivt, y = wt, color = mod_simplified)) +
  geom_point(size = 4, alpha = 0.75) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey50") +
  
  # Add hline for for wt cutoff of greater 10%
  geom_hline(yintercept = 10, linetype = "dotted", color = "red") +
  # Add vline for ivt cutoff of greater 5%
  geom_vline(xintercept = 5, linetype = "dotted", color = "red") +
  
  scale_color_manual(values = c("Ψ" = "#e78429", "U_mods_nonΨ" = "#32CD32", "mod_adjacent" = "#C71585", "unmodified" = "grey60")) +
  xlim(0, 100) +
  ylim(0, 100) +
  theme_minimal() +
  # Prevent alpha in legends
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = "IVT rep1 [%]",
       y = "WT rep1 [%]") +
  theme_pubr() + t

ggsave(filename = here("..","__results_Gregor","rRNA","wt_vs_ivt", "pseU","scatterplots","ecoli_28S_18S_rep1.pdf"),       
       plot = last_plot(), 
       device = "pdf", 
       width = 6, height = 6, units = "in")

# Boxplot Quantity:
library(EnvStats)

ecoli_filt %>%
  filter(chr %in% c("23", "16")) %>%
  
  # remove sites introduced by join and used to viz in scatter
  filter(wt != 0 & ivt != 0) %>%
  
  ggplot(aes(y = wt, x = reorder(modification, -wt), fill = modification)) +
  geom_boxplot() +  # no outliers, whisker caps shown by default
  #scale_fill_manual(values = c("Ψ" = "#e78429", "U_mods_nonΨ" = "#32CD32", "mod_adjacent" = "#C71585", "unmodified" = "grey60")) +
  stat_n_text() +
  coord_flip() +
  theme_pubr() + t

ggsave(filename = here("..","__results_Gregor","rRNA","wt_vs_ivt", "pseU","boxlpots_wt_amounts","ecoli_23S_16S_rep1_simplified_all_mods.pdf"),       
       plot = last_plot(), 
       device = "pdf", 
       width = 6, height = 6, units = "in")


## Barplots relative abundance 


human_summary <- human_filt %>% 
  filter(wt_threshold == 1 & ivt_threshold == 1) %>% 
  # filter out 0s introduced by including ivt sites!
  filter(wt > 0 & ivt > 0) %>%  
  count(mod_simplified) %>% 
  mutate(proportion = (n / sum(n)),
         label = paste0(round(proportion * 100, 1), "%")) %>% 
  mutate(species = "human")

mouse_summary <- mouse_filt %>% 
  filter(wt_threshold == 1 & ivt_threshold == 1) %>%
  filter(wt > 0 & ivt > 0) %>%  
  count(mod_simplified) %>% 
  mutate(proportion = (n / sum(n)),
         label = paste0(round(proportion * 100, 1), "%")) %>% 
  mutate(species = "mouse")

ecoli_summary <- ecoli_filt %>% 
  filter(wt_threshold == 1 & ivt_threshold == 1) %>% 
  filter(wt > 0 & ivt > 0) %>%  
  count(mod_simplified) %>% 
  mutate(proportion = (n / sum(n)),
         label = paste0(round(proportion * 100, 1), "%")) %>% 
  mutate(species = "ecoli")

master_df <- rbind(human_summary, mouse_summary, ecoli_summary)

master_df$mod_simplified <- factor(master_df$mod_simplified, levels = c("unmodified", "mod_adjacent","U_mods_nonΨ","Ψ"))

master_df$species <- factor(master_df$species, levels = c("human", "mouse","ecoli"))


master_df %>% 
ggplot(aes(x = species, y = proportion, fill = mod_simplified)) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), color = "white", size = 4) +
  scale_fill_manual(values = c("Ψ" = "#e78429", "U_mods_nonΨ" = "#32CD32", "mod_adjacent" = "#C71585", "unmodified" = "grey60")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = NULL,
    y = "Proportion",
    fill = "Modkit",
    title = "Unfiltered Modkit [> 10% wt-cutoff | <5% ivt-cutoff]"
  ) + theme_pubr() + t


ggsave(filename = here("..","__results_Gregor","rRNA","wt_vs_ivt", "pseU","stacked_barchart","all_species_wt+ivt_filter_28S-23S_18S-16S_rep1_simplified.pdf"),       
       plot = last_plot(), 
       device = "pdf", 
       width = 6, height = 4, units = "in")












###############
## Test Zone ##
###############


# Create the label column for plotting
ecoli_filt <- ecoli_filt %>%
  mutate(label = paste(chr, start, sep = "_"))


# Plot
ecoli_filt %>%
  filter(chr %in% c("23", "16")) %>%
  ggplot(aes(x = ivt, y = wt, color = mod_simplified)) +
  geom_point(size = 4, alpha = 0.75) +
  
  # Label points with "chr_start"
  geom_text(aes(label = label), 
           hjust = -0.1, vjust = 0.5, size = 3, check_overlap = TRUE) +
  
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey50") +
  geom_hline(yintercept = 10, linetype = "dotted", color = "red") +
  geom_vline(xintercept = 5, linetype = "dotted", color = "red") +
  
  scale_color_manual(values = c("Ψ" = "#e78429", "U_mods_nonΨ" = "#32CD32", "mod_adjacent" = "#C71585", "unmodified" = "grey60")) +
  xlim(0, 100) +
  ylim(0, 100) +
  theme_minimal() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  labs(x = "IVT rep1 [%]",
       y = "WT rep1 [%]") +
  theme_pubr() + t
