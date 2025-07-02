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
  unite(col = overlap,c("V1", "V2"),  sep = "_", remove = F) %>% 
  # Remove single FP Am
  filter(!(V1 == "28s" & V2 == 1880 & V4 == "Am"))

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

dir_human <- setwd(here("..","__data","rRNA","modkit","m5C","C_only","human"))

dir_mouse <- setwd(here("..","__data","rRNA","modkit","m5C","C_only","mouse"))

dir_ecoli <- setwd(here("..","__data","rRNA","modkit","m5C","C_only","ecoli"))


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
    #change m to m5C to match ref annotation
    mutate(X4 = if_else(X4 == "m", "m5C", X4)) %>% 
    separate(sample, into = c("condition","replicate")) %>% 
    filter(X5 > 500) %>%
    select("chr" = X1, "start" = X2, "end" = X3, "Nval" = X5, "percent_mod" = X11, condition, replicate) %>% 
    mutate(chr = as.character(chr)) %>%
    unite(col = overlap,c("chr", "start"), sep = "_", remove = F) %>% 
    left_join(ref, by = c("overlap" = "overlap")) %>% 
    select(-V1, -V2, -V3, -overlap, -Nval, "modification" = V4) %>% 
    pivot_wider(names_from = condition, values_from = percent_mod, values_fill = 0) %>%
    
    mutate(mod_simplified = case_when(
      modification == "m5C" ~ "m5C",
      is.na(modification) ~ "unmodified",
      TRUE ~ "C_mod_non_m5C"
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



###############################################
### Boxplots of average modification stoich ###
###############################################


library(EnvStats)

mouse_filt %>% 
  ### Filter main rRNAs ####
  filter(chr %in% c("28", "18")) %>% 
  ### Filter for rep1 ####
  filter(replicate == "rep1") %>%
  
  ggplot(aes(x = reorder(modification, -wt), y = wt, fill = modification)) + 
  #scale_fill_manual(values = c("m5C" = "#f9d14a", "C_mod_non_m5C" = "#32CD32", "mod_adjacent" = "#C71585", "unmodified" = "grey60")) +
  geom_boxplot() +
  coord_flip() + 
  stat_n_text() +
  labs(y= "precent modified [%]",
       x = "") +
  theme_pubr() + t



ggsave(
  filename = here::here("..", "__results_Gregor", "rRNA", "wt_vs_ivt", "m5C", "boxplots_wt_amounts", "mouse_28s_18s.pdf"),
  plot = last_plot(),
  device = "pdf",
  width = 6,
  height = 8,
  units = "in"
)
