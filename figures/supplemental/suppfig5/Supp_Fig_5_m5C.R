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




############################################
### Generate scatterplots of correlation ###
############################################

plot_correlation_scatter <- function(df, condition, rrna, labelx, labely, filename, modification){
  
  # Convert string to symbol for tidy evaluation
  condition_sym <- rlang::sym(condition)
  
  df_plot <- df %>% 
    select(chr, start, end, replicate, !!condition_sym) %>% 
    tidyr::pivot_wider(names_from = replicate, values_from = !!condition_sym, values_fill = 0)
  
  # Ensure rep1 and rep2 exist after pivot_wider
  if (!all(c("rep1", "rep2") %in% names(df_plot))) {
    stop("Columns 'rep1' and 'rep2' must exist after pivot_wider. Check your replicate values.")
  }
  
  # Compute Pearson correlation
  cor_val <- df_plot %>%
    filter(chr %in% rrna) %>%
    summarise(r = cor(rep1, rep2, method = "pearson")) %>%
    pull(r)
  
  # Create label with italic r
  corr_label <- tibble(
    x = 5,
    y = 90,
    label = paste0("italic(r) == ", round(cor_val, 2))
  )
  
  # Generate plot object
  p <- df_plot %>%
    filter(chr %in% rrna) %>%
    ggplot(aes(x = rep1, y = rep2)) +
    geom_point(size = 4, alpha = 0.75) +
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey50") +
    xlim(0, 100) +
    ylim(0, 100) +
    theme_minimal() +
    geom_text(data = corr_label, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, hjust = 0, size = 6, parse = TRUE) +
    labs(x = labelx, y = labely) +
    ggpubr::theme_pubr() + t
  
  # Save the plot using the provided modification name
  ggsave(
    filename = here::here("..", "__results_Gregor", "rRNA", "wt_vs_ivt", modification, "correlation_scatter", filename),
    plot = p,
    device = "pdf",
    width = 4,
    height = 4,
    units = "in"
  )
}



###### HUMAN pseU ######

plot_correlation_scatter(
  df = human_filt,
  condition = "ivt",
  rrna = c("28s", "18s"),
  labelx = "IVT rep1 [%]",
  labely = "IVT rep2 [%]",
  filename = "human_28s_18s_IVT.pdf",
  modification = "m5C"  
)


plot_correlation_scatter(
  df = human_filt,
  condition = "wt",
  rrna = c("28s", "18s"),
  labelx = "WT rep1 [%]",
  labely = "WT rep2 [%]",
  filename = "human_28s_18s_WT.pdf",
  modification = "m5C"  
)


###### Mouse pseU ######

plot_correlation_scatter(
  df = mouse_filt,
  condition = "ivt",
  rrna = c("28", "18"),
  labelx = "IVT rep1 [%]",
  labely = "IVT rep2 [%]",
  filename = "mouse_28s_18s_IVT.pdf",
  modification = "m5C"  
)


plot_correlation_scatter(
  df = mouse_filt,
  condition = "wt",
  rrna = c("28", "18"),
  labelx = "WT rep1 [%]",
  labely = "WT rep2 [%]",
  filename = "mouse_28s_18s_WT.pdf",
  modification = "m5C"  
)


###### Ecoli pseU ######

plot_correlation_scatter(
  df = ecoli_filt,
  condition = "ivt",
  rrna = c("23", "16"),
  labelx = "IVT rep1 [%]",
  labely = "IVT rep2 [%]",
  filename = "ecoli_23s_16s_IVT.pdf",
  modification = "m5C"  
)


plot_correlation_scatter(
  df = ecoli_filt,
  condition = "wt",
  rrna = c("23", "16"),
  labelx = "WT rep1 [%]",
  labely = "WT rep2 [%]",
  filename = "ecoli_23s_16s_WT.pdf",
  modification = "m5C"  
)


