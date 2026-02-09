#########################################################
## Plotting Densities for per site and modification   ##
#########################################################

# Task 0: Load packages

pkgs <- c("openintro","tidyverse","backports","here","skimr","dplyr", "ggplot2", "ggsci","ggforce","viridis",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr","scales", "rstatix", "ggdist")

lapply(pkgs, library, character.only = TRUE)

#####################################
# Task 1: Import and annotate Data  #
#####################################


inosine_m6A <- read.table(file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","inosine_m6A_all_kmers.tsv"), header = T, sep = "\t") %>% 
  filter(ratio == 0.5) %>% 
  mutate(mod_status = case_when(
    is_true == 1 ~ "m6A_100",
    T ~ "UNM"
  ))

m6A_DRACH <- read.table(file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","m6A_DRACH_all_kmers.tsv"), header = T, sep = "\t") %>% 
  filter(ratio == 0.5) %>% 
  mutate(mod_status = case_when(
    is_true == 1 ~ "m6A_100",
    T ~ "UNM"
  ))

Y_pseU_all_kmers <- read.table(file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","pseU_all_kmers.tsv"), header = T, sep = "\t") %>% 
  filter(ratio == 0.5) %>% 
  mutate(mod_status = case_when(
    is_true == 1 ~ "Y_100",
    T ~ "UNM"
  ))

m1Y_pseU_all_kmers <- read.table(file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","m1Y_all_kmers.tsv"), header = T, sep = "\t") %>% 
  filter(ratio == 0.5) %>% 
  mutate(mod_status = case_when(
    is_true == 1 ~ "m1Y_100",
    T ~ "UNM"
  ))

m5U_pseU_all_kmers <- read.table(file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","m5U_all_kmers.tsv"), header = T, sep = "\t") %>% 
  filter(ratio == 0.5) %>% 
  mutate(mod_status = case_when(
    is_true == 1 ~ "m5U_100",
    T ~ "UNM"
  ))



m5C_m5C_all_kmers <- read.table(file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","m5C_all_kmers.tsv"), header = T, sep = "\t") %>% 
  filter(ratio == 0.5) %>% 
  mutate(mod_status = case_when(
    is_true == 1 ~ "m5C_100",
    T ~ "UNM"
  ))

hm5C_m5C_all_kmers <- read.table(file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","hm5C_model_all_kmers.tsv"), header = T, sep = "\t") %>% 
  filter(ratio == 0.5) %>% 
  mutate(mod_status = case_when(
    is_true == 1 ~ "hm5C_100",
    T ~ "UNM"
  ))

ac4C_m5C_all_kmers <- read.table(file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","ac4C_model_all_kmers.tsv"), header = T, sep = "\t") %>% 
  filter(ratio == 0.5) %>% 
  mutate(mod_status = case_when(
    is_true == 1 ~ "ac4C_100",
    T ~ "UNM"
  ))

#####################################################
### add columns describing is_single_central_5mer ###
#####################################################

annotate <- function(df, central_kmer){
  
  df <- df %>%
    # Annotate single central modified 5mers
    mutate(is_single_central_5mer = case_when(
      str_detect(query_kmer, central_kmer) ~ 1,
      T ~ 0
    )) 
  
  return(df)
  
}


m6A_DRACH_filtered <- annotate(df = m6A_DRACH, central_kmer = "^[GT][G][A][C][CT]$")

inosine_m6A_filtered <- annotate(df = inosine_m6A, central_kmer = "^[GTC][GTC][A][GTC][GTC]$")



Y_filtered <- annotate(df = Y_pseU_all_kmers, central_kmer = "^[GAC][GAC][T][GAC][GAC]$") 

m1Y_filtered <- annotate(df = m1Y_pseU_all_kmers, central_kmer = "^[GAC][GAC][T][GAC][GAC]$") 

m5U_filtered <- annotate(df = m5U_pseU_all_kmers, central_kmer = "^[GAC][GAC][T][GAC][GAC]$") 

pseU_filtered <- rbind(Y_filtered, m1Y_filtered, m5U_filtered)



m5C_filtered <- annotate(df = m5C_m5C_all_kmers, central_kmer = "^[GTA][GTA][C][GTA][GTA]$")

hm5C_filtered <- annotate(df = hm5C_m5C_all_kmers, central_kmer = "^[GTA][GTA][C][GTA][GTA]$")

ac4C_filtered <- annotate(df = ac4C_m5C_all_kmers, central_kmer = "^[GTA][GTA][C][GTA][GTA]$")

m5C_filtered <- rbind(m5C_filtered, hm5C_filtered, ac4C_filtered)


###########################################
## Plot modification probability per kmer #
###########################################

## Custom theme

t <- theme(
  legend.title = element_blank(),
  legend.text = element_text(size = 15),
  legend.position = "bottom",
  axis.text = element_text(size=15),
  axis.title=element_text(size=15),
  panel.grid.major.y = element_line(color = "grey90"),
  panel.grid.minor.y = element_line(color = "grey90"),
  panel.border = element_rect(colour = "black", fill = NA),
  strip.text.x = element_text(size = 15, face = "bold")
)


plot_density <- function(df,fill,name, mod){
  
  # Task 2.3: Compute median value
  df <- df %>% 
    filter(is_single_central_5mer == 1)
  
  median_mod <- df %>%
    filter(mod_status == mod) %>% 
    summarize(median_mod_qual = median(mod_qual)) %>%
    pull(median_mod_qual)
  
  median_unm <- df %>%
    filter(mod_status == 'UNM') %>%
    summarize(median_mod_qual = median(mod_qual)) %>%
    pull(median_mod_qual)
  
  
  # Task 2.4: Density Plot of predicted stoichiometries and associated median value
  
  
  df %>% 
    ggplot(aes(x = mod_qual, fill = mod_status)) +
    geom_density(aes(y = ..scaled..), alpha = 0.5, colour = NA) +
    # Add the median mod and unmod as a vertical line
    geom_vline(xintercept = median_mod, linetype = "dashed", color = "black", size = 0.75) +
    geom_vline(xintercept = median_unm, linetype = "dashed", color = "black", size = 0.75) +
    scale_x_continuous(limits = c(0,1)) +
    # Add consitend mod coloring
    scale_fill_manual(values = fill) +
    # Add median mod and unmod text label
    ggplot2::annotate("text", x = median_mod, y = Inf, label = paste0("Median: ", round(median_mod, 3)),
             vjust = 2, hjust = 1.1, color = "black") +
    ggplot2::annotate("text", x = median_unm, y = Inf, label = paste0("Median: ", round(median_unm, 3)),
             vjust = 2, hjust = -0.1, color = "black") +
    labs(x = "Modification Probability [%]" ,y = "") +
    theme_pubr() + t
  
  
  
  ggsave(filename = here("..","__results_Gregor","curlcakes","per_read","hac", "density_per_position", "tmp", paste0(name,"_all_kmer_central_mod_density_per_read.pdf")),       
         plot = last_plot(), 
         device = "pdf", 
         width = 7, height = 6, units = "in") 
  
}


plot_density(df = m6A_DRACH_filtered, fill = c("m6A_100" ='#381a61', "UNM" = "grey75") , name = "m6A_DRACH", mod = 'm6A_100')

plot_density(df = inosine_m6A_filtered, fill = c("m6A_100" ='#381a61', "UNM" = "grey75") , name = "inosine_m6A", mod = 'm6A_100')

plot_density(df = pseU_filtered, fill = c("Y_100" ='#e78429', "m1Y_100" ='#EE9A91', "m5U_100" ='#82476F',"UNM" = "grey75") , name = "pseU", mod = 'Y_100')

plot_density(df = m5C_filtered, fill = c("m5C_100" ='#f9d14a', "hm5C_100" ='#AD2D1F', "ac4C_100" ='#96A6D5', "UNM" = "grey75") , name = "m5C", mod = 'm5C_100')


