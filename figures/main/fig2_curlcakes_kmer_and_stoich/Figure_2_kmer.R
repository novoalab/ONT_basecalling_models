###########################################################
## Plotting Densities for different kmer compositionts   ##
###########################################################

# Task 0: Load packages

pkgs <- c("openintro","tidyverse","backports","here","skimr","dplyr", "ggplot2", "ggsci","ggforce",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr","scales", "rstatix", "ggdist")

lapply(pkgs, library, character.only = TRUE)

##############################
# m6A_DRACH and inosine_m6A ##
##############################

#####################################
# Task 1: Import and annotate Data  #
#####################################


m6A_DRACH_100 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m6A_DRACH", "pod5---bc_31_s.txt"), header = T, sep = "\t") 

inosine_m6A_100 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_inosine_m6A", "pod5---bc_31_s.txt"), header = T, sep = "\t")

Y_pseU_100 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_pseU", "pod5---bc_46_s.txt"), header = T, sep = "\t")

m5C_m5C_100 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m5C", "pod5---bc_13_s.txt"), header = T, sep = "\t")


filter_and_annotate <- function(df, mod_code, is_true, central_kmer, minus_2_C, minus_1_C){
  
  df <- df %>%
    # Remove potential "-" strand predictions
    filter(ref_mod_strand == "+") %>%
    # Filter for desired mod code
    filter(mod_code == !!mod_code) %>% 
    # Filter out baseQ >= 10
    filter(base_qual >= 10) %>% 
    # Add ground-truth column
    mutate(is_true = is_true) %>% 
    # Annotate single central modified 5mers
    mutate(is_single_central_5mer = case_when(
      str_detect(query_kmer, central_kmer) ~ 1,
      T ~ 0
    )) %>% 
    # report only single A 5mers
    filter(is_single_central_5mer == 1) %>% 
    # assign flag for lowly and highly modified sites
    mutate(which_mod_prob = case_when(
      mod_qual >= 0.75 ~ 1,
      T ~ 0
    )) %>% 
    # Group by kmer and calculate median mod probability
    group_by(query_kmer) %>% 
    mutate(n = n(),# Number of observations per group
           median_mod_prob_per_kmer = median(mod_qual),
           median_base_qual_per_kmer = median(base_qual)) %>% 
    ungroup() %>% 
    # filter out any incomplete 5kmers ("-" containing)
    filter(!grepl("-", query_kmer)) %>% 
    mutate(is_leading_C = case_when(
      str_detect(query_kmer, minus_2_C) ~ 1,
      str_detect(query_kmer, minus_1_C) ~ 1,
      T ~ 0
    )) %>% 
    # add column with padding
    mutate(padded_kmer = paste0("NNN", query_kmer, "NNN"))
  
  # Modify x-axis labels to include 'n'
  df$query_kmer_n <- paste0(df$query_kmer, " [n=", df$n, "]")
  
  # Factorize leading C column
  
  df$is_leading_C <- as.factor(df$is_leading_C)
  
  return(df)
    
}


m6A_DRACH_100_filtered <- filter_and_annotate(df = m6A_DRACH_100, mod_code = "a", is_true = 1, central_kmer = "^[GTC][GTC][A][GTC][GTC]$", minus_2_C ="^[C][GTC][T][GTC][GTC]$", minus_1_C ="^[GTC][C][T][GTC][GTC]$")

inosine_m6A_100_filtered <- filter_and_annotate(df = inosine_m6A_100, mod_code = "a", is_true = 1, central_kmer = "^[GTC][GTC][A][GTC][GTC]$", minus_2_C ="^[C][GTC][A][GTC][GTC]$", minus_1_C ="^[GTC][C][A][GTC][GTC]$")

Y_pseU_100_filtered <- filter_and_annotate(df = Y_pseU_100, mod_code = "17802", is_true = 1, central_kmer = "^[GAC][GAC][T][GAC][GAC]$", minus_2_C ="^[C][GAC][T][GAC][GAC]$", minus_1_C ="^[GAC][C][T][GAC][GAC]$") 

m5C_m5C_100_filtered <- filter_and_annotate(df = m5C_m5C_100, mod_code = "m", is_true = 1, central_kmer = "^[GTA][GTA][C][GTA][GTA]$", minus_2_C ="^[C][GAC][C][GAC][GAC]$", minus_1_C ="^[GAC][C][C][GAC][GAC]$")


###########################################
## Plot modification probability per kmer #
###########################################

## Custom theme

t <- theme(
  legend.title = element_blank(),
  legend.text = element_text( size = 15),
  legend.position = "none",
  axis.text = element_text(size=15),
  axis.title=element_text(size=15),
  axis.text.x = element_text(angle = 90, hjust = 1),
  panel.grid.major.y = element_line(color = "grey90"),
  panel.grid.minor.y = element_line(color = "grey90"),
  panel.border = element_rect(colour = "black", fill = NA),
  strip.text.x = element_text(size = 15, face = "bold")
)


## Plot modification probability per kmer

list_of_plots  <- list(c(inosine_m6A_100_filtered,'inosine_m6A'),
                       c(Y_pseU_100_filtered,'pseU'),
                       c(m5C_m5C_100_filtered,'m5C'),
                       c(m6A_DRACH_100_filtered,'m6A_DRACH'))



plot_kmer_boxplots <- function(list_of_plots){
  for (i in 1:length(list_of_plots)) {
   
    df <- as.data.frame(list_of_plots[[i]])
    
    df %>%
      ggplot(aes(x = reorder(query_kmer_n, +median_mod_qual_per_kmer), y = mod_qual, fill = is_leading_C)) +
      stat_boxplot(geom = "errorbar",
                   lwd = 0.4,
                   width = 0.6) +
      # Central boxplot
      geom_boxplot(aes(fill = is_leading_C),
                   width=0.6, color="grey20",
                   alpha = 1, lwd = 0.4, notch = F, outlier.shape = NA) +
      scale_fill_manual(values = c('0' = 'grey80', '1' = '#AFAED7')) +
      labs(x= "", y = "Modification Probability") +
      theme_pubr() + t
    
    
    ggsave(filename = here("..","__results_Gregor","curlcakes","per_read","hac","per_kmer_analysis",paste0(list_of_plots[[i]][26],"_boxplots_per_kmer_for_baseQ_baseQ_greater_10_is_true_1.pdf")),       
           plot = last_plot(), 
           device = "pdf", 
           width = 15, height = 7.5, units = "in")
     
  }
}

plot_kmer_boxplots(list_of_plots = list_of_plots)


######################################
## Summarize to get median per kmer ##
######################################


calculate_summary <- function(df, name){
  
  df <- df %>% 
    group_by(query_kmer) %>% 
    summarise(median_kmer = median(mod_qual)) %>% 
    mutate(prob_group = case_when(
      median_kmer <= 0.25 ~ "low",
      median_kmer > 0.25 & median_kmer <= 0.50 ~ "medium low",
      median_kmer > 0.50 & median_kmer <= 0.75 ~ "medium high",
      median_kmer > 0.75 ~ "high"
    )) %>% 
    ungroup() %>% 
    count(prob_group) %>% 
    mutate(ratio = n / sum(n) * 100) %>% 
    mutate(model = name)
  
  df$prob_group <- factor(df$prob_group, levels = c("high","medium high","medium low","low"))
  
  return(df)
  
}



m6A_DRACH_summary <- calculate_summary(df = m6A_DRACH_100_filtered, name = 'm6A_DRACH')

inosine_m6A_summary <- calculate_summary(df = inosine_m6A_100_filtered, name = 'inosine_m6A')

pseU_summary <- calculate_summary(df = Y_pseU_100_filtered, name = 'pseU')

m5C_summary <- calculate_summary(df = m5C_m5C_100_filtered, name = 'm5C')

master_df <- rbind (m6A_DRACH_summary, inosine_m6A_summary, pseU_summary, m5C_summary)

master_df$model <- factor(master_df$model, levels = c('m6A_DRACH','inosine_m6A','pseU','m5C'))


 master_df%>% 
  ggplot(aes(x = model, y = ratio, fill = prob_group)) +
  geom_col() +
  geom_text(aes(label = paste0(round(ratio, 1), "%")),  # Add percentage labels
             position = position_stack(vjust = 0.5),  # Center labels in each section
             color = "black", size = 5) +  # Adjust text color and size
  scale_fill_manual(values = c('high' = '#7D82B8','medium high' = '#9FA3D0','medium low' = '#C1C4E5','low' = '#E3E5F5')) +
  labs(x = '', y = 'relative contribution [%]') +
  theme_pubr() + t +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = 'bottom'
  )

ggsave(filename = here("..","__results_Gregor","curlcakes","per_read","hac","enriched_kmers","relative_contribution_barplot.pdf"),       
       plot = last_plot(), 
       device = "pdf", 
       width = 6, height = 8, units = "in")

#########################################################
# Denisty plots of baseQ and modQ for identified kmers ##
#########################################################




###############
# inosine_m6A #
###############

inosine_m6A_100_filtered_density <- inosine_m6A_100_filtered %>%
  mutate(in_low = case_when(
    str_detect(query_kmer, "^[CT][C][A][T][G]$") ~ 1,
    T ~ 0
  ))

inosine_m6A_100_filtered_density$in_low <- as.factor(inosine_m6A_100_filtered_density$in_low)

inosine_m6A_100_filtered_density %>% count(in_low)

inosine_m6A_100_filtered_density %>% 
  ggplot(aes(x = mod_qual, fill = in_low)) +
  geom_density(aes(y = ..scaled..), alpha = 0.5) + 
  theme_pubr() + t +
  theme(
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.minor.x = element_line(color = "grey90"))


ggsave(filename = here("..","__results_Gregor","curlcakes","motifs","inosine_m6A_m6A_100","YCAUG_mod_qual.pdf"),       
       plot = last_plot(), 
       device = "pdf", 
       width = 8, height = 5, units = "in")


###############
# pseU #
###############

Y_pseU_100_filtered_density <- Y_pseU_100_filtered %>%
  mutate(in_low = case_when(
    str_detect(query_kmer, "^[C][CG][T][G][AC]$") ~ 1,
    T ~ 0
  ))

Y_pseU_100_filtered_density$in_low <- as.factor(Y_pseU_100_filtered_density$in_low)

Y_pseU_100_filtered_density %>% count(in_low)

Y_pseU_100_filtered_density %>% 
  ggplot(aes(x = mod_qual, fill = in_low)) +
  geom_density(aes(y = ..scaled..), alpha = 0.5) + 
  theme_pubr() + t +
  theme(
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.minor.x = element_line(color = "grey90"))


ggsave(filename = here("..","__results_Gregor","curlcakes","motifs","pseU_Y_100","CSUGM_mod_qual.pdf"),       
       plot = last_plot(), 
       device = "pdf", 
       width = 8, height = 5, units = "in")

#######
# m5C #
#######

m5C_m5C_100_filtered_density <- m5C_m5C_100_filtered %>%
  mutate(in_low = case_when(
    str_detect(query_kmer, "^[AT][A][C][AT][GT]$") ~ 1,
    T ~ 0
  ))

m5C_m5C_100_filtered_density$in_low <- as.factor(m5C_m5C_100_filtered_density$in_low)

m5C_m5C_100_filtered_density %>% count(in_low)

m5C_m5C_100_filtered_density %>% 
  ggplot(aes(x = mod_qual, fill = in_low)) +
  geom_density(aes(y = ..scaled..), alpha = 0.5) + 
  theme_pubr() + t +
  theme(
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.minor.x = element_line(color = "grey90"))


ggsave(filename = here("..","__results_Gregor","curlcakes","motifs","m5C_m5C_100","WACWK_mod_qual.pdf"),       
       plot = last_plot(), 
       device = "pdf", 
       width = 8, height = 5, units = "in")



#####################
# Kmer analysis UNM #
#####################

UNM_m6A_DRACH <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m6A_DRACH", "pod5---bc_1_s.txt"), header = T, sep = "\t") 

UNM_inosine_m6A <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_inosine_m6A", "pod5---bc_1_s.txt"), header = T, sep = "\t") 

UNM_pseU <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_pseU", "pod5---bc_1_s.txt"), header = T, sep = "\t") 

UNM_m5C <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m5C", "pod5---bc_1_s.txt"), header = T, sep = "\t") 


m6A_DRACH_UNM_fitlered <- filter_and_annotate(df = UNM_m6A_DRACH, mod_code = "a", is_true = 1, central_kmer = "^[GTC][GTC][A][GTC][GTC]$", minus_2_C ="^[C][GTC][T][GTC][GTC]$", minus_1_C ="^[GTC][C][T][GTC][GTC]$")

inosine_m6A_UNM_filtered <- filter_and_annotate(df = UNM_inosine_m6A, mod_code = "a", is_true = 1, central_kmer = "^[GTC][GTC][A][GTC][GTC]$", minus_2_C ="^[C][GTC][A][GTC][GTC]$", minus_1_C ="^[GTC][C][A][GTC][GTC]$")

pseU_UNM_filtered <- filter_and_annotate(df = UNM_pseU, mod_code = "17802", is_true = 1, central_kmer = "^[GAC][GAC][T][GAC][GAC]$", minus_2_C ="^[C][GAC][T][GAC][GAC]$", minus_1_C ="^[GAC][C][T][GAC][GAC]$") 

m5C_UNM_filtered <- filter_and_annotate(df = UNM_m5C, mod_code = "m", is_true = 1, central_kmer = "^[GTA][GTA][C][GTA][GTA]$", minus_2_C ="^[C][GAC][C][GAC][GAC]$", minus_1_C ="^[GAC][C][C][GAC][GAC]$")

###################################
## Save tables for motif analysis #
###################################

#################
## inosine_m6A ##
#################

inosine_m6A_UNM_filtered_high <- inosine_m6A_UNM_filtered %>% 
  # Filter for FP
  filter(mod_qual >= 0.05) %>% 
  select(padded_kmer)

inosine_m6A_UNM_filtered_low <- inosine_m6A_UNM_filtered %>% 
  # Filter for FP
  filter(mod_qual < 0.05) %>% 
  select(padded_kmer)

write.table(inosine_m6A_UNM_filtered_high, here("..", "__data","curlcakes", "processed","modkit_extract_full_for_motif", "0%_modified", "inosine_m6A_UNM_higher_0.05.tsv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(inosine_m6A_UNM_filtered_low, here("..", "__data","curlcakes", "processed","modkit_extract_full_for_motif", "0%_modified", "inosine_m6A_UNM_lower_0.05.tsv"), row.names = FALSE, col.names = FALSE, quote = FALSE)

#########
## m5C ##
#########

m5C_UNM_filtered_high <- m5C_UNM_filtered %>% 
  # Filter for FP
  filter(mod_qual >= 0.05) %>% 
  select(padded_kmer)

m5C_UNM_filtered_low <- m5C_UNM_filtered %>% 
  # Filter for FP
  filter(mod_qual < 0.05) %>% 
  select(padded_kmer)

write.table(m5C_UNM_filtered_high, here("..", "__data","curlcakes", "processed","modkit_extract_full_for_motif", "0%_modified", "m5C_UNM_higher_0.05.tsv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(m5C_UNM_filtered_low, here("..", "__data","curlcakes", "processed","modkit_extract_full_for_motif", "0%_modified", "m5C_UNM_lower_0.05.tsv"), row.names = FALSE, col.names = FALSE, quote = FALSE)

##########
## pseU ##
##########

pseU_UNM_filtered_high <- pseU_UNM_filtered %>% 
  # Filter for FP
  filter(mod_qual >= 0.05) %>% 
  select(padded_kmer)

pseU_UNM_filtered_low <- pseU_UNM_filtered %>% 
  # Filter for FP
  filter(mod_qual < 0.05) %>% 
  select(padded_kmer)

write.table(pseU_UNM_filtered_high, here("..", "__data","curlcakes", "processed","modkit_extract_full_for_motif", "0%_modified", "pseU_UNM_higher_0.05.tsv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(pseU_UNM_filtered_low, here("..", "__data","curlcakes", "processed","modkit_extract_full_for_motif", "0%_modified", "pseU_UNM_lower_0.05.tsv"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  

#######################
## Plot densities UNM #
#######################

pseU_UNM_filtered %>% 
  # Filter for FP
  #filter(mod_qual > 0.05) %>% 
  # Filter for single central mod
  mutate(is_target_kmer = factor(case_when(
    str_detect(query_kmer,"^[C][A][T][A][G]$") ~ 1,
    T ~ 0
  ))) %>% 
  ggplot(aes(x = mod_qual, fill = is_target_kmer)) +
  geom_density(aes( y = ..scaled..), alpha = 0.5) +
  theme_pubr() + t +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.major.x = element_line(color = "grey90"),
    panel.grid.minor.x = element_line(color = "grey90"))

ggsave(filename = here("..","__results_Gregor","curlcakes","motifs","0_modified", "pseU_Y_100","CAUAG_mod_qual.pdf"),       
       plot = last_plot(), 
       device = "pdf", 
       width = 8, height = 5, units = "in")







