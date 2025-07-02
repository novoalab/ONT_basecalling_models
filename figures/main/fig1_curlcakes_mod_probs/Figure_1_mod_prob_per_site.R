#########################################################
## Plotting Densities for per site and modification   ##
#########################################################

# Task 0: Load packages

pkgs <- c("openintro","tidyverse","backports","here","skimr","dplyr", "ggplot2", "ggsci","ggforce","viridis",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr","scales", "rstatix", "ggdist")

lapply(pkgs, library, character.only = TRUE)

##############################
# m6A_DRACH and inosine_m6A ##
##############################

#####################################
# Task 1: Import and annotate Data  #
#####################################


m6A_DRACH_m6A_100_rep1 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m6A_DRACH", "pod5---bc_31_s.txt"), header = T, sep = "\t") 

m6A_DRACH_m6A_100_rep1$mod_status <- 'm6A_100'

m6A_DRACH_UNM_rep1 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m6A_DRACH", "pod5---bc_1_s.txt"), header = T, sep = "\t") 

m6A_DRACH_UNM_rep1$mod_status <- 'UNM'

m6A_DRACH <- rbind(m6A_DRACH_m6A_100_rep1, m6A_DRACH_UNM_rep1)



inosine_m6A_m6A_100_rep1 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_inosine_m6A", "pod5---bc_31_s.txt"), header = T, sep = "\t")

inosine_m6A_m6A_100_rep1$mod_status <- 'm6A_100'

inosine_m6A_UNM_rep1 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_inosine_m6A", "pod5---bc_1_s.txt"), header = T, sep = "\t")

inosine_m6A_UNM_rep1$mod_status <- 'UNM'

inosine_m6A <- rbind(inosine_m6A_m6A_100_rep1, inosine_m6A_UNM_rep1)



pseU_Y_100 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_pseU", "pod5---bc_46_s.txt"), header = T, sep = "\t")

pseU_Y_100$mod_status <- 'Y_100'

pseU_UNM <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_pseU", "pod5---bc_1_s.txt"), header = T, sep = "\t")

pseU_UNM$mod_status <- 'UNM'

m1Y_pseU <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_pseU", "pod5---bc_10_s.txt"), header = T, sep = "\t")

m1Y_pseU$mod_status <- 'm1Y_100'

m5U_pseU <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241016_Pool3_rna004_130bps_hac@v5.1.0_pseU", "pod5---bc_49_s.txt"), header = T, sep = "\t")

m5U_pseU$mod_status <- 'm5U_100'


pseU <- rbind(pseU_Y_100, pseU_UNM, m1Y_pseU, m5U_pseU)



m5C_m5C_100 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m5C", "pod5---bc_13_s.txt"), header = T, sep = "\t")

m5C_m5C_100$mod_status <- 'm5C_100'

m5C_UNM <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m5C", "pod5---bc_1_s.txt"), header = T, sep = "\t")

m5C_UNM$mod_status <- 'UNM'

hm5C_m5C_model <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m5C", "pod5---bc_7_s.txt"), header = T, sep = "\t")

hm5C_m5C_model$mod_status <- "hm5C_100"

ac4C_m5C_model <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m5C", "pod5---bc_4_s.txt"), header = T, sep = "\t")

ac4C_m5C_model$mod_status <- "ac4C_100"

m5C <- rbind(m5C_m5C_100, m5C_UNM, hm5C_m5C_model, ac4C_m5C_model)


rm(m6A_DRACH_m6A_100_rep1, m6A_DRACH_UNM_rep1, inosine_m6A_m6A_100_rep1, inosine_m6A_UNM_rep1, pseU_Y_100, pseU_UNM, m5C_m5C_100, m5C_UNM)



filter_and_annotate <- function(df, mod_code, is_true, central_kmer, minus_2_C, minus_1_C){
  
  df <- df %>%
    # Remove potential "-" strand predictions
    filter(ref_mod_strand == "+") %>%
    # Filter out any incomplete 5mers (e.g: "GCAT-")
    filter(!grepl("-", query_kmer)) %>%
    # Filter for desired mod code
    filter(mod_code == !!mod_code) %>% 
    # Turn mod_qual into percent
    mutate(mod_qual = mod_qual * 100) %>%
    # Filter out baseQ >= 10
    filter(base_qual >= 10) %>% 
    # Add ground-truth column
    mutate(is_true = is_true) %>% 
    # Annotate single central modified 5mers
    mutate(is_single_central_5mer = case_when(
      str_detect(query_kmer, central_kmer) ~ 1,
      T ~ 0
    )) 
  
  return(df)
  
}


m6A_DRACH_filtered <- filter_and_annotate(df = m6A_DRACH, mod_code = "a", is_true = 1, central_kmer = "^[GT][G][A][C][CT]$")

inosine_m6A_filtered <- filter_and_annotate(df = inosine_m6A, mod_code = "a", is_true = 1, central_kmer = "^[GTC][GTC][A][GTC][GTC]$")

pseU_filtered <- filter_and_annotate(df = pseU, mod_code = "17802", is_true = 1, central_kmer = "^[GAC][GAC][T][GAC][GAC]$") 

m5C_filtered <- filter_and_annotate(df = m5C, mod_code = "m", is_true = 1, central_kmer = "^[GTA][GTA][C][GTA][GTA]$")


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
      scale_x_continuous(limits = c(0,100)) +
      # Add consitend mod coloring
      scale_fill_manual(values = fill) +
      # Add median mod and unmod text label
      annotate("text", x = median_mod, y = Inf, label = paste0("Median: ", round(median_mod, 1)),
               vjust = 2, hjust = 1.1, color = "black") +
      annotate("text", x = median_unm, y = Inf, label = paste0("Median: ", round(median_unm, 1)),
               vjust = 2, hjust = -0.1, color = "black") +
      labs(x = "Modification Probability [%]" ,y = "") +
      theme_pubr() + t
    
    
    
    ggsave(filename = here("..","__results_Gregor","curlcakes","per_read","hac", "density_per_position",paste0(name,"single_central_mod_density_per_read.pdf")),       
           plot = last_plot(), 
           device = "pdf", 
           width = 7, height = 6, units = "in") 
    
}


plot_density(df = m6A_DRACH_filtered, fill = c("m6A_100" ='#381a61', "UNM" = "grey75") , name = "m6A_DRACH", mod = 'm6A_100')

plot_density(df = inosine_m6A_filtered, fill = c("m6A_100" ='#381a61', "UNM" = "grey75") , name = "inosine_m6A", mod = 'm6A_100')

plot_density(df = pseU_filtered, fill = c("Y_100" ='#e78429', "m1Y_100" ='#EE9A91', "m5U_100" ='#82476F',"UNM" = "grey75") , name = "pseU", mod = 'Y_100')

plot_density(df = m5C_filtered, fill = c("m5C_100" ='#f9d14a', "hm5C_100" ='#AD2D1F', "ac4C_100" ='#96A6D5', "UNM" = "grey75") , name = "m5C", mod = 'm5C_100')












###########################################################
## Which kmers cause the bump in lowly mod prob regions ###
###########################################################

test <- inosine_m6A_filtered %>%
  filter(mod_qual > 0) %>% 
  group_by(mod_status) %>% 
  #filter(mod_qual >= 4.5 & mod_qual <= 7.5) %>% 
  count(query_kmer)
  #slice_max(n = 20, order_by = n)

pseU_filtered %>% 
  #filter(mod_status == "UNM") %>% 
  filter(query_kmer == 'CATCC') %>% 
  filter(mod_qual > 80) %>% 
  count(query_kmer)
  ggplot(aes(x = mod_qual, fill = mod_status)) +
  geom_density(aes(y = ..scaled..), alpha = 0.75) +
  theme_pubr() + t + 
  labs(x = 'Modification Probability [%]', y = '')


test %>% 
  filter(mod_status == 'UNM') %>% 
  ggplot(aes(x = reorder(query_kmer, +n), y = n)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "", y = "counts [n]") +
  coord_flip() +
  theme_pubr() + t 
  




##########################################################
### Check for overrepresented kmers in UNM > 0% modified #
##########################################################



plot_enrichment_score <- function(df, name){
  
  # Step 1: Compute total proportion of each k-mer in the full dataset
  kmer_total_counts <- df %>%
    filter(is_single_central_5mer == 1) %>% 
    filter(mod_status == 'UNM') %>% 
    group_by(query_kmer) %>%
    summarise(total_count = n(), .groups = "drop") %>%
    mutate(total_prop = total_count / sum(total_count))  # Convert to proportions
  
  # Step 2: Filter data where mod_qual > 0 and compute proportion in filtered dataset
  filtered_data <- df %>%
    filter(is_single_central_5mer == 1) %>% 
    filter(mod_status == 'UNM') %>% 
    filter(mod_qual > 0)  
  
  kmer_filtered_counts <- filtered_data %>%
    group_by(query_kmer) %>%
    summarise(filtered_count = n(), .groups = "drop") %>%
    mutate(filtered_prop = filtered_count / sum(filtered_count))  # Convert to proportions
  
  # Step 3: Merge both proportions and compute enrichment ratio
  enrichment <- kmer_total_counts %>%
    left_join(kmer_filtered_counts, by = "query_kmer") %>%
    #replace_na(list(filtered_count = 0, filtered_prop = 0)) %>%  # Fill missing values
    mutate(enrichment_ratio = filtered_prop / total_prop) %>% # Compute enrichment 
    slice_max(n=20 , order_by = enrichment_ratio)
  
  
  # Step 4: Plot enrichment ratio
  ggplot(enrichment, aes(x = reorder(query_kmer, enrichment_ratio), y = enrichment_ratio, fill = enrichment_ratio)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey75", size = 1) +  # Reference line at 1
    scale_fill_viridis_c(option = "viridis", direction = 1) +# Fill viridis
    labs(title = "",
         x = "",
         y = "Enrichment Score",
         fill = "Enrichment") +
    coord_flip() +
    theme_pubr() + t +
    theme(
      legend.text = element_text(size = 15, angle = 270)
    )
  
  ggsave(filename = here("..","__results_Gregor","curlcakes","per_read","hac", "enriched_kmers",paste0(name,"_central_mod_5mers_enriched_in_UNM_greater_0.pdf")),       
         plot = last_plot(), 
         device = "pdf", 
         width = 5, height = 7, units = "in")
  
}


plot_enrichment_score(df = inosine_m6A_filtered, name = 'inosine_m6A')

plot_enrichment_score(df = m6A_DRACH_filtered, name = 'm6A_DRACH')

plot_enrichment_score(df = pseU_filtered, name = 'pseU')

plot_enrichment_score(df = m5C_filtered, name = 'm5C')



inosine_m6A_filtered %>%
  #filter(mod_qual > 0.75) %>% 
  group_by(mod_status, read_id) %>%
  summarize(median_mod_qual = median(mod_qual)) %>% 
  ungroup() %>% 
  filter(median_mod_qual > 0.75, mod_status == 'UNM') %>% 
  count()
  ggplot(aes(x = mod_qual, fill = mod_status)) +
  geom_density(aes(y = ..scaled..))

  
  
test <- inosine_m6A_filtered %>% filter(query_kmer == 'CCACG') %>% filter(mod_status == "UNM") %>% group_by(read_id) %>% count()
  
  
