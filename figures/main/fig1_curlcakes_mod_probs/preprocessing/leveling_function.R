############################################################################################
## Leveling and annotation function to obtains different mixes of UNM/MOD CC1-4  ##
############################################################################################

# Task 0: Load packages

pkgs <- c("openintro","tidyverse","backports","here","skimr","dplyr", "ggplot2", "ggsci","ggforce",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr","scales")

lapply(pkgs, library, character.only = TRUE)

###############################
# Sampling Function for CC1-4 #
###############################


sample_by_ratios <- function(focal_df, sample_df, ratios) {
  
  if (any(ratios < 0 | ratios > 1)) {
    stop("All ratios must be between 0 and 1")
  }
  
  results_list <- list()
  
  for (r in ratios) {
    
    TP_max <- nrow(sample_df)
    TN_max <- nrow(focal_df)
    
    # ---- Compute largest feasible total_n that preserves the ratio ----
    total_max <- floor(min(TP_max / r, TN_max / (1 - r)))
    
    # If impossible, skip or warn
    if (total_max <= 0) {
      warning("Not enough data to sample ratio = ", r)
      next
    }
    
    sample_n <- round(total_max * r)      # TP
    focal_n  <- total_max - sample_n      # TN
    
    # ---- Sampling ----
    set.seed(333)
    focal_sample  <- focal_df[sample(nrow(focal_df),  focal_n,  replace = FALSE), ]
    sample_sample <- sample_df[sample(nrow(sample_df), sample_n, replace = FALSE), ]
    
    # ---- Sanity checks ----
    focal_unique  <- length(unique(focal_sample$query_kmer))
    sample_unique <- length(unique(sample_sample$query_kmer))
    
    message("\nSanity check for ratio = ", r)
    message("  - total_max rows sampled: ", total_max)
    message("  - focal_n:  ", focal_n,  "  (TN)")
    message("  - sample_n: ", sample_n, "  (TP)")
    message("  - Unique kmers in focal_sample:  ", focal_unique)
    message("  - Unique kmers in sample_sample: ", sample_unique)
    
    # ---- Combine and shuffle ----
    combined <- rbind(focal_sample, sample_sample)
    combined <- combined[sample(nrow(combined)), ]
    
    combined$ratio <- r
    
    results_list[[as.character(r)]] <- combined
  }
  
  final_df <- do.call(rbind, results_list)
  rownames(final_df) <- NULL
  
  return(final_df)
}


################################
####### Inosine_m6A ############
################################

# Task 1: Data import

UNM_inosine_m6A_DRACH_rep1 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_inosine_m6A", "pod5---bc_1_s.txt"), header = T, sep = "\t")

UNM_inosine_m6A_DRACH_rep2 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_inosine_m6A", "pod5---bc_2_s.txt"), header = T, sep = "\t")

UNM_inosine_m6A_DRACH_rep3 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_inosine_m6A", "pod5---bc_3_s.txt"), header = T, sep = "\t")

UNM_inosine_m6A_DRACH <- rbind(UNM_inosine_m6A_DRACH_rep1, UNM_inosine_m6A_DRACH_rep2, UNM_inosine_m6A_DRACH_rep3)

rm(UNM_inosine_m6A_DRACH_rep1, UNM_inosine_m6A_DRACH_rep2, UNM_inosine_m6A_DRACH_rep3)


inosine_m6A_100_rep1 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_inosine_m6A", "pod5---bc_31_s.txt"), header = T, sep = "\t")

inosine_m6A_100_rep2 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_inosine_m6A", "pod5---bc_32_s.txt"), header = T, sep = "\t")

inosine_m6A_100_rep3 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_inosine_m6A", "pod5---bc_33_s.txt"), header = T, sep = "\t")

inosine_m6A_100 <- rbind(inosine_m6A_100_rep1,inosine_m6A_100_rep2,inosine_m6A_100_rep3)

rm(inosine_m6A_100_rep1,inosine_m6A_100_rep2,inosine_m6A_100_rep3)


# Task 1: Produce functions for annotation / filtering

filter_and_annotate <- function(df, mod_code, is_true){
  
  df <- df %>%
    # Remove potential "-" strand predictions
    filter(ref_mod_strand == "+") %>%
    # Filter out any incomplete 5mers (e.g: "GCAT-")
    filter(!grepl("-", query_kmer)) %>%
    # Filter out any rows in which reference is not matching read
    filter(query_kmer == ref_kmer) %>%
    # Filter out baseQ >= 10
    filter(base_qual >= 10) %>% 
    # Filter for desired mod code
    filter(mod_code == !!mod_code) %>% 
    # Add ground-truth column
    mutate(is_true = is_true) %>% 
    # Annotate DRACH 
    mutate(is_DRACH = case_when(
      str_detect(query_kmer, "^[AGT][GA][A][C][ACT]$") ~ 1,
      T ~ 0
    )) %>%
    # Annotate single A-in 5mer DRACH 
    mutate(is_DRACH_single_A = case_when(
      str_detect(query_kmer, "^[GT][G][A][C][CT]$") ~ 1,
      T ~ 0
    )) %>%
    # Annotate single-A 5mers
    mutate(is_single_A_5mer = case_when(
      str_detect(query_kmer, "^[GTC][GTC][A][GTC][GCT]$") ~ 1,
      T ~ 0
    ))
  
}


UNM_inosine_m6A_DRACH <- filter_and_annotate(df = UNM_inosine_m6A_DRACH, mod_code = "a", is_true = 0)

inosine_m6A_100 <- filter_and_annotate(df = inosine_m6A_100, mod_code = "a", is_true = 1)

inosine_m6A_all_kmers <- sample_by_ratios(UNM_inosine_m6A_DRACH, inosine_m6A_100, ratios = c(0.5, 0.1, 0.01, 0.001))

write.table(inosine_m6A_all_kmers, file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","inosine_m6A_all_kmers.tsv"), col.names = T, row.names = F ,sep = "\t")


##############################
####### m6A_DRACH ############
##############################

# Task 1: Data import

UNM_rep1_m6A_DRACH <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m6A_DRACH", "pod5---bc_1_s.txt"), header = T, sep = "\t")

UNM_rep2_m6A_DRACH <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m6A_DRACH", "pod5---bc_2_s.txt"), header = T, sep = "\t")

UNM_rep3_m6A_DRACH <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m6A_DRACH", "pod5---bc_3_s.txt"), header = T, sep = "\t")

UNM_m6A_DRACH <- rbind(UNM_rep1_m6A_DRACH, UNM_rep2_m6A_DRACH, UNM_rep3_m6A_DRACH)

rm(UNM_rep1_m6A_DRACH, UNM_rep2_m6A_DRACH, UNM_rep3_m6A_DRACH)


m6A_DRACH_100_rep1 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m6A_DRACH", "pod5---bc_31_s.txt"), header = T, sep = "\t")

m6A_DRACH_100_rep2 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m6A_DRACH", "pod5---bc_32_s.txt"), header = T, sep = "\t")

m6A_DRACH_100_rep3 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m6A_DRACH", "pod5---bc_33_s.txt"), header = T, sep = "\t")

m6A_DRACH_100 <- rbind(m6A_DRACH_100_rep1,m6A_DRACH_100_rep2,m6A_DRACH_100_rep3)

rm(m6A_DRACH_100_rep1,m6A_DRACH_100_rep2,m6A_DRACH_100_rep3)

# Task 1: Produce functions for annotation / filtering 

UNM_m6A_DRACH <- filter_and_annotate(df = UNM_m6A_DRACH, mod_code = "a", is_true = 0)

m6A_DRACH_100 <- filter_and_annotate(df = m6A_DRACH_100, mod_code = "a", is_true = 1)

m6A_DRACH_all_kmers <- sample_by_ratios(UNM_m6A_DRACH, m6A_DRACH_100, ratios = c(0.5, 0.1, 0.01, 0.001))

write.table(m6A_DRACH_all_kmers, file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","m6A_DRACH_all_kmers.tsv"), col.names = T, row.names = F ,sep = "\t")


#########################
####### pseU ############
#########################

# Task 1: Data import

UNM_rep1_pseU <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_pseU", "pod5---bc_1_s.txt"), header = T, sep = "\t")

UNM_rep2_pseU <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_pseU", "pod5---bc_2_s.txt"), header = T, sep = "\t")

UNM_rep3_pseU <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_pseU", "pod5---bc_3_s.txt"), header = T, sep = "\t")

UNM_pseU <- rbind(UNM_rep1_pseU,UNM_rep2_pseU,UNM_rep3_pseU)

rm(UNM_rep1_pseU,UNM_rep2_pseU,UNM_rep3_pseU)


Y_pseU_rep1 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_pseU", "pod5---bc_46_s.txt"), header = T, sep = "\t")

Y_pseU_rep2 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_pseU", "pod5---bc_47_s.txt"), header = T, sep = "\t")

Y_pseU_rep3 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_pseU", "pod5---bc_48_s.txt"), header = T, sep = "\t")

Y_pseU <- rbind(Y_pseU_rep1, Y_pseU_rep2, Y_pseU_rep3)

rm(Y_pseU_rep1,Y_pseU_rep2,Y_pseU_rep3)


m1Y_pseU_rep1 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_pseU", "pod5---bc_10_s.txt"), header = T, sep = "\t")

m1Y_pseU_rep2 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_pseU", "pod5---bc_11_s.txt"), header = T, sep = "\t")

m1Y_pseU_rep3 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_pseU", "pod5---bc_12_s.txt"), header = T, sep = "\t")

m1Y_pseU <- rbind(m1Y_pseU_rep1,m1Y_pseU_rep2,m1Y_pseU_rep3)

rm(m1Y_pseU_rep1,m1Y_pseU_rep2,m1Y_pseU_rep3)


m5U_pseU_rep1 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool3_rna004_130bps_hac@v5.1.0_pseU", "pod5---bc_49_s.txt"), header = T, sep = "\t")

m5U_pseU_rep2 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool3_rna004_130bps_hac@v5.1.0_pseU", "pod5---bc_50_s.txt"), header = T, sep = "\t")

m5U_pseU_rep3 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool3_rna004_130bps_hac@v5.1.0_pseU", "pod5---bc_51_s.txt"), header = T, sep = "\t")

m5U_pseU <- rbind(m5U_pseU_rep1, m5U_pseU_rep2, m5U_pseU_rep3)

rm(m5U_pseU_rep1, m5U_pseU_rep2, m5U_pseU_rep3)
  

filter_and_annotate_pseU <- function(df, mod_code, is_true){
  
  df <- df %>%
    # Remove potential "-" strand predictions
    filter(ref_mod_strand == "+") %>%
    # Filter out any incomplete 5mers (e.g: "GCAT-")
    filter(!grepl("-", query_kmer)) %>%
    # Filter out any rows in which reference is not matching read
    filter(query_kmer == ref_kmer) %>%
    # Filter out baseQ >= 10
    filter(base_qual >= 10) %>%
    # Filter for desired mod code
    filter(mod_code == !!mod_code) %>% 
    # Add ground-truth column
    mutate(is_true = is_true) %>% 
    # Annotate single_T_5mer 
    mutate(is_single_T_5mer = case_when(
      str_detect(query_kmer, "^[AGC][AGC][T][AGC][AGC]$") ~ 1,
      T ~ 0
    ))
  
}

UNM_pseU <- filter_and_annotate_pseU(df = UNM_pseU, mod_code = 17802, is_true = 0)

Y_pseU <- filter_and_annotate_pseU(df = Y_pseU, mod_code = 17802, is_true = 1)

m1Y_pseU <- filter_and_annotate_pseU(df = m1Y_pseU, mod_code = 17802, is_true = 1)

m5U_pseU <- filter_and_annotate_pseU(df = m5U_pseU, mod_code = 17802, is_true = 1)

Y_pseU_all_kmers <- sample_by_ratios(UNM_pseU, Y_pseU, ratios = c(0.5, 0.1, 0.01, 0.001))

m1Y_pseU_all_kmers <- sample_by_ratios(UNM_pseU, m1Y_pseU, ratios = c(0.5, 0.1, 0.01, 0.001))

m5U_pseU <- sample_by_ratios(UNM_pseU, m5U_pseU, ratios = c(0.5, 0.1, 0.01, 0.001))


write.table(Y_pseU_all_kmers, file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","pseU_all_kmers.tsv"), col.names = T, row.names = F ,sep = "\t")

write.table(m1Y_pseU_all_kmers, file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","m1Y_all_kmers.tsv"), col.names = T, row.names = F ,sep = "\t")

write.table(m5U_pseU, file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","m5U_all_kmers.tsv"), col.names = T, row.names = F ,sep = "\t")


########################
####### m5C ############
########################

UNM_rep1_m5C <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m5C", "pod5---bc_1_s.txt"), header = T, sep = "\t")

UNM_rep2_m5C <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m5C", "pod5---bc_2_s.txt"), header = T, sep = "\t")

UNM_rep3_m5C <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m5C", "pod5---bc_3_s.txt"), header = T, sep = "\t")

UNM_m5C <- rbind(UNM_rep1_m5C,UNM_rep2_m5C,UNM_rep3_m5C)

rm(UNM_rep1_m5C,UNM_rep2_m5C,UNM_rep3_m5C) 


m5C_m5C_model_rep1 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m5C", "pod5---bc_13_s.txt"), header = T, sep = "\t")

m5C_m5C_model_rep2 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m5C", "pod5---bc_13_s.txt"), header = T, sep = "\t")

m5C_m5C_model_rep3 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m5C", "pod5---bc_13_s.txt"), header = T, sep = "\t")

m5C_m5C_model <- rbind(m5C_m5C_model_rep1, m5C_m5C_model_rep2, m5C_m5C_model_rep3)

rm(m5C_m5C_model_rep1, m5C_m5C_model_rep2, m5C_m5C_model_rep3)


hm5C_m5C_model_rep1 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m5C", "pod5---bc_7_s.txt"), header = T, sep = "\t")

hm5C_m5C_model_rep2 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m5C", "pod5---bc_8_s.txt"), header = T, sep = "\t")

hm5C_m5C_model_rep3 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m5C", "pod5---bc_9_s.txt"), header = T, sep = "\t")

hm5C_m5C_model <- rbind(hm5C_m5C_model_rep1,hm5C_m5C_model_rep2,hm5C_m5C_model_rep3)

rm(hm5C_m5C_model_rep1,hm5C_m5C_model_rep2,hm5C_m5C_model_rep3)


ac4C_m5C_model_rep1 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m5C", "pod5---bc_4_s.txt"), header = T, sep = "\t")

ac4C_m5C_model_rep2 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m5C", "pod5---bc_5_s.txt"), header = T, sep = "\t")

ac4C_m5C_model_rep3 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full_ref_kmer", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m5C", "pod5---bc_6_s.txt"), header = T, sep = "\t")

ac4C_m5C_model <- rbind(ac4C_m5C_model_rep1,ac4C_m5C_model_rep2,ac4C_m5C_model_rep3)

rm(ac4C_m5C_model_rep1,ac4C_m5C_model_rep2,ac4C_m5C_model_rep3)


filter_and_annotate_m5C <- function(df, mod_code, is_true){
  
  df <- df %>%
    # Remove potential "-" strand predictions
    filter(ref_mod_strand == "+") %>%
    # Filter out any incomplete 5mers (e.g: "GCAT-")
    filter(!grepl("-", query_kmer)) %>%
    # Filter out any rows in which reference is not matching read
    filter(query_kmer == ref_kmer) %>%
    # Filter out baseQ >= 10
    filter(base_qual >= 10) %>%
    # Filter for desired mod code
    filter(mod_code == !!mod_code) %>% 
    # Add ground-truth column
    mutate(is_true = is_true) %>% 
    # Annotate single_C_5mer 
    mutate(is_single_C_5mer = case_when(
      str_detect(query_kmer, "^[AGT][AGT][C][AGT][AGT]$") ~ 1,
      T ~ 0
    ))
  
}

UNM_m5C <- filter_and_annotate_m5C(df = UNM_m5C, mod_code = "m", is_true = 0)

m5C_m5C_model <- filter_and_annotate_m5C(df = m5C_m5C_model, mod_code = "m", is_true = 1)

hm5C_m5C_model <- filter_and_annotate_m5C(df = hm5C_m5C_model, mod_code = "m", is_true = 1)

ac4C_m5C_model <- filter_and_annotate_m5C(df = ac4C_m5C_model, mod_code = "m", is_true = 1)

m5C_m5C_model_all_kmers <- sample_by_ratios(UNM_m5C, m5C_m5C_model, ratios = c(0.5, 0.1, 0.01, 0.001))

hm5C_m5C_model_all_kmers <- sample_by_ratios(UNM_m5C, hm5C_m5C_model, ratios = c(0.5, 0.1, 0.01, 0.001))

ac4C_m5C_model_all_kmers  <- sample_by_ratios(UNM_m5C, ac4C_m5C_model, ratios = c(0.5, 0.1, 0.01, 0.001))


write.table(m5C_m5C_model_all_kmers, file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","m5C_all_kmers.tsv"), col.names = T, row.names = F ,sep = "\t")

write.table(hm5C_m5C_model_all_kmers, file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","hm5C_model_all_kmers.tsv"), col.names = T, row.names = F ,sep = "\t")

write.table(ac4C_m5C_model_all_kmers, file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","ac4C_model_all_kmers.tsv"), col.names = T, row.names = F ,sep = "\t")

