#############################################################
## Generate master table with Threshold-Dependant Metrics  ##
#############################################################

# Task 0: Load packages

pkgs <- c("openintro","tidyverse","backports","here","skimr","dplyr", "ggplot2", "ggsci","ggforce",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr","scales")

lapply(pkgs, library, character.only = TRUE)



#####################################
# Task 1: Import and annotate Data  #
#####################################

inosine_m6A <- read.table(file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","inosine_m6A_all_kmers.tsv"), header = T, sep = "\t") %>% 
  mutate(
    model = "inosine_m6A",
    modification = "m6A",
    motif = "NNANN"
  ) %>% 
  relocate(model, modification, motif, .before = 1)

m6A_DRACH <- read.table(file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","m6A_DRACH_all_kmers.tsv"), header = T, sep = "\t")  %>% 
  mutate(
    model = "m6A_DRACH",
    modification = "m6A",
    motif = "NNANN"
  ) %>% 
  relocate(model, modification, motif, .before = 1)



Y_pseU_all_kmers <- read.table(file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","pseU_all_kmers.tsv"), header = T, sep = "\t") %>% 
  mutate(
    model = "pseU",
    modification = "Y",
    motif = "NNTNN"
  ) %>% 
  relocate(model, modification, motif, .before = 1)

m1Y_pseU_all_kmers <- read.table(file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","m1Y_all_kmers.tsv"), header = T, sep = "\t") %>% 
  mutate(
    model = "pseU",
    modification = "m1Y",
    motif = "NNTNN"
  ) %>% 
  relocate(model, modification, motif, .before = 1)

m5U_pseU_all_kmers <- read.table(file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","m5U_all_kmers.tsv"), header = T, sep = "\t") %>% 
  mutate(
    model = "pseU",
    modification = "m5U",
    motif = "NNTNN"
  ) %>% 
  relocate(model, modification, motif, .before = 1)



m5C_m5C_all_kmers <- read.table(file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","m5C_all_kmers.tsv"), header = T, sep = "\t") %>% 
  mutate(
    model = "m5C",
    modification = "m5C",
    motif = "NNCNN"
  ) %>% 
  relocate(model, modification, motif, .before = 1)

hm5C_m5C_all_kmers <- read.table(file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","hm5C_model_all_kmers.tsv"), header = T, sep = "\t") %>% 
  mutate(
    model = "m5C",
    modification = "hm5C",
    motif = "NNCNN"
  ) %>% 
  relocate(model, modification, motif, .before = 1)

ac4C_m5C_all_kmers <- read.table(file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","ac4C_model_all_kmers.tsv"), header = T, sep = "\t") %>% 
  mutate(
    model = "m5C",
    modification = "ac4C",
    motif = "NNCNN"
  ) %>% 
  relocate(model, modification, motif, .before = 1)


###########################################
## Calculate Threshold-Dependent Metrics ##
###########################################


### Calculate Threshold-Dependent Metrics

compute_metrics_fast <- function(
    df,
    score_col  = "mod_qual",   # continuous model score in [0,1]
    true_col   = "is_true",    # 0/1 ground truth
    thresholds = c(seq(0, 0.9, 0.1), 0.99)
) {
  score_sym <- sym(score_col)
  true_sym  <- sym(true_col)
  
  # group by model metadata + ratio, but rename ratio in output
  df %>%
    group_by(
      model,
      modification,
      motif,
      mod_to_unmod_ratio = ratio
    ) %>%
    group_modify(~ {
      g <- .x
      scores <- g[[score_col]]
      y      <- g[[true_col]]
      
      # logical vectors for positives/negatives
      y_pos <- (y == 1)
      y_neg <- !y_pos
      
      # matrix: rows = samples, cols = thresholds
      # pred_mat[i, j] = (score_i >= thresholds[j])
      pred_mat <- outer(scores, thresholds, `>=`)
      
      # confusion matrix per threshold (vectorised over columns)
      TP <- colSums(pred_mat  & y_pos)
      FP <- colSums(pred_mat  & y_neg)
      TN <- colSums(!pred_mat & y_neg)
      FN <- colSums(!pred_mat & y_pos)
      
      total <- TP + FP + TN + FN
      
      # percentages
      TP_pct <- TP / total
      FP_pct <- FP / total
      TN_pct <- TN / total
      FN_pct <- FN / total
      
      # metrics
      precision <- ifelse(TP + FP == 0, 0, TP / (TP + FP))
      recall    <- ifelse(TP + FN == 0, 0, TP / (TP + FN))        # TPR
      f1        <- ifelse(precision + recall == 0, 0,
                          2 * precision * recall / (precision + recall))
      
      TNR <- ifelse(TN + FP == 0, 0, TN / (TN + FP))              # specificity
      FPR <- ifelse(FP + TN == 0, 0, FP / (FP + TN))              # false positive rate
      
      tibble(
        mod_prob_cutoff = thresholds,
        TP  = TP,
        FP  = FP,
        TN  = TN,
        FN  = FN,
        TP_pct = TP_pct,
        FP_pct = FP_pct,
        TN_pct = TN_pct,
        FN_pct = FN_pct,
        precision = precision,
        recall    = recall,
        f1        = f1,
        TNR       = TNR,
        FPR       = FPR
      )
    }) %>%
    ungroup()
}



inosine_m6A_summary <- compute_metrics_fast(inosine_m6A)


m6A_DRACH_summary <- compute_metrics_fast(m6A_DRACH)


Y_pseU_summary <- compute_metrics_fast(Y_pseU_all_kmers)


m1Y_pseU_summary <- compute_metrics_fast(m1Y_pseU_all_kmers)


m5U_pseU_summary <- compute_metrics_fast(m5U_pseU_all_kmers)


m5C_m5C_summary <- compute_metrics_fast(m5C_m5C_all_kmers)


hm5C_m5C_summary <- compute_metrics_fast(hm5C_m5C_all_kmers)


ac4C_m5C_summary <- compute_metrics_fast(ac4C_m5C_all_kmers)


master_table <- rbind(inosine_m6A_summary, m6A_DRACH_summary, Y_pseU_summary, m1Y_pseU_summary, m5U_pseU_summary, m5C_m5C_summary, hm5C_m5C_summary, ac4C_m5C_summary)


write.table(master_table, file = here("..","__results_Gregor","curlcakes","per_read","hac","tables","model_metrics_per_ratio.tsv"), col.names = T, row.names = F ,sep = "\t")








