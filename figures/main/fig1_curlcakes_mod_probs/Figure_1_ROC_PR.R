###################################################################
## Plotting ROC/PR curves on sampled modification probabilities  ##
###################################################################

# Task 0: Load packages

pkgs <- c("openintro","tidyverse","backports","here","skimr","dplyr", "ggplot2", "ggsci","ggforce",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr","scales")

lapply(pkgs, library, character.only = TRUE)


##############################
# m6A_DRACH and inosine_m6A ##
##############################

#####################################
# Task 1: Import and annotate Data  #
#####################################


UNM_rep1_m6A_DRACH <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m6A_DRACH", "pod5---bc_1_s.txt"), header = T, sep = "\t")

UNM_rep2_m6A_DRACH <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m6A_DRACH", "pod5---bc_2_s.txt"), header = T, sep = "\t")

UNM_rep3_m6A_DRACH <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m6A_DRACH", "pod5---bc_3_s.txt"), header = T, sep = "\t")

UNM_m6A_DRACH <- rbind(UNM_rep1_m6A_DRACH,UNM_rep2_m6A_DRACH,UNM_rep3_m6A_DRACH)

rm(UNM_rep1_m6A_DRACH,UNM_rep2_m6A_DRACH,UNM_rep3_m6A_DRACH)


UNM_rep1_inosine_m6A_DRACH <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_inosine_m6A", "pod5---bc_1_s.txt"), header = T, sep = "\t")

UNM_rep2_inosine_m6A_DRACH <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_inosine_m6A", "pod5---bc_2_s.txt"), header = T, sep = "\t")

UNM_rep3_inosine_m6A_DRACH <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_inosine_m6A", "pod5---bc_3_s.txt"), header = T, sep = "\t")

UNM_inosine_m6A_DRACH <- rbind(UNM_rep1_inosine_m6A_DRACH,UNM_rep2_inosine_m6A_DRACH,UNM_rep3_inosine_m6A_DRACH)

rm(UNM_rep1_inosine_m6A_DRACH,UNM_rep2_inosine_m6A_DRACH,UNM_rep3_inosine_m6A_DRACH)


m6A_DRACH_100 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m6A_DRACH", "pod5---bc_31_s.txt"), header = T, sep = "\t") 


inosine_m6A_100 <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_inosine_m6A", "pod5---bc_31_s.txt"), header = T, sep = "\t")



filter_and_annotate <- function(df, mod_code, is_true){
  
  df <- df %>%
    # Remove potential "-" strand predictions
    filter(ref_mod_strand == "+") %>%
    # Filter out any incomplete 5mers (e.g: "GCAT-")
    filter(!grepl("-", query_kmer)) %>%
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

UNM_m6A_DRACH <- filter_and_annotate(df = UNM_m6A_DRACH, mod_code = "a", is_true = 0)

UNM_inosine_m6A_DRACH <- filter_and_annotate(df = UNM_inosine_m6A_DRACH, mod_code = "a", is_true = 0)

m6A_DRACH_100 <- filter_and_annotate(df = m6A_DRACH_100, mod_code = "a", is_true = 1)

inosine_m6A_100 <- filter_and_annotate(df = inosine_m6A_100, mod_code = "a", is_true = 1)


##########################
# Task 2: Data Wrangling #
##########################


# Task : Create a balanced datasets

# Set seed to make sampling reads reproducible

################
### All kmers ##
################

set.seed(123)

list_of_dfs <- list(UNM_inosine_m6A_DRACH, inosine_m6A_100)

min_rows <- min(sapply(list_of_dfs, nrow))

UNM_inosine_m6A_DRACH_all_kmers <- slice_sample(UNM_inosine_m6A_DRACH, n = min_rows)

inosine_m6A_100_all_kmers <- slice_sample(inosine_m6A_100, n = min_rows)

inosine_m6A_all_kmers <- rbind(UNM_inosine_m6A_DRACH_all_kmers, inosine_m6A_100_all_kmers)

rm(UNM_inosine_m6A_DRACH_all_kmers,inosine_m6A_100_all_kmers)

###############################
### All kmers central A only ##
###############################

UNM_inosine_m6A_all_kmers_singleA <- UNM_inosine_m6A_DRACH %>%
  filter(is_single_A_5mer == 1)

inosine_m6A_all_kmers_singleA <- inosine_m6A_100 %>% 
  filter(is_single_A_5mer == 1)

list_of_dfs <- list(UNM_inosine_m6A_all_kmers_singleA, inosine_m6A_all_kmers_singleA)

min_rows <- min(sapply(list_of_dfs, nrow))

UNM_inosine_m6A_all_kmers_singleA <- slice_sample(UNM_inosine_m6A_all_kmers_singleA, n = min_rows)

inosine_m6A_all_kmers_singleA <- slice_sample(inosine_m6A_all_kmers_singleA, n = min_rows)

inosine_m6A_all_kmers_singleA_final <- rbind(UNM_inosine_m6A_all_kmers_singleA, inosine_m6A_all_kmers_singleA)

rm(UNM_inosine_m6A_all_kmers_singleA,inosine_m6A_all_kmers_singleA)


############
### DRACH ##
############

UNM_DRACH_100_kmer_DRACH <- UNM_m6A_DRACH %>%
  filter(is_DRACH == 1)

m6A_DRACH_100_kmer_DRACH <- m6A_DRACH_100 %>% 
  filter(is_DRACH == 1)

UNM_inosine_m6A_kmer_DRACH <- UNM_inosine_m6A_DRACH %>%
  filter(is_DRACH == 1)

inosine_m6A_kmer_DRACH <- inosine_m6A_100 %>% 
  filter(is_DRACH == 1)

list_of_dfs <- list(UNM_DRACH_100_kmer_DRACH, m6A_DRACH_100_kmer_DRACH, UNM_inosine_m6A_kmer_DRACH, inosine_m6A_kmer_DRACH)

min_rows <- min(sapply(list_of_dfs, nrow))


UNM_DRACH_100_kmer_DRACH <- slice_sample(UNM_DRACH_100_kmer_DRACH, n = min_rows)

m6A_DRACH_100_kmer_DRACH <- slice_sample(m6A_DRACH_100_kmer_DRACH, n = min_rows)

UNM_inosine_m6A_kmer_DRACH <- slice_sample(UNM_inosine_m6A_kmer_DRACH, n = min_rows)

inosine_m6A_kmer_DRACH <- slice_sample(inosine_m6A_kmer_DRACH, n = min_rows)


m6A_DRACH_kmer_DRACH <- rbind(UNM_DRACH_100_kmer_DRACH, m6A_DRACH_100_kmer_DRACH)

inosine_m6A_kmer_DRACH_final <- rbind(UNM_inosine_m6A_kmer_DRACH, inosine_m6A_kmer_DRACH)

rm(UNM_DRACH_100_kmer_DRACH,m6A_DRACH_100_kmer_DRACH,UNM_inosine_m6A_kmer_DRACH,inosine_m6A_kmer_DRACH)


############################
### DRACH central A only ##
###########################

UNM_DRACH_100_singleA_DRACH <- UNM_m6A_DRACH %>%
  filter(is_DRACH_single_A == 1)

m6A_DRACH_100_singleA_DRACH <- m6A_DRACH_100 %>% 
  filter(is_DRACH_single_A == 1)

UNM_inosine_m6A_singleA_DRACH <- UNM_inosine_m6A_DRACH %>%
  filter(is_DRACH_single_A == 1)

inosine_m6A_singleA_DRACH <- inosine_m6A_100 %>% 
  filter(is_DRACH_single_A == 1)

list_of_dfs <- list(UNM_DRACH_100_singleA_DRACH, m6A_DRACH_100_singleA_DRACH, UNM_inosine_m6A_singleA_DRACH, inosine_m6A_singleA_DRACH)

min_rows <- min(sapply(list_of_dfs, nrow))


UNM_DRACH_100_singleA_DRACH <- slice_sample(UNM_DRACH_100_singleA_DRACH, n = min_rows)

m6A_DRACH_100_singleA_DRACH <- slice_sample(m6A_DRACH_100_singleA_DRACH, n = min_rows)

UNM_inosine_m6A_singleA_DRACH <- slice_sample(UNM_inosine_m6A_singleA_DRACH, n = min_rows)

inosine_m6A_singleA_DRACH <- slice_sample(inosine_m6A_singleA_DRACH, n = min_rows)


m6A_DRACH_singleA_DRACH_final <- rbind(UNM_DRACH_100_singleA_DRACH, m6A_DRACH_100_singleA_DRACH)

inosine_m6A_singleA_DRACH_final <- rbind(UNM_inosine_m6A_singleA_DRACH, inosine_m6A_singleA_DRACH)

rm(UNM_DRACH_100_singleA_DRACH,m6A_DRACH_100_singleA_DRACH,UNM_inosine_m6A_singleA_DRACH,inosine_m6A_singleA_DRACH)




####################
# Task 3: Plot ROC #
####################


library(pROC)


### All context ###

roc_inosine_m6A_all_kmers <- roc(response = inosine_m6A_all_kmers$is_true, predictor = inosine_m6A_all_kmers$mod_qual)

### All context single A ###

roc_inosine_m6A_all_kmers_singleA_final <- roc(response = inosine_m6A_all_kmers_singleA_final$is_true, predictor = inosine_m6A_all_kmers_singleA_final$mod_qual)

### DRACH context ### 

roc_m6A_DRACH_kmer_DRACH <- roc(response = m6A_DRACH_kmer_DRACH$is_true, predictor = m6A_DRACH_kmer_DRACH$mod_qual)

roc_inosine_m6A_kmer_DRACH_final <- roc(response = inosine_m6A_kmer_DRACH_final$is_true, predictor = inosine_m6A_kmer_DRACH_final$mod_qual)

### DRACH single A context ###

roc_m6A_DRACH_singleA_DRACH_final <- roc(response = m6A_DRACH_singleA_DRACH_final$is_true, predictor = m6A_DRACH_singleA_DRACH_final$mod_qual)

roc_inosine_m6A_singleA_DRACH_final <- roc(response = inosine_m6A_singleA_DRACH_final$is_true, predictor = inosine_m6A_singleA_DRACH_final$mod_qual)



########################################
## m6A_DRACH and inosine_m6A combined ##
########################################

# Open a PDF device
pdf(file = here("..","__results_Gregor","curlcakes","per_read","hac","ROC","m6A_DRACH_inosine_m6A_combined.pdf"), 
    width = 8, height = 8)

# Plot the first ROC curve to set up the plot
plot(roc_inosine_m6A_all_kmers, col = "#299189", main = "ROC Curve for Percentage Modified", lwd = 4)

# Add the remaining ROC curves using lines()
lines(roc_inosine_m6A_all_kmers_singleA_final, col = "#67BEB4", lwd = 4)
lines(roc_inosine_m6A_kmer_DRACH_final, col = "#8CD5CA", lwd = 4)
lines(roc_inosine_m6A_singleA_DRACH_final, col = "#B3ECE0", lwd = 4)
lines(roc_m6A_DRACH_kmer_DRACH, col = "#381a61", lwd = 4)
lines(roc_m6A_DRACH_singleA_DRACH_final, col = "#E6D4FF", lwd = 4)


# Compute AUC values
auc_values <- c(
  auc(roc_inosine_m6A_all_kmers),
  auc(roc_inosine_m6A_all_kmers_singleA_final),
  auc(roc_inosine_m6A_kmer_DRACH_final),
  auc(roc_inosine_m6A_singleA_DRACH_final),
  auc(roc_m6A_DRACH_kmer_DRACH),
  auc(roc_m6A_DRACH_singleA_DRACH_final)
)

# Add a legend to identify curves and AUC values
legend("bottomright", legend = c(
  paste("inosine_m6A | NNANN = ", round(auc_values[1], 3)),
  paste("inosine_m6A | BBABB = ", round(auc_values[2], 3)),
  paste("inosine_m6A | DRACH = ", round(auc_values[3], 3)),
  paste("inosine_m6A | KGACY = ", round(auc_values[4], 3)),
  paste("m6A_DRACH | DRACH = ", round(auc_values[5], 3)),
  paste("m6A_DRACH | KGACY = ", round(auc_values[6], 3))
), col = c("#299189", "#67BEB4", "#8CD5CA", "#B3ECE0", "#381a61","#E6D4FF"), lwd = 4, cex = 1, title = "Area under the curve [AUC]:")


# Close the PDF device
dev.off()


################
## m6A_DRACH  ##
################

# Open a PDF device
pdf(file = here("..","__results_Gregor","curlcakes","per_read","hac","ROC","m6A_DRACH.pdf"), 
    width = 8, height = 8)

# Plot the first ROC curve to set up the plot
plot(roc_m6A_DRACH_kmer_DRACH, col = "#381a61", lwd = 4, main = "ROC Curve for Percentage Modified")

# Add the remaining ROC curves using lines()
lines(roc_m6A_DRACH_singleA_DRACH_final, col = "#E6D4FF", lwd = 4)


# Compute AUC values
auc_values <- c(
  auc(roc_m6A_DRACH_kmer_DRACH),
  auc(roc_m6A_DRACH_singleA_DRACH_final)
)

# Add a legend to identify curves and AUC values
legend("bottomright", legend = c(
  paste("m6A_DRACH | DRACH = ", round(auc_values[1], 3)),
  paste("m6A_DRACH | KGACY = ", round(auc_values[2], 3))
), col = c("#381a61","#E6D4FF"), lwd = 4, cex = 1, title = "Area under the curve [AUC]:")

# Close the PDF device
dev.off()


################
## inosine_m6A  ##
################

# Open a PDF device
pdf(file = here("..","__results_Gregor","curlcakes","per_read","hac","ROC","inosine_m6A.pdf"), 
    width = 8, height = 8)

# Plot the first ROC curve to set up the plot
plot(roc_inosine_m6A_all_kmers, col = "#299189", main = "ROC Curve for Percentage Modified", lwd = 4)

# Add the remaining ROC curves using lines()
lines(roc_inosine_m6A_all_kmers_singleA_final, col = "#67BEB4", lwd = 4)
lines(roc_inosine_m6A_kmer_DRACH_final, col = "#8CD5CA", lwd = 4)
lines(roc_inosine_m6A_singleA_DRACH_final, col = "#B3ECE0", lwd = 4)


# Compute AUC values
auc_values <- c(
  auc(roc_inosine_m6A_all_kmers),
  auc(roc_inosine_m6A_all_kmers_singleA_final),
  auc(roc_inosine_m6A_kmer_DRACH_final),
  auc(roc_inosine_m6A_singleA_DRACH_final)
)

# Add a legend to identify curves and AUC values
legend("bottomright", legend = c(
  paste("inosine_m6A | NNANN = ", round(auc_values[1], 3)),
  paste("inosine_m6A | BBABB = ", round(auc_values[2], 3)),
  paste("inosine_m6A | DRACH = ", round(auc_values[3], 3)),
  paste("inosine_m6A | KGACY = ", round(auc_values[4], 3))
), col = c("#299189", "#67BEB4", "#8CD5CA", "#B3ECE0"), lwd = 4, cex = 1, title = "Area under the curve [AUC]:")


# Close the PDF device
dev.off()


###################
# Task 4: Plot PR #
###################

# Load necessary libraries
library(PRROC)
library(ggplot2)



# Compute Precision-Recall curves for each dataset
pr_inosine_m6A_all_kmers <- pr.curve(scores.class0 = inosine_m6A_all_kmers$mod_qual, 
                                     weights.class0 = inosine_m6A_all_kmers$is_true, curve = TRUE)

pr_inosine_m6A_all_kmers_singleA_final <- pr.curve(scores.class0 = inosine_m6A_all_kmers_singleA_final$mod_qual, 
                                                   weights.class0 = inosine_m6A_all_kmers_singleA_final$is_true, curve = TRUE)

pr_m6A_DRACH_kmer_DRACH <- pr.curve(scores.class0 = m6A_DRACH_kmer_DRACH$mod_qual, 
                                    weights.class0 = m6A_DRACH_kmer_DRACH$is_true, curve = TRUE)

pr_inosine_m6A_kmer_DRACH_final <- pr.curve(scores.class0 = inosine_m6A_kmer_DRACH_final$mod_qual, 
                                            weights.class0 = inosine_m6A_kmer_DRACH_final$is_true, curve = TRUE)

pr_m6A_DRACH_singleA_DRACH_final <- pr.curve(scores.class0 = m6A_DRACH_singleA_DRACH_final$mod_qual, 
                                             weights.class0 = m6A_DRACH_singleA_DRACH_final$is_true, curve = TRUE)

pr_inosine_m6A_singleA_DRACH_final <- pr.curve(scores.class0 = inosine_m6A_singleA_DRACH_final$mod_qual, 
                                               weights.class0 = inosine_m6A_singleA_DRACH_final$is_true, curve = TRUE)


########################################
## m6A_DRACH and inosine_m6A combined ##
########################################

# Convert PR curve data into a dataframe for ggplot

pr_data <- list(
  pr_inosine_m6A_all_kmers = pr_inosine_m6A_all_kmers$curve,
  pr_inosine_m6A_all_kmers_singleA_final = pr_inosine_m6A_all_kmers_singleA_final$curve,
  pr_inosine_m6A_kmer_DRACH_final = pr_inosine_m6A_kmer_DRACH_final$curve,
  pr_inosine_m6A_singleA_DRACH_final = pr_inosine_m6A_singleA_DRACH_final$curve,
  pr_m6A_DRACH_kmer_DRACH = pr_m6A_DRACH_kmer_DRACH$curve,
  pr_m6A_DRACH_singleA_DRACH_final = pr_m6A_DRACH_singleA_DRACH_final$curve
)

# Assign colors matching your ROC plot
colors <- c("#299189", "#67BEB4", "#8CD5CA", "#B3ECE0", "#381a61", "#E6D4FF")


# Open a EPS device
postscript(file = here("..","__results_Gregor","curlcakes","per_read","hac","precision_recall","m6A_DRACH_inosine_m6A_combined.eps"),
           width = 8, height = 10)

# Create an empty plot
plot(NA, xlim = c(0, 1), ylim = c(0.2, 1), xlab = "Recall", ylab = "Precision", 
     main = "Precision-Recall Curve", type = "n")

# Add a dashed black horizontal baseline at precision = 0.5 
abline(h = 0.5, col = "black", lty = 2, lwd = 3)

# Add lines for each PR curve
i <- 1
for (curve_name in names(pr_data)) {
  lines(pr_data[[curve_name]][,1], pr_data[[curve_name]][,2], col = colors[i], lwd = 4)
  i <- i + 1
}

# Compute AUC-PR values
pr_auc_values <- sapply(pr_data, function(pr) pr.curve(scores.class0 = pr[,1], weights.class0 = pr[,2])$auc.integral)

# Add a legend
legend("bottomleft", legend = c(
  paste("inosine_m6A | NNANN = ", round(pr_auc_values[1], 3)),
  paste("inosine_m6A | BBABB = ", round(pr_auc_values[2], 3)),
  paste("inosine_m6A | DRACH = ", round(pr_auc_values[3], 3)),
  paste("inosine_m6A | KGACY = ", round(pr_auc_values[4], 3)),
  paste("m6A_DRACH | DRACH = ", round(pr_auc_values[5], 3)),
  paste("m6A_DRACH | KGACY = ", round(pr_auc_values[6], 3))
), col = colors, lwd = 4, cex = 1, title = "Area under the PR curve [AUC-PR]:")

# Close the PDF device
dev.off()

################
## m6A_DRACH  ##
################

# Convert PR curve data into a dataframe for ggplot

pr_data <- list(
  pr_m6A_DRACH_kmer_DRACH = pr_m6A_DRACH_kmer_DRACH$curve,
  pr_m6A_DRACH_singleA_DRACH_final = pr_m6A_DRACH_singleA_DRACH_final$curve
)

# Assign colors matching your ROC plot
colors <- c("#381a61", "#E6D4FF")


# Open a PDF device
pdf(file = here("..","__results_Gregor","curlcakes","per_read","hac","precision_recall","m6A_DRACH.pdf"), 
    width = 8, height = 8)


# Create an empty plot
plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "Recall", ylab = "Precision", 
     main = "Precision-Recall Curve", type = "n")

# Add lines for each PR curve
i <- 1
for (curve_name in names(pr_data)) {
  lines(pr_data[[curve_name]][,1], pr_data[[curve_name]][,2], col = colors[i], lwd = 4)
  i <- i + 1
}

# Compute AUC-PR values
pr_auc_values <- sapply(pr_data, function(pr) pr.curve(scores.class0 = pr[,1], weights.class0 = pr[,2])$auc.integral)

# Add a legend
legend("bottomleft", legend = c(
  paste("m6A_DRACH | DRACH = ", round(pr_auc_values[1], 3)),
  paste("m6A_DRACH | KGACY = ", round(pr_auc_values[2], 3))
), col = colors, lwd = 4, cex = 1, title = "Area under the PR curve [AUC-PR]:")

# Close the PDF device
dev.off()

##################
## inosine_m6A  ##
##################

# Convert PR curve data into a dataframe for ggplot

pr_data <- list(
  pr_inosine_m6A_all_kmers = pr_inosine_m6A_all_kmers$curve,
  pr_inosine_m6A_all_kmers_singleA_final = pr_inosine_m6A_all_kmers_singleA_final$curve,
  pr_inosine_m6A_kmer_DRACH_final = pr_inosine_m6A_kmer_DRACH_final$curve,
  pr_inosine_m6A_singleA_DRACH_final = pr_inosine_m6A_singleA_DRACH_final$curve
)

# Assign colors matching your ROC plot
colors <- c("#299189", "#67BEB4", "#8CD5CA", "#B3ECE0")


# Open a PDF device
pdf(file = here("..","__results_Gregor","curlcakes","per_read","hac","precision_recall","inosine_m6A.pdf"), 
    width = 8, height = 8)


# Create an empty plot
plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "Recall", ylab = "Precision", 
     main = "Precision-Recall Curve", type = "n")

# Add lines for each PR curve
i <- 1
for (curve_name in names(pr_data)) {
  lines(pr_data[[curve_name]][,1], pr_data[[curve_name]][,2], col = colors[i], lwd = 4)
  i <- i + 1
}

# Compute AUC-PR values
pr_auc_values <- sapply(pr_data, function(pr) pr.curve(scores.class0 = pr[,1], weights.class0 = pr[,2])$auc.integral)

# Add a legend
legend("bottomleft", legend = c(
  paste("inosine_m6A | NNANN = ", round(pr_auc_values[1], 3)),
  paste("inosine_m6A | BBABB = ", round(pr_auc_values[2], 3)),
  paste("inosine_m6A | DRACH = ", round(pr_auc_values[3], 3)),
  paste("inosine_m6A | KGACY = ", round(pr_auc_values[4], 3))),
  col = colors, lwd = 4, cex = 1, title = "Area under the PR curve [AUC-PR]:")

# Close the PDF device
dev.off()



#################
### pseU model ##
#################

#####################################
# Task 1: Import and annotate Data  #
#####################################


UNM_rep1_pseU <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_pseU", "pod5---bc_1_s.txt"), header = T, sep = "\t")

UNM_rep2_pseU <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_pseU", "pod5---bc_2_s.txt"), header = T, sep = "\t")

UNM_rep3_pseU <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_pseU", "pod5---bc_3_s.txt"), header = T, sep = "\t")

UNM_pseU <- rbind(UNM_rep1_pseU,UNM_rep2_pseU,UNM_rep3_pseU)

rm(UNM_rep1_pseU,UNM_rep2_pseU,UNM_rep3_pseU)


Y_pseU <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_pseU", "pod5---bc_46_s.txt"), header = T, sep = "\t")

m1Y_pseU <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_pseU", "pod5---bc_10_s.txt"), header = T, sep = "\t")

m5U_pseU <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241016_Pool3_rna004_130bps_hac@v5.1.0_pseU", "pod5---bc_49_s.txt"), header = T, sep = "\t")



filter_and_annotate_pseU <- function(df, mod_code, is_true){
  
  df <- df %>%
    # Remove potential "-" strand predictions
    filter(ref_mod_strand == "+") %>%
    # Filter out any incomplete 5mers (e.g: "GCAT-")
    filter(!grepl("-", query_kmer)) %>%
    # Filter out baseQ >= 10
    filter(base_qual >= 10) %>%
    # Filter for desired mod code
    filter(mod_code == !!mod_code) %>% 
    # Add ground-truth column
    mutate(is_true = is_true) %>% 
    # Annotate DRACH 
    mutate(is_single_T_5mer = case_when(
      str_detect(query_kmer, "^[AGC][AGC][T][AGC][AGC]$") ~ 1,
      T ~ 0
    ))
  
}

UNM_pseU <- filter_and_annotate_pseU(df = UNM_pseU, mod_code = 17802, is_true = 0)

Y_pseU <- filter_and_annotate_pseU(df = Y_pseU, mod_code = 17802, is_true = 1)

m1Y_pseU <- filter_and_annotate_pseU(df = m1Y_pseU, mod_code = 17802, is_true = 1)

m5U_pseU <- filter_and_annotate_pseU(df = m5U_pseU, mod_code = 17802, is_true = 1)


##########################
# Task 2: Data Wrangling #
##########################


# Task : Create a balanced datasets

# Set seed to make sampling reads reproducible

################
### All kmers ##
################

set.seed(123)

list_of_dfs <- list(UNM_pseU, Y_pseU, m1Y_pseU, m5U_pseU)

min_rows <- min(sapply(list_of_dfs, nrow))

UNM_pseU_final <- slice_sample(UNM_pseU, n = min_rows)

Y_pseU_final <- slice_sample(Y_pseU, n = min_rows)

m1Y_pseU_final <- slice_sample(m1Y_pseU, n = min_rows)

m5U_pseU_final <- slice_sample(m5U_pseU, n = min_rows)


Y_all_kmers <- rbind(Y_pseU_final, UNM_pseU_final)

m1Y_all_kmers <- rbind(m1Y_pseU_final, UNM_pseU_final)

m5U_all_kmers <- rbind(m5U_pseU_final, UNM_pseU_final)


rm(UNM_pseU_final, Y_pseU_final, m1Y_pseU_final, m5U_pseU_final)


###############################
### All kmers central T only ##
###############################


UNM_pseU_single_T <- UNM_pseU %>% 
  filter(is_single_T_5mer == 1)

Y_pseU_single_T <- Y_pseU %>% 
  filter(is_single_T_5mer == 1)

m1Y_pseU_single_T <- m1Y_pseU %>% 
  filter(is_single_T_5mer == 1)

m5U_pseU_single_T <- m5U_pseU %>% 
  filter(is_single_T_5mer == 1)


list_of_dfs <- list(UNM_pseU_single_T, Y_pseU_single_T, m1Y_pseU_single_T, m5U_pseU_single_T)

min_rows <- min(sapply(list_of_dfs, nrow))


UNM_pseU_single_T_final <- slice_sample(UNM_pseU_single_T, n = min_rows)

Y_pseU_single_T_final <- slice_sample(Y_pseU_single_T, n = min_rows)

m1Y_pseU_single_T_final <- slice_sample(m1Y_pseU_single_T, n = min_rows)

m5U_pseU_single_T_final <- slice_sample(m5U_pseU_single_T, n = min_rows)


Y_single_T <- rbind(Y_pseU_single_T_final, UNM_pseU_single_T_final)

m1Y_single_T <- rbind(m1Y_pseU_single_T_final, UNM_pseU_single_T_final)

m5U_single_T <- rbind(m5U_pseU_single_T_final, UNM_pseU_single_T_final)


rm(UNM_pseU_single_T_final,Y_pseU_single_T_final,m1Y_pseU_single_T_final,m5U_pseU_single_T_final)



####################
# Task 3: Plot ROC #
####################


library(pROC)


### All context ###

roc_Y_all_kmers <- roc(response = Y_all_kmers$is_true, predictor = Y_all_kmers$mod_qual)

roc_m1Y_all_kmers <- roc(response = m1Y_all_kmers$is_true, predictor = m1Y_all_kmers$mod_qual)

roc_m5U_all_kmers <- roc(response = m5U_all_kmers$is_true, predictor = m5U_all_kmers$mod_qual)

### All context single T ###

roc_Y_single_T <- roc(response = Y_single_T$is_true, predictor = Y_single_T$mod_qual)

roc_m1Y_single_T <- roc(response = m1Y_single_T$is_true, predictor = m1Y_single_T$mod_qual)

roc_m5U_single_T <- roc(response = m5U_single_T$is_true, predictor = m5U_single_T$mod_qual)


############################
## pseU all mods combined ##
############################

# Open a PDF device
pdf(file = here("..","__results_Gregor","curlcakes","per_read","hac","ROC","pseU_allmods_combined.pdf"), 
    width = 8, height = 8)

# Plot the first ROC curve to set up the plot
plot(roc_Y_all_kmers, col = "#E58112", main = "ROC Curve for Percentage Modified", lwd = 4)

# Add the remaining ROC curves using lines()
lines(roc_Y_single_T, col = "#F9CFA2", lwd = 4)
lines(roc_m1Y_all_kmers, col = "#EE9A91", lwd = 4)
lines(roc_m1Y_single_T, col = "#F9E0DB", lwd = 4)
lines(roc_m5U_all_kmers, col = "#82476F", lwd = 4)
lines(roc_m5U_single_T, col = "#D7AEC6", lwd = 4)


# Compute AUC values
auc_values <- c(
  auc(roc_Y_all_kmers),
  auc(roc_Y_single_T),
  auc(roc_m1Y_all_kmers),
  auc(roc_m1Y_single_T),
  auc(roc_m5U_all_kmers),
  auc(roc_m5U_single_T)
)

# Add a legend to identify curves and AUC values
legend("bottomright", legend = c(
  paste("Y | NNTNN = ", round(auc_values[1], 3)),
  paste("Y | VVTVV = ", round(auc_values[2], 3)),
  paste("m1Y | NNTNN = ", round(auc_values[3], 3)),
  paste("m1Y | VVTVV = ", round(auc_values[4], 3)),
  paste("m5U | NNTNN = ", round(auc_values[5], 3)),
  paste("m5U | VVTVV = ", round(auc_values[6], 3))
), col = c("#E58112", "#F9CFA2", "#EE9A91", "#F9E0DB", "#82476F","#D7AEC6"), lwd = 4, cex = 1, title = "Area under the curve [AUC]:")


# Close the PDF device
dev.off()


###############
## PR curve ###
###############

#############################
## pseU all mods combined  ##
#############################

# Compute Precision-Recall curves for each dataset
pr_Y_all_kmers <- pr.curve(scores.class0 = as.numeric(Y_all_kmers$mod_qual), 
                           weights.class0 = as.numeric(Y_all_kmers$is_true), curve = TRUE)

pr_Y_single_T <- pr.curve(scores.class0 = as.numeric(Y_single_T$mod_qual), 
                          weights.class0 = as.numeric(Y_single_T$is_true), curve = TRUE)

pr_m1Y_all_kmers <- pr.curve(scores.class0 = as.numeric(m1Y_all_kmers$mod_qual), 
                             weights.class0 = as.numeric(m1Y_all_kmers$is_true), curve = TRUE)

pr_m1Y_single_T <- pr.curve(scores.class0 = as.numeric(m1Y_single_T$mod_qual), 
                            weights.class0 = as.numeric(m1Y_single_T$is_true), curve = TRUE)

pr_m5U_all_kmers <- pr.curve(scores.class0 = as.numeric(m5U_all_kmers$mod_qual), 
                             weights.class0 = as.numeric(m5U_all_kmers$is_true), curve = TRUE)

pr_m5U_single_T <- pr.curve(scores.class0 = as.numeric(m5U_single_T$mod_qual), 
                            weights.class0 = as.numeric(m5U_single_T$is_true), curve = TRUE)


# Convert PR curve data into a dataframe for ggplot

pr_data <- list(
  pr_Y_all_kmers = pr_Y_all_kmers$curve,
  pr_Y_single_T = pr_Y_single_T$curve,
  pr_m1Y_all_kmers = pr_m1Y_all_kmers$curve,
  pr_m1Y_single_T = pr_m1Y_single_T$curve,
  pr_m5U_all_kmers = pr_m5U_all_kmers$curve,
  pr_m5U_single_T = pr_m5U_single_T$curve
)

# Assign colors matching your ROC plot
colors <- c("#E58112", "#F9CFA2", "#EE9A91", "#F9E0DB", "#82476F","#D7AEC6")


# Open a EPS device
postscript(file = here("..","__results_Gregor","curlcakes","per_read","hac","precision_recall","pseU_allmods_combined.eps"),
           width = 8, height = 8)

# Create an empty plot
plot(NA, xlim = c(0, 1), ylim = c(0.2, 1), xlab = "Recall", ylab = "Precision", 
     main = "Precision-Recall Curve", type = "n")

# Add a dashed black horizontal baseline at precision = 0.5 
abline(h = 0.5, col = "black", lty = 2, lwd = 3)  

# Add lines for each PR curve
i <- 1
for (curve_name in names(pr_data)) {
  lines(pr_data[[curve_name]][,1], pr_data[[curve_name]][,2], col = colors[i], lwd = 4)
  i <- i + 1
}

# Compute AUC-PR values
pr_auc_values <- sapply(pr_data, function(pr) pr.curve(scores.class0 = pr[,1], weights.class0 = pr[,2])$auc.integral)

legend("bottomleft", legend = c(
  paste("Y | NNTNN = ", round(pr_auc_values[1], 3)),
  paste("Y | VVTVV = ", round(pr_auc_values[2], 3)),
  paste("m1Y | NNTNN = ", round(pr_auc_values[3], 3)),
  paste("m1Y | VVTVV = ", round(pr_auc_values[4], 3)),
  paste("m5U | NNTNN = ", round(pr_auc_values[5], 3)),
  paste("m5U | VVTVV = ", round(pr_auc_values[6], 3))
), col = c("#E58112", "#F9CFA2", "#EE9A91", "#F9E0DB", "#82476F","#D7AEC6"), lwd = 4, cex = 1, title = "Area under the PR curve [AUC-PR]:")

# Close the PDF device
dev.off()




################
### m5C model ##
################

#####################################
# Task 1: Import and annotate Data  #
#####################################


UNM_rep1_m5C <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m5C", "pod5---bc_1_s.txt"), header = T, sep = "\t")

UNM_rep2_m5C <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m5C", "pod5---bc_2_s.txt"), header = T, sep = "\t")

UNM_rep3_m5C <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m5C", "pod5---bc_3_s.txt"), header = T, sep = "\t")

UNM_m5C <- rbind(UNM_rep1_m5C,UNM_rep2_m5C,UNM_rep3_m5C)

rm(UNM_rep1_m5C,UNM_rep2_m5C,UNM_rep3_m5C) 


m5C_m5C_model <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m5C", "pod5---bc_13_s.txt"), header = T, sep = "\t")

hm5C_m5C_model <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m5C", "pod5---bc_7_s.txt"), header = T, sep = "\t")

ac4C_m5C_model <- read.table(file = here("..", "__data","curlcakes", "raw", "hac_modkit_extract_full", "RNA241003_Pool1_rna004_130bps_hac@v5.1.0_m5C", "pod5---bc_4_s.txt"), header = T, sep = "\t")


filter_and_annotate_m5C <- function(df, mod_code, is_true){
  
  df <- df %>%
    # Remove potential "-" strand predictions
    filter(ref_mod_strand == "+") %>%
    # Filter out any incomplete 5mers (e.g: "GCAT-")
    filter(!grepl("-", query_kmer)) %>%
    # Filter out baseQ >= 10
    filter(base_qual >= 10) %>%
    # Filter for desired mod code
    filter(mod_code == !!mod_code) %>% 
    # Add ground-truth column
    mutate(is_true = is_true) %>% 
    # Annotate DRACH 
    mutate(is_single_C_5mer = case_when(
      str_detect(query_kmer, "^[AGT][AGT][C][AGT][AGT]$") ~ 1,
      T ~ 0
    ))
  
}

UNM_m5C <- filter_and_annotate_m5C(df = UNM_m5C, mod_code = "m", is_true = 0)

m5C_m5C_model <- filter_and_annotate_m5C(df = m5C_m5C_model, mod_code = "m", is_true = 1)

hm5C_m5C_model <- filter_and_annotate_m5C(df = hm5C_m5C_model, mod_code = "m", is_true = 1)

ac4C_m5C_model <- filter_and_annotate_m5C(df = ac4C_m5C_model, mod_code = "m", is_true = 1)


##########################
# Task 2: Data Wrangling #
##########################


# Task : Create a balanced datasets

# Set seed to make sampling reads reproducible

################
### All kmers ##
################

set.seed(123)

list_of_dfs <- list(UNM_m5C,m5C_m5C_model,hm5C_m5C_model,ac4C_m5C_model)

min_rows <- min(sapply(list_of_dfs, nrow))

UNM_m5C_final <- slice_sample(UNM_m5C, n = min_rows)

m5C_m5C_model_final <- slice_sample(m5C_m5C_model, n = min_rows)

hm5C_m5C_model_final <- slice_sample(hm5C_m5C_model, n = min_rows)

ac4C_m5C_model_final <- slice_sample(ac4C_m5C_model, n = min_rows)


m5C_all_kmers <- rbind(m5C_m5C_model_final, UNM_m5C_final)

hm5C_all_kmers <- rbind(hm5C_m5C_model_final, UNM_m5C_final)

ac4C_all_kmers <- rbind(ac4C_m5C_model_final, UNM_m5C_final)


rm(UNM_m5C_final, m5C_m5C_model_final, hm5C_m5C_model_final, ac4C_m5C_model_final)


###############################
### All kmers central C only ##
###############################


UNM_m5C_single_C_tmp <- UNM_m5C %>% 
  filter(is_single_C_5mer == 1)

m5C_m5C_model_single_C_tmp <- m5C_m5C_model %>% 
  filter(is_single_C_5mer == 1)

hm5C_m5C_model_single_C_tmp <- hm5C_m5C_model %>% 
  filter(is_single_C_5mer == 1)

ac4C_m5C_model_single_C_tmp <- ac4C_m5C_model %>% 
  filter(is_single_C_5mer == 1)



list_of_dfs <- list(UNM_m5C_single_C_tmp, m5C_m5C_model_single_C_tmp, hm5C_m5C_model_single_C_tmp, ac4C_m5C_model_single_C_tmp)

min_rows <- min(sapply(list_of_dfs, nrow))


UNM_m5C_single_C_final <- slice_sample(UNM_m5C_single_C_tmp, n = min_rows)

m5C_m5C_model_single_C_final <- slice_sample(m5C_m5C_model_single_C_tmp, n = min_rows)

hm5C_m5C_model_single_C_final <- slice_sample(hm5C_m5C_model_single_C_tmp, n = min_rows)

ac4C_m5C_model_single_C_final <- slice_sample(ac4C_m5C_model_single_C_tmp, n = min_rows)


m5C_m5C_model_single_C <- rbind(m5C_m5C_model_single_C_final,UNM_m5C_single_C_final)

hm5C_m5C_model_single_C <- rbind(hm5C_m5C_model_single_C_final,UNM_m5C_single_C_final)

ac4C_m5C_model_single_C <- rbind(ac4C_m5C_model_single_C_final,UNM_m5C_single_C_final)


rm(UNM_m5C_single_C_final,m5C_m5C_model_single_C_final,hm5C_m5C_model_single_C_final,ac4C_m5C_model_single_C_final)


###########################
## m5C all mods combined ##
###########################

### All context ###

roc_m5C_all_kmers <- roc(response = m5C_all_kmers$is_true, predictor = m5C_all_kmers$mod_qual)

roc_hm5C_all_kmers <- roc(response = hm5C_all_kmers$is_true, predictor = hm5C_all_kmers$mod_qual)

roc_ac4C_all_kmers <- roc(response = ac4C_all_kmers$is_true, predictor = ac4C_all_kmers$mod_qual)

### All context single C ###

roc_m5C_single_C <- roc(response = m5C_m5C_model_single_C$is_true, predictor = m5C_m5C_model_single_C$mod_qual)

roc_hm5C_single_C <- roc(response = hm5C_m5C_model_single_C$is_true, predictor = hm5C_m5C_model_single_C$mod_qual)

roc_ac4C_single_C <- roc(response = ac4C_m5C_model_single_C$is_true, predictor = ac4C_m5C_model_single_C$mod_qual)


# Open a PDF device
pdf(file = here("..","__results_Gregor","curlcakes","per_read","hac","ROC","m5C_allmods_combined.pdf"), 
    width = 8, height = 8)

# Plot the first ROC curve to set up the plot
plot(roc_m5C_all_kmers, col = "#FBD33D", main = "ROC Curve for Percentage Modified", lwd = 4)

# Add the remaining ROC curves using lines()
lines(roc_m5C_single_C, col = "#FEEEB9", lwd = 4)
lines(roc_hm5C_all_kmers, col = "#AD2D1F", lwd = 4)
lines(roc_hm5C_single_C, col = "#E8A399", lwd = 4)
lines(roc_ac4C_all_kmers, col = "#96A6D5", lwd = 4)
lines(roc_ac4C_single_C, col = "#E1E6F5", lwd = 4)


# Compute AUC values
auc_values <- c(
  auc(roc_m5C_all_kmers),
  auc(roc_m5C_single_C),
  auc(roc_hm5C_all_kmers),
  auc(roc_hm5C_single_C),
  auc(roc_ac4C_all_kmers),
  auc(roc_ac4C_single_C)
)

# Add a legend to identify curves and AUC values
legend("bottomright", legend = c(
  paste("m5C | NNCNN = ", round(auc_values[1], 3)),
  paste("m5C | DDCDD = ", round(auc_values[2], 3)),
  paste("hm5C | NNCNN = ", round(auc_values[3], 3)),
  paste("hm5C | DDCDD = ", round(auc_values[4], 3)),
  paste("ac4C | NNCNN = ", round(auc_values[5], 3)),
  paste("ac4C | DDCDD = ", round(auc_values[6], 3))
), col = c("#FBD33D", "#FEEEB9", "#AD2D1F", "#E8A399", "#96A6D5","#E1E6F5"), lwd = 4, cex = 1, title = "Area under the curve [AUC]:")


# Close the PDF device
dev.off()


###############
## PR curve ###
###############

#############################
## pseU all mods combined  ##
#############################

# Compute Precision-Recall curves for each dataset

pr_m5C_all_kmers <- pr.curve(scores.class0 = as.numeric(m5C_all_kmers$mod_qual), 
                             weights.class0 = as.numeric(m5C_all_kmers$is_true), curve = TRUE)

pr_roc_m5C_single_C <- pr.curve(scores.class0 = as.numeric(m5C_m5C_model_single_C$mod_qual), 
                                weights.class0 = as.numeric(m5C_m5C_model_single_C$is_true), curve = TRUE)

pr_hm5C_all_kmers <- pr.curve(scores.class0 = as.numeric(hm5C_all_kmers$mod_qual), 
                              weights.class0 = as.numeric(hm5C_all_kmers$is_true), curve = TRUE)

pr_roc_hm5C_single_C <- pr.curve(scores.class0 = as.numeric(hm5C_m5C_model_single_C$mod_qual), 
                                 weights.class0 = as.numeric(hm5C_m5C_model_single_C$is_true), curve = TRUE)

pr_ac4C_all_kmers <- pr.curve(scores.class0 = as.numeric(ac4C_all_kmers$mod_qual), 
                              weights.class0 = as.numeric(ac4C_all_kmers$is_true), curve = TRUE)

pr_roc_ac4C_single_C <- pr.curve(scores.class0 = as.numeric(ac4C_m5C_model_single_C$mod_qual), 
                                 weights.class0 = as.numeric(ac4C_m5C_model_single_C$is_true), curve = TRUE)



# Convert PR curve data into a dataframe for ggplot

pr_data <- list(
  pr_m5C_all_kmers = pr_m5C_all_kmers$curve,
  pr_roc_m5C_single_C = pr_roc_m5C_single_C$curve,
  pr_hm5C_all_kmers = pr_hm5C_all_kmers$curve,
  pr_roc_hm5C_single_C = pr_roc_hm5C_single_C$curve,
  pr_ac4C_all_kmers = pr_ac4C_all_kmers$curve,
  pr_roc_ac4C_single_C = pr_roc_ac4C_single_C$curve
)

# Assign colors matching your ROC plot
colors <- c("#FBD33D", "#FEEEB9", "#AD2D1F", "#E8A399", "#96A6D5","#E1E6F5")



# Open a EPS device
postscript(file = here("..","__results_Gregor","curlcakes","per_read","hac","precision_recall","m5C_allmods_combined.eps"),
           width = 8, height = 8)

# Create an empty plot
plot(NA, xlim = c(0, 1), ylim = c(0.2, 1), xlab = "Recall", ylab = "Precision", 
     main = "Precision-Recall Curve", type = "n")

# Add a dashed black horizontal baseline at precision = 0.5 
abline(h = 0.5, col = "black", lty = 2, lwd = 3) 

# Add lines for each PR curve
i <- 1
for (curve_name in names(pr_data)) {
  lines(pr_data[[curve_name]][,1], pr_data[[curve_name]][,2], col = colors[i], lwd = 4)
  i <- i + 1
}

# Compute AUC-PR values
pr_auc_values <- sapply(pr_data, function(pr) pr.curve(scores.class0 = pr[,1], weights.class0 = pr[,2])$auc.integral)

legend("bottomleft", legend = c(
  paste("m5C | NNCNN = ", round(pr_auc_values[1], 3)),
  paste("m5C | DDCDD = ", round(pr_auc_values[2], 3)),
  paste("hm5C | NNCNN = ", round(pr_auc_values[3], 3)),
  paste("hm5C | DDCDD = ", round(pr_auc_values[4], 3)),
  paste("ac4C | NNCNN = ", round(pr_auc_values[5], 3)),
  paste("ac4C | DDCDD = ", round(pr_auc_values[6], 3))
), col = c("#FBD33D", "#FEEEB9", "#AD2D1F", "#E8A399", "#96A6D5","#E1E6F5"), lwd = 4, cex = 1, title = "Area under the PR curve [AUC - PR]:")

# Close the PDF device
dev.off()



