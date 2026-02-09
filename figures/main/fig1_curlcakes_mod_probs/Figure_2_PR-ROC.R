###################################################################
## Plotting PR-ROC curves on sampled modification probabilities  ##
###################################################################

# Task 0: Load packages

pkgs <- c("openintro","tidyverse","backports","here","skimr","dplyr", "ggplot2", "ggsci","ggforce",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr","scales")

lapply(pkgs, library, character.only = TRUE)

# Script specific package

library(PRROC)

# Note: pROC was plotted solely for 0.5 ratio as Sensitivity vs. 1-Specificity is insensitive to class-imbalance hence redundant

##############################
# m6A_DRACH and inosine_m6A ##
##############################

#####################################
# Task 1: Import and annotate Data  #
#####################################

inosine_m6A <- read.table(file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","inosine_m6A_all_kmers.tsv"), header = T, sep = "\t")

m6A_DRACH <- read.table(file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","m6A_DRACH_all_kmers.tsv"), header = T, sep = "\t") 

#######################
# Task 3: Plot PR-ROC #
#######################

###################
### Inosine_m6A ###
###################

# Compute Precision-Recall curves for each dataset

pr_inosine_m6A_0.5 <- inosine_m6A %>%
  filter(ratio == 0.5) %>%
  filter(is_single_A_5mer == 1) %>% 
  { pr.curve(
    scores.class0 = .$mod_qual,
    weights.class0 = .$is_true,
    curve = TRUE
  ) }

pr_inosine_m6A_0.1 <- inosine_m6A %>%
  filter(ratio == 0.1) %>%
  filter(is_single_A_5mer == 1) %>% 
  { pr.curve(
    scores.class0 = .$mod_qual,
    weights.class0 = .$is_true,
    curve = TRUE
  ) }

pr_inosine_m6A_0.01 <- inosine_m6A %>%
  filter(ratio == 0.01) %>%
  filter(is_single_A_5mer == 1) %>% 
  { pr.curve(
    scores.class0 = .$mod_qual,
    weights.class0 = .$is_true,
    curve = TRUE
  ) }

pr_inosine_m6A_0.001 <- inosine_m6A %>%
  filter(ratio == 0.001) %>%
  filter(is_single_A_5mer == 1) %>%  
  { pr.curve(
    scores.class0 = .$mod_qual,
    weights.class0 = .$is_true,
    curve = TRUE
  ) }

# Convert PR curve data into a dataframe for ggplot

pr_data <- list(
  pr_inosine_m6A_0.5 = pr_inosine_m6A_0.5$curve,
  pr_inosine_m6A_0.1 = pr_inosine_m6A_0.1$curve,
  pr_inosine_m6A_0.01 = pr_inosine_m6A_0.01$curve,
  pr_inosine_m6A_0.001 = pr_inosine_m6A_0.001$curve
)


# Assign colors matching your ROC plot
colors <- c("#381a61", "#7150AD", "#9A80D1", "#D8D0ED")

# Open a EPS device
postscript(file = here("..","__results_Gregor","curlcakes","per_read","hac","precision_recall", "tmp","inosine_m6A_combined_all_kmers.eps"),
           width = 8, height = 10)

# Create an empty plot
plot(NA, xlim = c(0, 1), ylim = c(0.0, 1), xlab = "Recall", ylab = "Precision", 
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
  paste("m6A/A = 0.5 |", round(pr_auc_values[1], 2)),
  paste("m6A/A = 0.1 |", round(pr_auc_values[2], 2)),
  paste("m6A/A = 0.01 |", round(pr_auc_values[3], 2)),
  paste("m6A/A = 0.001 |", round(pr_auc_values[4], 2))
), col = colors, lwd = 4, cex = 1, title = "Area under the PR curve [AUC-PR]:")

# Close the PDF device
dev.off()


###################
### m6A_DRACH #####
###################

# Compute Precision-Recall curves for each dataset

pr_m6A_DRACH_0.5 <- m6A_DRACH %>%
  filter(ratio == 0.5) %>%
  #filter(is_DRACH_single_A == 1) %>% 
  { pr.curve(
    scores.class0 = .$mod_qual,
    weights.class0 = .$is_true,
    curve = TRUE
  ) }

pr_m6A_DRACH_0.1 <- m6A_DRACH %>%
  filter(ratio == 0.1) %>%
  #filter(is_DRACH_single_A == 1) %>% 
  { pr.curve(
    scores.class0 = .$mod_qual,
    weights.class0 = .$is_true,
    curve = TRUE
  ) }

pr_m6A_DRACH_0.01 <- m6A_DRACH %>%
  filter(ratio == 0.01) %>%
  #filter(is_DRACH_single_A == 1) %>% 
  { pr.curve(
    scores.class0 = .$mod_qual,
    weights.class0 = .$is_true,
    curve = TRUE
  ) }

pr_m6A_DRACH_0.001 <- m6A_DRACH %>%
  filter(ratio == 0.001) %>%
  #filter(is_DRACH_single_A == 1) %>% 
  { pr.curve(
    scores.class0 = .$mod_qual,
    weights.class0 = .$is_true,
    curve = TRUE
  ) }

# Convert PR curve data into a dataframe for ggplot

pr_data <- list(
  pr_m6A_DRACH_0.5 = pr_m6A_DRACH_0.5$curve,
  pr_m6A_DRACH_0.1 = pr_m6A_DRACH_0.1$curve,
  pr_m6A_DRACH_0.01 = pr_m6A_DRACH_0.01$curve,
  pr_m6A_DRACH_0.001 = pr_m6A_DRACH_0.001$curve
)


# Assign colors matching your ROC plot
colors <- c("#381a61", "#7150AD", "#9A80D1", "#D8D0ED")

# Open a EPS device
postscript(file = here("..","__results_Gregor","curlcakes","per_read","hac","precision_recall", "tmp","m6A_DRACH_combined_all_kmers_test.eps"),
           width = 8, height = 10)

# Create an empty plot
plot(NA, xlim = c(0, 1), ylim = c(0.0, 1), xlab = "Recall", ylab = "Precision", 
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
  paste("m6A/A = 0.5 |", round(pr_auc_values[1], 2)),
  paste("m6A/A = 0.1 |", round(pr_auc_values[2], 2)),
  paste("m6A/A = 0.01 |", round(pr_auc_values[3], 2)),
  paste("m6A/A = 0.001 |", round(pr_auc_values[4], 2))
), col = colors, lwd = 4, cex = 1, title = "Area under the PR curve [AUC-PR]:")

# Close the PDF device
dev.off()



#########
# pseU ##
#########

#####################################
# Task 1: Import and annotate Data  #
#####################################

Y_pseU_all_kmers <- read.table(file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","pseU_all_kmers.tsv"), header = T, sep = "\t")

#######################
# Task 3: Plot PR-ROC #
#######################

# Compute Precision-Recall curves for each dataset

pr_Y_pseU_0.5 <- Y_pseU_all_kmers %>%
  filter(ratio == 0.5) %>%
  filter(is_single_T_5mer == 1) %>% 
  { pr.curve(
    scores.class0 = .$mod_qual,
    weights.class0 = .$is_true,
    curve = TRUE
  ) }

pr_Y_pseU_0.1 <- Y_pseU_all_kmers %>%
  filter(ratio == 0.1) %>%
  filter(is_single_T_5mer == 1) %>% 
  { pr.curve(
    scores.class0 = .$mod_qual,
    weights.class0 = .$is_true,
    curve = TRUE
  ) }

pr_Y_pseU_0.01 <- Y_pseU_all_kmers %>%
  filter(ratio == 0.01) %>%
  filter(is_single_T_5mer == 1) %>% 
  { pr.curve(
    scores.class0 = .$mod_qual,
    weights.class0 = .$is_true,
    curve = TRUE
  ) }

pr_Y_pseU_0.001 <- Y_pseU_all_kmers %>%
  filter(ratio == 0.001) %>%
  filter(is_single_T_5mer == 1) %>% 
  { pr.curve(
    scores.class0 = .$mod_qual,
    weights.class0 = .$is_true,
    curve = TRUE
  ) }

# Convert PR curve data into a dataframe for ggplot

pr_data <- list(
  pr_Y_pseU_0.5 = pr_Y_pseU_0.5$curve,
  pr_Y_pseU_0.1 = pr_Y_pseU_0.1$curve,
  pr_Y_pseU_0.01 = pr_Y_pseU_0.01$curve,
  pr_Y_pseU_0.001 = pr_Y_pseU_0.001$curve
)


# Assign colors matching your ROC plot
colors <- c("#e78429", "#f1b275", "#f6c99c", "#fbe1c2")

# Open a EPS device
postscript(file = here("..","__results_Gregor","curlcakes","per_read","hac","precision_recall", "tmp","pseU_combined_single_central_kmer.eps"),
           width = 8, height = 10)

# Create an empty plot
plot(NA, xlim = c(0, 1), ylim = c(0.0, 1), xlab = "Recall", ylab = "Precision", 
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
  paste("Y/U = 0.5 |", round(pr_auc_values[1], 2)),
  paste("Y/U = 0.1 |", round(pr_auc_values[2], 2)),
  paste("Y/U = 0.01 |", round(pr_auc_values[3], 2)),
  paste("Y/U = 0.001 |", round(pr_auc_values[4], 2))
), col = colors, lwd = 4, cex = 1, title = "Area under the PR curve [AUC-PR]:")

# Close the PDF device
dev.off()


########
# m5C ##
########

#####################################
# Task 1: Import and annotate Data  #
#####################################

m5C_m5C_all_kmers <- read.table(file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","m5C_all_kmers.tsv"), header = T, sep = "\t")

#######################
# Task 3: Plot PR-ROC #
#######################

# Compute Precision-Recall curves for each dataset

pr_m5C_m5C_0.5 <- m5C_m5C_all_kmers %>%
  filter(is_single_C_5mer == 1) %>% 
  filter(ratio == 0.5) %>%
  { pr.curve(
    scores.class0 = .$mod_qual,
    weights.class0 = .$is_true,
    curve = TRUE
  ) }

pr_m5C_m5C_0.1 <- m5C_m5C_all_kmers %>%
  filter(is_single_C_5mer == 1) %>% 
  filter(ratio == 0.1) %>%
  { pr.curve(
    scores.class0 = .$mod_qual,
    weights.class0 = .$is_true,
    curve = TRUE
  ) }

pr_m5C_m5C_0.01 <- m5C_m5C_all_kmers %>%
  filter(ratio == 0.01) %>%
  { pr.curve(
    scores.class0 = .$mod_qual,
    weights.class0 = .$is_true,
    curve = TRUE
  ) }

pr_m5C_m5C_0.001 <- m5C_m5C_all_kmers %>%
  filter(is_single_C_5mer == 1) %>% 
  filter(ratio == 0.001) %>%
  { pr.curve(
    scores.class0 = .$mod_qual,
    weights.class0 = .$is_true,
    curve = TRUE
  ) }

# Convert PR curve data into a dataframe for ggplot

pr_data <- list(
  pr_m5C_m5C_0.5 = pr_m5C_m5C_0.5$curve,
  pr_m5C_m5C_0.1 = pr_m5C_m5C_0.1$curve,
  pr_m5C_m5C_0.01 = pr_m5C_m5C_0.01$curve,
  pr_m5C_m5C_0.001 = pr_m5C_m5C_0.001$curve
)


# Assign colors matching your ROC plot
colors <- c("#f9d14a", "#fadb66", "#fde47f", "#ffed99")

# Open a EPS device
postscript(file = here("..","__results_Gregor","curlcakes","per_read","hac","precision_recall", "tmp","m5C_combined_single_central_kmer.eps"),
           width = 8, height = 10)

# Create an empty plot
plot(NA, xlim = c(0, 1), ylim = c(0.0, 1), xlab = "Recall", ylab = "Precision", 
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
  paste("m5C/C = 0.5 |", round(pr_auc_values[1], 2)),
  paste("m5C/C = 0.1 |", round(pr_auc_values[2], 2)),
  paste("m5C/C = 0.01 |", round(pr_auc_values[3], 2)),
  paste("m5C/C = 0.001 |", round(pr_auc_values[4], 2))
), col = colors, lwd = 4, cex = 1, title = "Area under the PR curve [AUC-PR]:")

# Close the PDF device
dev.off()

