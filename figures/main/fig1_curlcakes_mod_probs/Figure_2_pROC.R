####################################################################
## Plotting pROC-curves on sampled modification probabilities  ##
####################################################################

# Task 0: Load packages

pkgs <- c("openintro","tidyverse","backports","here","skimr","dplyr", "ggplot2", "ggsci","ggforce",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr","scales")

lapply(pkgs, library, character.only = TRUE)

# Script specific package

library(pROC)

# Note: pROC was plotted solely for 0.5 ratio as Sensitivity vs. 1-Specificity is insensitive to class-imbalance hence redundant

##############################
# m6A_DRACH and inosine_m6A ##
##############################

#####################################
# Task 1: Import and annotate Data  #
#####################################

inosine_m6A <- read.table(file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","inosine_m6A_all_kmers.tsv"), header = T, sep = "\t") %>% 
  filter(ratio == 0.5)

m6A_DRACH <- read.table(file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","m6A_DRACH_all_kmers.tsv"), header = T, sep = "\t") %>% 
  filter(ratio == 0.5)

####################
# Task 3: Plot ROC #
####################

### inosine_m6A | NNANN ###

roc_inosine_m6A_NNANN <- inosine_m6A %>%
  roc(is_true ~ mod_qual, data = .)

### inosine_m6A | BBABB ###

roc_inosine_m6A_BBABB <- inosine_m6A %>%
  filter(is_single_A_5mer == 1) %>% 
  roc(is_true ~ mod_qual, data = .)

### inosine_m6A | DRACH ###

roc_inosine_m6A_DRACH <- inosine_m6A %>%
  filter(is_DRACH == 1) %>% 
  roc(is_true ~ mod_qual, data = .)

### inosine_m6A | KGACY ###

roc_inosine_m6A_KGACY <- inosine_m6A %>%
  filter(is_DRACH_single_A == 1) %>% 
  roc(is_true ~ mod_qual, data = .)

### m6A_DRACH | DRACH ### 

roc_m6A_DRACH_DRACH <- m6A_DRACH %>%
  filter(is_DRACH == 1) %>% 
  roc(is_true ~ mod_qual, data = .)

### m6A_DRACH | KGACY ###

roc_m6A_DRACH_KGACY <- m6A_DRACH %>%
  filter(is_DRACH_single_A == 1) %>% 
  roc(is_true ~ mod_qual, data = .)


########################################
## m6A_DRACH and inosine_m6A combined ##
########################################

# Open a PDF device
pdf(file = here("..","__results_Gregor","curlcakes","per_read","hac","ROC","tmp","m6A_DRACH_inosine_m6A_combined.pdf"), 
    width = 8, height = 8)

# Plot the first ROC curve to set up the plot
plot(roc_inosine_m6A_NNANN, col = "#299189", main = "ROC Curve for Percentage Modified", lwd = 4, legacy.axes = T)

# Add the remaining ROC curves using lines()
lines(roc_inosine_m6A_BBABB, col = "#67BEB4", lwd = 4)
lines(roc_inosine_m6A_DRACH, col = "#8CD5CA", lwd = 4)
lines(roc_inosine_m6A_KGACY, col = "#B3ECE0", lwd = 4)
lines(roc_m6A_DRACH_DRACH, col = "#381a61", lwd = 4)
lines(roc_m6A_DRACH_KGACY, col = "#E6D4FF", lwd = 4)


# Compute AUC values
auc_values <- c(
  auc(roc_inosine_m6A_NNANN),
  auc(roc_inosine_m6A_BBABB),
  auc(roc_inosine_m6A_DRACH),
  auc(roc_inosine_m6A_KGACY),
  auc(roc_m6A_DRACH_DRACH),
  auc(roc_m6A_DRACH_KGACY)
)

# Add a legend to identify curves and AUC values
legend("bottomright", legend = c(
  paste("inosine_m6A | NNANN = ", round(auc_values[1], 2)),
  paste("inosine_m6A | BBABB = ", round(auc_values[2], 2)),
  paste("inosine_m6A | DRACH = ", round(auc_values[3], 2)),
  paste("inosine_m6A | KGACY = ", round(auc_values[4], 2)),
  paste("m6A_DRACH | DRACH = ", round(auc_values[5], 2)),
  paste("m6A_DRACH | KGACY = ", round(auc_values[6], 2))
), col = c("#299189", "#67BEB4", "#8CD5CA", "#B3ECE0", "#381a61","#E6D4FF"), lwd = 4, cex = 1, title = "Area under the curve [AUC]:")


# Close the PDF device
dev.off()


#################
### pseU model ##
#################

########################
# Task 1: Import Data  #
########################


Y_pseU_all_kmers <- read.table(file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","pseU_all_kmers.tsv"), header = T, sep = "\t") %>% 
  filter(ratio == 0.5)

m1Y_pseU_all_kmers <- read.table(file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","m1Y_all_kmers.tsv"), header = T, sep = "\t") %>% 
  filter(ratio == 0.5)

m5U_pseU_all_kmers <- read.table(file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","m5U_all_kmers.tsv"), header = T, sep = "\t") %>% 
  filter(ratio == 0.5)


####################
# Task 3: Plot ROC #
####################

### Y | pseU | NNTNN ###

roc_Y_pseU_NNANN <- Y_pseU_all_kmers %>%
  roc(is_true ~ mod_qual, data = .)

### m1Y | pseU | NNTNN ###

roc_m1Y_pseU_NNANN <- m1Y_pseU_all_kmers %>%
  roc(is_true ~ mod_qual, data = .)

### m5U | pseU | NNTNN ###

roc_m5U_pseU_NNANN <- m5U_pseU_all_kmers %>%
  roc(is_true ~ mod_qual, data = .)

### Y | pseU | VVTVV ###

roc_Y_pseU_VVTVV <- Y_pseU_all_kmers %>%
  filter(is_single_T_5mer == 1) %>% 
  roc(is_true ~ mod_qual, data = .)

### m1Y | pseU | VVTVV ###

roc_m1Y_pseU_VVTVV <- m1Y_pseU_all_kmers %>%
  filter(is_single_T_5mer == 1) %>%
  roc(is_true ~ mod_qual, data = .)

### m5U | pseU | VVTVV ###

roc_m5U_pseU_VVTVV <- m5U_pseU_all_kmers %>%
  filter(is_single_T_5mer == 1) %>%
  roc(is_true ~ mod_qual, data = .)


############################
## pseU all mods combined ##
############################

# Open a PDF device
pdf(file = here("..","__results_Gregor","curlcakes","per_read","hac","ROC","tmp","pseU_allmods_combined.pdf"), 
    width = 8, height = 8)

# Plot the first ROC curve to set up the plot
plot(roc_Y_pseU_NNANN, col = "#E58112", main = "ROC Curve for Percentage Modified", lwd = 4)

# Add the remaining ROC curves using lines()
lines(roc_Y_pseU_VVTVV, col = "#F9CFA2", lwd = 4)
lines(roc_m1Y_pseU_NNANN, col = "#EE9A91", lwd = 4)
lines(roc_m1Y_pseU_VVTVV, col = "#F9E0DB", lwd = 4)
lines(roc_m5U_pseU_NNANN, col = "#82476F", lwd = 4)
lines(roc_m5U_pseU_VVTVV, col = "#D7AEC6", lwd = 4)


# Compute AUC values
auc_values <- c(
  auc(roc_Y_pseU_NNANN),
  auc(roc_Y_pseU_VVTVV),
  auc(roc_m1Y_pseU_NNANN),
  auc(roc_m1Y_pseU_VVTVV),
  auc(roc_m5U_pseU_NNANN),
  auc(roc_m5U_pseU_VVTVV)
)

# Add a legend to identify curves and AUC values
legend("bottomright", legend = c(
  paste("Y | NNTNN = ", round(auc_values[1], 2)),
  paste("Y | VVTVV = ", round(auc_values[2], 2)),
  paste("m1Y | NNTNN = ", round(auc_values[3], 2)),
  paste("m1Y | VVTVV = ", round(auc_values[4], 2)),
  paste("m5U | NNTNN = ", round(auc_values[5], 2)),
  paste("m5U | VVTVV = ", round(auc_values[6], 2))
), col = c("#E58112", "#F9CFA2", "#EE9A91", "#F9E0DB", "#82476F","#D7AEC6"), lwd = 4, cex = 1, title = "Area under the curve [AUC]:")


# Close the PDF device
dev.off()


################
### m5C model ##
################

#####################################
# Task 1: Import and annotate Data  #
#####################################

m5C_m5C_all_kmers <- read.table(file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","m5C_all_kmers.tsv"), header = T, sep = "\t") %>% 
  filter(ratio == 0.5)

hm5C_m5C_all_kmers <- read.table(file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","hm5C_model_all_kmers.tsv"), header = T, sep = "\t") %>% 
  filter(ratio == 0.5)

ac4C_m5C_all_kmers <- read.table(file = here("..","__data","curlcakes","processed","modkit_extract_full_tables_at_ratios","ac4C_model_all_kmers.tsv"), header = T, sep = "\t") %>% 
  filter(ratio == 0.5)


###########################
## m5C all mods combined ##
###########################

### m5C | m5C | NNCNN ###

roc_m5C_m5C_NNCNN <- m5C_m5C_all_kmers %>%
  roc(is_true ~ mod_qual, data = .)

### hm5C | hm5C | NNCNN ###

roc_hm5C_m5C_NNCNN <- hm5C_m5C_all_kmers %>%
  roc(is_true ~ mod_qual, data = .)

### ac4C | ac4C | NNCNN ###

roc_ac4C_m5C_NNCNN <- ac4C_m5C_all_kmers %>%
  roc(is_true ~ mod_qual, data = .)

### m5C | m5C | DDCDD ###

roc_m5C_m5C_DDCDD <- m5C_m5C_all_kmers %>%
  filter(is_single_C_5mer == 1) %>% 
  roc(is_true ~ mod_qual, data = .)

### hm5C | hm5C | DDCDD ###

roc_hm5C_m5C_DDCDD <- hm5C_m5C_all_kmers %>%
  filter(is_single_C_5mer == 1) %>% 
  roc(is_true ~ mod_qual, data = .)

### ac4C | ac4C | DDCDD ###

roc_ac4C_m5C_DDCDD <- ac4C_m5C_all_kmers %>%
  filter(is_single_C_5mer == 1) %>% 
  roc(is_true ~ mod_qual, data = .)


# Open a PDF device
pdf(file = here("..","__results_Gregor","curlcakes","per_read","hac","ROC","tmp","m5C_allmods_combined.pdf"), 
    width = 8, height = 8)

# Plot the first ROC curve to set up the plot
plot(roc_m5C_m5C_NNCNN, col = "#FBD33D", main = "ROC Curve for Percentage Modified", lwd = 4)

# Add the remaining ROC curves using lines()
lines(roc_m5C_m5C_DDCDD, col = "#FEEEB9", lwd = 4)
lines(roc_hm5C_m5C_NNCNN, col = "#AD2D1F", lwd = 4)
lines(roc_hm5C_m5C_DDCDD, col = "#E8A399", lwd = 4)
lines(roc_ac4C_m5C_NNCNN, col = "#96A6D5", lwd = 4)
lines(roc_ac4C_m5C_DDCDD, col = "#E1E6F5", lwd = 4)


# Compute AUC values
auc_values <- c(
  auc(roc_m5C_m5C_NNCNN),
  auc(roc_m5C_m5C_DDCDD),
  auc(roc_hm5C_m5C_NNCNN),
  auc(roc_hm5C_m5C_DDCDD),
  auc(roc_ac4C_m5C_NNCNN),
  auc(roc_ac4C_m5C_DDCDD)
)

# Add a legend to identify curves and AUC values
legend("bottomright", legend = c(
  paste("m5C | NNCNN = ", round(auc_values[1], 2)),
  paste("m5C | DDCDD = ", round(auc_values[2], 2)),
  paste("hm5C | NNCNN = ", round(auc_values[3], 2)),
  paste("hm5C | DDCDD = ", round(auc_values[4], 2)),
  paste("ac4C | NNCNN = ", round(auc_values[5], 2)),
  paste("ac4C | DDCDD = ", round(auc_values[6], 2))
), col = c("#FBD33D", "#FEEEB9", "#AD2D1F", "#E8A399", "#96A6D5","#E1E6F5"), lwd = 4, cex = 1, title = "Area under the curve [AUC]:")


# Close the PDF device
dev.off()
