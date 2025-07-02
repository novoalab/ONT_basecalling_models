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

master_df <- read.table(file = here("..", "__data","curlcakes","processed","hac_0.1_filter-pass","master_table.tsv"), header = T, sep = "\t")


#################################
# Task 2: Plot Denisty profiles #
#################################

# Task 2.0: Save a custom theme

t <- theme(
  legend.title = element_blank(),
  legend.text = element_text( size = 15),
  legend.position = "none",
  axis.text = element_text(size=15),
  axis.title=element_text(size=15),
  panel.grid.major.y = element_line(color = "grey90"),
  panel.grid.minor.y = element_line(color = "grey90"),
  panel.border = element_rect(colour = "black", fill = NA),
  strip.text.x = element_text(size = 15, face = "bold")
)

# Task 2.2: Generate index to iterate over

list_of_plots  <- list(c("RNA241003_Pool1_rna004_130bps","A","m6A_DRACH","m6A_100",'#381a61'),
                       c("RNA241003_Pool1_rna004_130bps","A","m6A_DRACH","m6A_75",'#381a61'),
                       c("RNA241003_Pool1_rna004_130bps","A","m6A_DRACH","m6A_50",'#381a61'),
                       c("RNA241003_Pool1_rna004_130bps","A","m6A_DRACH","m6A_25",'#381a61'),
                       c("RNA241003_Pool1_rna004_130bps","A","m6A_DRACH","m6A_12.5",'#381a61'),
                       
                       c("RNA241003_Pool1_rna004_130bps","A","inosine_m6A-a","m6A_100",'#381a61'),
                       c("RNA241003_Pool1_rna004_130bps","A","inosine_m6A-a","m6A_75",'#381a61'),
                       c("RNA241003_Pool1_rna004_130bps","A","inosine_m6A-a","m6A_50",'#381a61'),
                       c("RNA241003_Pool1_rna004_130bps","A","inosine_m6A-a","m6A_25",'#381a61'),
                       c("RNA241003_Pool1_rna004_130bps","A","inosine_m6A-a","m6A_12.5",'#381a61'),
                       
                       c("RNA241003_Pool1_rna004_130bps","A","inosine_m6A-17596","m6A_100",'#381a61'),
                       c("RNA241003_Pool1_rna004_130bps","A","inosine_m6A-17596","m6A_75",'#381a61'),
                       c("RNA241003_Pool1_rna004_130bps","A","inosine_m6A-17596","m6A_50",'#381a61'),
                       c("RNA241003_Pool1_rna004_130bps","A","inosine_m6A-17596","m6A_25",'#381a61'),
                       c("RNA241003_Pool1_rna004_130bps","A","inosine_m6A-17596","m6A_12.5",'#381a61'),
                       
                       c("RNA241003_Pool1_rna004_130bps","T","pseU","Y_100",'#e78429'),
                       c("RNA241003_Pool1_rna004_130bps","T","pseU","Y_75",'#e78429'),
                       c("RNA241003_Pool1_rna004_130bps","T","pseU","Y_50",'#e78429'),
                       c("RNA241003_Pool1_rna004_130bps","T","pseU","Y_25",'#e78429'),
                       c("RNA241003_Pool1_rna004_130bps","T","pseU","Y_12.5",'#e78429'),
                       c("RNA241003_Pool1_rna004_130bps","T","pseU","m1Y",'#ed968c'),
                       c("RNA241016_Pool3_rna004_130bps","T","pseU","m5U_100",'#7c4b73'),
                       
                       c("RNA241003_Pool1_rna004_130bps","C","m5C","m5C",'#f9d14a'),
                       c("RNA241003_Pool1_rna004_130bps","C","m5C","hm5C","#ab3329"),
                       c("RNA241003_Pool1_rna004_130bps","C","m5C","ac4C","#88a0dc"),
                       
                       c("RNA241003_Pool1_rna004_130bps","A","m6A_DRACH","UNM","grey75"),
                       c("RNA241003_Pool1_rna004_130bps","A","inosine_m6A-a","UNM","grey75"),
                       c("RNA241003_Pool1_rna004_130bps","A","inosine_m6A-17596","UNM","grey75"),
                       c("RNA241003_Pool1_rna004_130bps","T","pseU","UNM","grey75"),
                       c("RNA241003_Pool1_rna004_130bps","C","m5C","UNM","grey75"))
                       
                    
########################################
# Task 3: Plot Boxplots / Half Density #
########################################

# Task 3.1: Filter for m6A and factorize

m6A_df <- master_df %>% 
  filter(sample == "RNA241003_Pool1_rna004_130bps") %>% 
  filter(model == "inosine_m6A-a") %>% 
  filter(base == "A") %>% 
  filter(modification %in% c("m6A_100","m6A_75","m6A_50","m6A_25","m6A_12.5", "UNM")) %>% 
  filter(coverage >= 5)

m6A_df$modification <- factor(m6A_df$modificatio, levels = c("m6A_100","m6A_75","m6A_50","m6A_25","m6A_12.5","UNM"))

# Calculate median per group
medians <- m6A_df %>%
  group_by(modification) %>%
  summarize(med = median(stoichiometry, na.rm = TRUE))


### Flipped boxplot of individual errors

library(ggdist)
library(EnvStats)

t <- theme(
  legend.title = element_blank(),
  legend.text = element_text( size = 15),
  legend.position = "none",
  axis.text = element_text(size=15),
  axis.title=element_text(size=15),
  panel.grid.major.x = element_line(color = "grey90"),
  panel.grid.minor.x = element_line(color = "grey90"),
  panel.border = element_rect(colour = "black", fill = NA),
  strip.text.x = element_text(size = 15, face = "bold")
)

m6A_df %>%
  ggplot(aes(x = modification, y = stoichiometry)) +
  # Stat Boxplot to add errorbar optic
  stat_boxplot(geom = "errorbar",
               lwd = 0.4,
               width = 0.1) +
  # Central boxplot
  geom_boxplot(aes(fill = modification),
               width=0.25, color="grey20",position = position_dodge(width =0.95),
               alpha = 1, lwd = 0.4, notch = T, outlier.shape = NA) +
  # Distribution as half plot
  stat_halfeye(aes(fill = modification), normalize = "xy", slab_colour = "black", slab_size = 0.5,
               alpha = 0.75, adjust = .5, width = .3, .width = 0, justification = -.7, point_colour = NA, show.legend = F) +
  scale_fill_manual(values = c('m6A_100' = '#381a61', 'm6A_75' = '#553487', 'm6A_50' = "#7150AD", 'm6A_25' = "#9A80D1", 'm6A_12.5' = "#D8D0ED","UNM" = "grey50")) +
  # Include median
  geom_text(data = medians, aes(x = modification, y = med, label = round(med, 1)), 
            vjust = -0.5, size = 4.5, color = "black") +
  scale_x_discrete(labels= c("100.0","75.0","50.0","25.0","12.5","0")) +
  labs(x="true stoichiometry [%]",y=expression("predicted stoichiometry [%]")) +
  stat_n_text(size = 5) +
  theme_pubr() + t +
  theme(axis.text.x = element_text(size = 15),
        legend.text = element_text( size = 15)) + 
  coord_flip()

ggsave(filename = here("..","__results_Gregor","curlcakes","per_position","boxplots_combined_stoichiometries_0.1_filter_pass","inosine_m6A_stoichiometries_combined_stat_n_median.pdf"),       
       plot = last_plot(), 
       device = "pdf", 
       width = 4, height = 8, units = "in")

    

# Task 3.1: Filter for m6A and factorize

m6A_df <- master_df %>% 
  filter(sample == "RNA241003_Pool1_rna004_130bps") %>% 
  filter(model == "m6A_DRACH") %>% 
  filter(base == "A") %>% 
  filter(modification %in% c("m6A_100","m6A_75","m6A_50","m6A_25","m6A_12.5", "UNM")) %>% 
  filter(coverage >= 5)

m6A_df$modification <- factor(m6A_df$modificatio, levels = c("m6A_100","m6A_75","m6A_50","m6A_25","m6A_12.5","UNM"))

# Calculate median per group
medians <- m6A_df %>%
  group_by(modification) %>%
  summarize(med = median(stoichiometry, na.rm = TRUE))


### Flipped boxplot of individual errors

library(ggdist)
library(EnvStats)

t <- theme(
  legend.title = element_blank(),
  legend.text = element_text( size = 15),
  legend.position = "none",
  axis.text = element_text(size=15),
  axis.title=element_text(size=15),
  panel.grid.major.x = element_line(color = "grey90"),
  panel.grid.minor.x = element_line(color = "grey90"),
  panel.border = element_rect(colour = "black", fill = NA),
  strip.text.x = element_text(size = 15, face = "bold")
)

m6A_df %>%
  ggplot(aes(x = modification, y = stoichiometry)) +
  # Stat Boxplot to add errorbar optic
  stat_boxplot(geom = "errorbar",
               lwd = 0.4,
               width = 0.1) +
  # Central boxplot
  geom_boxplot(aes(fill = modification),
               width=0.25, color="grey20",position = position_dodge(width =0.95),
               alpha = 1, lwd = 0.4, notch = T, outlier.shape = NA) +
  # Distribution as half plot
  stat_halfeye(aes(fill = modification), normalize = "xy", slab_colour = "black", slab_size = 0.5,
               alpha = 0.75, adjust = .5, width = .3, .width = 0, justification = -.7, point_colour = NA, show.legend = F) +
  scale_fill_manual(values = c('m6A_100' = '#381a61', 'm6A_75' = '#553487', 'm6A_50' = "#7150AD", 'm6A_25' = "#9A80D1", 'm6A_12.5' = "#D8D0ED","UNM" = "grey50")) +
  # Include median
  geom_text(data = medians, aes(x = modification, y = med, label = round(med, 1)), 
            vjust = -0.5, size = 4.5, color = "black") +
  scale_x_discrete(labels= c("100.0","75.0","50.0","25.0","12.5","0")) +
  labs(x="true stoichiometry [%]",y=expression("predicted stoichiometry [%]")) +
  stat_n_text(size = 5) +
  theme_pubr() + t +
  theme(axis.text.x = element_text(size = 15),
        legend.text = element_text( size = 15)) + 
  coord_flip()

ggsave(filename = here("..","__results_Gregor","curlcakes","per_position","boxplots_combined_stoichiometries_0.1_filter_pass","m6A_DRACH_stoichiometries_combined_stat_n_median.pdf"),       
       plot = last_plot(), 
       device = "pdf", 
       width = 4, height = 8, units = "in")






# Task 3.1: Filter for pseU and factorize

pseU_df <- master_df %>% 
  filter(sample == "RNA241003_Pool1_rna004_130bps") %>% 
  filter(model == "pseU") %>% 
  filter(base == "T") %>% 
  filter(modification %in% c("Y_100","Y_75","Y_50","Y_25","Y_12.5", "UNM")) %>% 
  filter(coverage >= 5)

pseU_df$modification <- factor(pseU_df$modification, levels = c("Y_100","Y_75","Y_50","Y_25","Y_12.5", "UNM"))

# Calculate median per group
medians <- pseU_df %>%
  group_by(modification) %>%
  summarize(med = median(stoichiometry, na.rm = TRUE))


### Flipped boxplot of individual errors

library(ggdist)
library(EnvStats)

t <- theme(
  legend.title = element_blank(),
  legend.text = element_text( size = 15),
  legend.position = "none",
  axis.text = element_text(size=15),
  axis.title=element_text(size=15),
  panel.grid.major.x = element_line(color = "grey90"),
  panel.grid.minor.x = element_line(color = "grey90"),
  panel.border = element_rect(colour = "black", fill = NA),
  strip.text.x = element_text(size = 15, face = "bold")
)

pseU_df %>%
  ggplot(aes(x = modification, y = stoichiometry)) +
  # Stat Boxplot to add errorbar optic
  stat_boxplot(geom = "errorbar",
               lwd = 0.4,
               width = 0.1) +
  # Central boxplot
  geom_boxplot(aes(fill = modification),
               width=0.25, color="grey20",position = position_dodge(width =0.95),
               alpha = 1, lwd = 0.4, notch = T, outlier.shape = NA) +
  # Distribution as half plot
  stat_halfeye(aes(fill = modification), normalize = "xy", slab_colour = "black", slab_size = 0.5,
               alpha = 0.75, adjust = .5, width = .3, .width = 0, justification = -.7, point_colour = NA, show.legend = F) +
  scale_fill_manual(values = c('Y_100' = '#e78429', 'Y_75' = '#ec9b50', 'Y_50' = "#f1b275", 'Y_25' = "#f6c99c", 'Y_12.5' = "#fbe1c2","UNM" = "grey50")) +
  # Include median
  geom_text(data = medians, aes(x = modification, y = med, label = round(med, 1)), 
            vjust = -0.5, size = 4.5, color = "black") +
  scale_x_discrete(labels= c("100.0","75.0","50.0","25.0","12.5","0")) +
  labs(x="true stoichiometry [%]",y=expression("predicted stoichiometry [%]")) +
  stat_n_text(size = 5) +
  theme_pubr() + t +
  theme(axis.text.x = element_text(size = 15),
        legend.text = element_text( size = 15)) + 
  coord_flip()

ggsave(filename = here("..","__results_Gregor","curlcakes","per_position","boxplots_combined_stoichiometries_0.1_filter_pass","pseU_stoichiometries_combined_stat_n_median.pdf"),       
       plot = last_plot(), 
       device = "pdf", 
       width = 4, height = 8, units = "in")


# Task 3.1: Filter for pseU and factorize

m5C_df <- master_df %>% 
  filter(sample == "RNA241003_Pool1_rna004_130bps") %>% 
  filter(model == "m5C") %>% 
  filter(base == "C") %>% 
  filter(modification %in% c("m5C","UNM")) %>% 
  filter(coverage >= 5)

m5C_df$modification <- factor(m5C_df$modification, levels = c("m5C", "UNM"))

# Calculate median per group
medians <- m5C_df %>%
  group_by(modification) %>%
  summarize(med = median(stoichiometry, na.rm = TRUE))


### Flipped boxplot of individual errors

library(ggdist)
library(EnvStats)

t <- theme(
  legend.title = element_blank(),
  legend.text = element_text( size = 15),
  legend.position = "none",
  axis.text = element_text(size=15),
  axis.title=element_text(size=15),
  panel.grid.major.x = element_line(color = "grey90"),
  panel.grid.minor.x = element_line(color = "grey90"),
  panel.border = element_rect(colour = "black", fill = NA),
  strip.text.x = element_text(size = 15, face = "bold")
)

m5C_df %>%
  ggplot(aes(x = modification, y = stoichiometry)) +
  # Stat Boxplot to add errorbar optic
  stat_boxplot(geom = "errorbar",
               lwd = 0.4,
               width = 0.1) +
  # Central boxplot
  geom_boxplot(aes(fill = modification),
               width=0.25, color="grey20",position = position_dodge(width =0.95),
               alpha = 1, lwd = 0.4, notch = T, outlier.shape = NA) +
  # Distribution as half plot
  stat_halfeye(aes(fill = modification), normalize = "xy", slab_colour = "black", slab_size = 0.5,
               alpha = 0.75, adjust = .5, width = .3, .width = 0, justification = -.7, point_colour = NA, show.legend = F) +
  scale_fill_manual(values = c('m5C' = '#f9d14a',"UNM" = "grey50")) +
  # Include median
  geom_text(data = medians, aes(x = modification, y = med, label = round(med, 1)), 
            vjust = -0.5, size = 4.5, color = "black") +
  scale_x_discrete(labels= c("100.0","75.0","50.0","25.0","12.5","0")) +
  labs(x="true stoichiometry [%]",y=expression("predicted stoichiometry [%]")) +
  stat_n_text(size = 5) +
  theme_pubr() + t +
  theme(axis.text.x = element_text(size = 15),
        legend.text = element_text( size = 15)) + 
  coord_flip()

ggsave(filename = here("..","__results_Gregor","curlcakes","per_position","boxplots_combined_stoichiometries_0.1_filter_pass","m5C_stoichiometries_combined_stat_n_median.pdf"),       
       plot = last_plot(), 
       device = "pdf", 
       width = 4, height = 4, units = "in")

