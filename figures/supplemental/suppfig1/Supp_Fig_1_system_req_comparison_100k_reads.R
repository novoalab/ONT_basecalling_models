#########################################
### Comparison of System requirements ###
#########################################

######################
### Load Packages ####
######################
 
pkgs <- c("backports","tidyverse","here","skimr","dplyr", "ggplot2", "ggsci","ggforce","rstatix",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr", "ggridges")

lapply(pkgs, library, character.only = TRUE)


##########################
#### 1: Import Data ######
##########################

# Task 1.1: Read Raw Data from google sheets

library(googlesheets4)

raw <- read_sheet("https://docs.google.com/spreadsheets/d/10dQev5qroKDvTcZ479MzsZSRu3xP34irR_zVJvch8ss/edit?gid=1238298005#gid=1238298005", sheet = 4)


##########################
#### 2: Data Wrangling ###
##########################

###################
### Real - Time ###
###################

# Task 2.1: Perform statistical analysis

### Calculate statistics to add manually:

stat.test_real_time <- raw %>%
  t_test(realtime ~ demux_model) %>% 
  add_significance("p")
stat.test_real_time

### Add x.y position for plotting

stat.test_real_time <- stat.test_real_time %>%
  add_xy_position(x = "demux_model", step.increase = 0.4)
stat.test_real_time$y.position <- stat.test_real_time$y.position



##############################
#### 3: Data Visualisation ###
##############################

### Specify a custom theme

t <- theme(
  legend.title = element_blank(),
  legend.text = element_text( size = 15),
  legend.position = "none",
  axis.text = element_text(size=15),
  axis.text.x = element_text(size=15, angle = 35,hjust=0.9),
  axis.title=element_text(size=15),
  panel.grid.major.y = element_line(color = "grey90"),
  panel.grid.minor.y = element_line(color = "grey90"),
  panel.border = element_rect(colour = "black", fill = NA),
  strip.text.x = element_text(size = 15, face = "bold")
)

# Set plotting order

raw$demux_model <- factor(raw$demux_model, levels = c( "b04_RNA004","b96_RNA004" ))
raw$replicate <- factor(raw$replicate, levels = c( "1", "2", "3"))

# Task 3.1: Compare Real Time on a per job basis 

library(ggdist)
library(EnvStats)

  

ggplot(raw, aes(x = demux_model, y = realtime)) +
  
  geom_bar(aes(fill = demux_model),stat = "summary", fun = "mean", width=0.9) +
  
  stat_summary(aes(color = demux_model),fun.data = "mean_sd", geom = "errorbar", width=.35) +
  
  geom_point(aes(x = demux_model, y = realtime, fill = demux_model), size=4, shape = 21,
             color = "black",
             position=position_jitter(width = 0.4),
             show.legend = F) +

  stat_pvalue_manual(stat.test_real_time,  label = "p.signif", tip.length = .02, size = 3, coord.flip = F) +
  
  scale_y_continuous(limits = c(0,1.5)) +
  scale_fill_manual(values = c('b04_RNA004' = '#9FD1C1', 'b96_RNA004' = '#6AAC9E')) +
  scale_color_manual(values = c('b04_RNA004' = '#9FD1C1', 'b96_RNA004' = '#6AAC9E')) +
  
  labs(x= "", y = "Real-Time [minutes/1e5 reads]") +
  
  #coord_flip() +
  theme_pubr() + t

ggsave(filename = "realtime_barplot_per100k_.pdf",
       plot = last_plot(),
       path = here("..","__results_Gregor", "b96_RNA004"),
       width = 4,
       height = 8,
       units = "in")


# Task 1.2: Read Raw Data from google sheets

library(googlesheets4)

raw <- read_sheet("https://docs.google.com/spreadsheets/d/10dQev5qroKDvTcZ479MzsZSRu3xP34irR_zVJvch8ss/edit?gid=1238298005#gid=1238298005", sheet = 4)

raw$barcode <- factor(raw$barcode)

### Specify a custom theme

t <- theme(
  legend.title = element_blank(),
  legend.text = element_text( size = 15),
  legend.position = "none",
  axis.text = element_text(size=15),
  axis.text.x = element_text(size=15, angle = 35,hjust=0.9),
  axis.title=element_text(size=15),
  panel.grid.major.x = element_line(color = "grey90"),
  panel.grid.minor.x = element_line(color = "grey90"),
  panel.border = element_rect(colour = "black", fill = NA),
  strip.text.x = element_text(size = 15, face = "bold")
)

# Set plotting order

ggplot(raw, aes(x = reorder(barcode, desc(barcode)), y = `reads_bq>50`)) +
  
  geom_col() +
  
  labs(x= "", y = "number of reads [baseQ > 50]") +
  
  coord_flip() +
  theme_pubr() + t

ggsave(filename = "number_of_reads__barplot_per100k_.pdf",
       plot = last_plot(),
       path = here("..","__results_Gregor", "b96_RNA004"),
       width = 4,
       height = 8,
       units = "in")

