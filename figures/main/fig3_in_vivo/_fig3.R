### The input files needed for this script are found in the fig3_in_vivo folder




### Panel B - replicability


### S. cerevisiae

# Load libraries
library(ggplot2)
library(dplyr)
library(readr)

# File paths
rep1_file <- "BY4741-1_rep1.bedgraph"
rep2_file <- "BY4741-1_rep2.bedgraph"

# Read data
rep1 <- read_tsv(rep1_file, col_names = FALSE, col_types = cols()) %>%
  setNames(c("rRNA", "start", "end", "stoich_rep1")) %>%
  mutate(position = paste(rRNA, end, sep = ":"))

rep2 <- read_tsv(rep2_file, col_names = FALSE, col_types = cols()) %>%
  setNames(c("rRNA", "start", "end", "stoich_rep2")) %>%
  mutate(position = paste(rRNA, end, sep = ":"))

# Merge and calculate R²
merged_data <- inner_join(rep1, rep2, by = "position")
r_squared <- round(summary(lm(stoich_rep2 ~ stoich_rep1, data = merged_data))$r.squared, 3)

# Plot
plot <- ggplot(merged_data, aes(x = stoich_rep1, y = stoich_rep2)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  coord_fixed() +
  labs(
    title = "Stoichiometry Comparison Between Replicates (BW25113_R1)",
    x = "Replicate 1 Stoichiometry (%)",
    y = "Replicate 2 Stoichiometry (%)"
  ) +
  annotate("text", x = Inf, y = -Inf, label = paste0("R² = ", r_squared),
           hjust = 1.1, vjust = -0.5, size = 4, color = "black") +
  theme_minimal()

plot

### E. coli


# File paths
rep1_file <- "BW25113_R1_rep1.bedgraph"
rep2_file <- "BW25113_R1_rep2.bedgraph"

# Read data
rep1 <- read_tsv(rep1_file, col_names = FALSE, col_types = cols()) %>%
  setNames(c("rRNA", "start", "end", "stoich_rep1")) %>%
  mutate(position = paste(rRNA, end, sep = ":"))

rep2 <- read_tsv(rep2_file, col_names = FALSE, col_types = cols()) %>%
  setNames(c("rRNA", "start", "end", "stoich_rep2")) %>%
  mutate(position = paste(rRNA, end, sep = ":"))

# Merge and calculate R²
merged_data <- inner_join(rep1, rep2, by = "position")
r_squared <- round(summary(lm(stoich_rep2 ~ stoich_rep1, data = merged_data))$r.squared, 3)

# Plot
plot <- ggplot(merged_data, aes(x = stoich_rep1, y = stoich_rep2)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  coord_fixed() +
  labs(
    title = "Stoichiometry Comparison Between Replicates (BW25113_R1)",
    x = "Replicate 1 Stoichiometry (%)",
    y = "Replicate 2 Stoichiometry (%)"
  ) +
  annotate("text", x = Inf, y = -Inf, label = paste0("R² = ", r_squared),
           hjust = 1.1, vjust = -0.5, size = 4, color = "black") +
  theme_minimal()

# Save plot
ggsave("BW25113_R1_replicability.pdf", plot = plot, width = 6, height = 6)



### Panel C

### Top - pseU, S. cerevisiae

# Read the additional reference file
reference_file <- "Yeast_rRNAmods_track.bed"
reference_data <- read.table(reference_file, header = FALSE)

# Process the reference data
reference_data <- reference_data %>% 
  filter(grepl("Ψ", V4)) %>%  # Keep rows containing 'Ψ'
  mutate(position = paste(V1, V3, sep=":")) %>%  # Merge columns 1 and 3 into 'position'
  select(position, modification = V4)  # Rename column 4 to 'modification'

# File path
file_path <- "BY4741-1_rep2.bedgraph"

# Ensure reference_data$position is character (not factor or int)
reference_data <- reference_data %>%
  mutate(position = as.character(position))

# Read and process the data
df <- read_tsv(file_path, col_names = c("rRNA", "start", "end", "fraction_modified"), col_types = cols()) %>%
  mutate(position = paste0(rRNA, "s:", end)) %>%
  select(position, fraction_modified) %>%
  filter(fraction_modified > 10) %>%
  left_join(reference_data, by = "position") %>%
  mutate(modification = ifelse(is.na(modification), "unmodified", modification))


# Check that at least some Ψ were matched
table(df$modification)

# Plot
plot <- ggplot(df, aes(x = reorder(position, -fraction_modified), y = fraction_modified, fill = modification)) +
  geom_bar(stat = "identity", color = "black", size = 0.2) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("unmodified" = "#4E79A7", "Ψ" = "#E15759")) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(
    title = "Fraction Modified for BY4741-1",
    x = "Position",
    y = "Fraction Modified",
    fill = "Modification"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

### Bottom - pseU, E. coli

# Read the additional reference file
reference_file <- "/no_backup/enovoa/reference_fasta/bacterial_references/rRNA_Ecoli_Mods.bed"
reference_data <- read.table(reference_file, header = FALSE)

# Process the reference data
reference_data <- reference_data %>% 
  filter(grepl("Y", V4)) %>%  # Keep rows containing 'Ψ'
  mutate(position = paste(V1, V3, sep=":")) %>%  # Merge columns 1 and 3 into 'position' 
  mutate(position = gsub("S", "s", position)) %>%  # Replace lowercase 's' with uppercase 'S'
  select(position, modification = V4)  # Rename column 4 to 'modification'

# File path
file_path <- "/no_backup/enovoa/nextflow_outputs/RNA004_mod_benchmarking/e_coli_pseudoU/modkit_analysis/filtered_pileups/BW25113_R1.bedgraph"

# Ensure reference_data$position is character (not factor or int)
reference_data <- reference_data %>%
  mutate(position = as.character(position))

# Read and process the data
df <- read_tsv(file_path, col_names = c("rRNA", "start", "end", "fraction_modified"), col_types = cols()) %>%
  mutate(position = paste0(rRNA, "s:", end)) %>%
  select(position, fraction_modified) %>%
  filter(fraction_modified > 10) %>%
  left_join(reference_data, by = "position") %>%
  mutate(modification = ifelse(is.na(modification), "unmodified", modification))



plot <- ggplot(df, aes(x = reorder(position, -fraction_modified), y = fraction_modified, fill = modification)) +
  geom_bar(stat = "identity", color = "black", size = 0.2) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("unmodified" = "#4E79A7", "Y" = "#E15759")) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(
    title = "Fraction Modified for BW25113",
    x = "Position",
    y = "Fraction Modified",
    fill = "Modification"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )








### Panel D


### Top - E. coli, m5C

# Read the additional reference file
reference_file <- "rRNA_Ecoli_Mods.bed"
reference_data <- read.table(reference_file, header = FALSE)

# Process the reference data
reference_data <- reference_data %>% 
  filter(grepl("m5C", V4)) %>%  # Keep rows containing 'm5C'
  mutate(position = paste(V1, V3, sep=":")) %>%  # Merge columns 1 and 3 into 'position' 
  mutate(position = gsub("S", "s", position)) %>%  # Replace lowercase 's' with uppercase 'S'
  select(position, modification = V4)  # Rename column 4 to 'modification'

# File path
file_path <- "BW25113_R1_rep1.bedgraph"

# Ensure reference_data$position is character (not factor or int)
reference_data <- reference_data %>%
  mutate(position = as.character(position))

# Read and process the data
df <- read_tsv(file_path, col_names = c("rRNA", "start", "end", "fraction_modified"), col_types = cols()) %>%
  mutate(position = paste0(rRNA, "s:", end)) %>%
  select(position, fraction_modified) %>%
  filter(fraction_modified > 10) %>%
  left_join(reference_data, by = "position") %>%
  mutate(modification = ifelse(is.na(modification), "unmodified", modification))



plot <- ggplot(df, aes(x = reorder(position, -fraction_modified), y = fraction_modified, fill = modification)) +
  geom_bar(stat = "identity", color = "black", size = 0.2) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("unmodified" = "#4E79A7", "m5C" = "#E15759")) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(
    title = "Fraction Modified for BY4741-1",
    x = "Position",
    y = "Fraction Modified",
    fill = "Modification"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )


### Bottom, S. cerevisiae, m5C

# Load libraries
library(ggplot2)
library(dplyr)


# Read the additional reference file
reference_file <- "Yeast_rRNAmods_track.bed"
reference_data <- read.table(reference_file, header = FALSE)

# Process the reference data
reference_data <- reference_data %>% 
  filter(grepl("m5C", V4)) %>%  # Keep rows containing 'm5C'
  mutate(position = paste(V1, V3, sep=":")) %>%  # Merge columns 1 and 3 into 'position' 
  mutate(position = gsub("S", "s", position)) %>%  # Replace lowercase 's' with uppercase 'S'
  select(position, modification = V4)  # Rename column 4 to 'modification'

# File path
file_path <- "BY4741-1_rep2.bedgraph"

# Ensure reference_data$position is character (not factor or int)
reference_data <- reference_data %>%
  mutate(position = as.character(position))

# Read and process the data
df <- read_tsv(file_path, col_names = c("rRNA", "start", "end", "fraction_modified"), col_types = cols()) %>%
  mutate(position = paste0(rRNA, "s:", end)) %>%
  select(position, fraction_modified) %>%
  filter(fraction_modified > 10) %>%
  left_join(reference_data, by = "position") %>%
  mutate(modification = ifelse(is.na(modification), "unmodified", modification))



plot <- ggplot(df, aes(x = reorder(position, -fraction_modified), y = fraction_modified, fill = modification)) +
  geom_bar(stat = "identity", color = "black", size = 0.2) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("unmodified" = "#4E79A7", "m5C" = "#E15759")) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(
    title = "Fraction Modified for BY4741",
    x = "Position",
    y = "Fraction Modified",
    fill = "Modification"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )













