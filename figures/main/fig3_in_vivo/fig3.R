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

# Load libraries
library(ggplot2)
library(dplyr)
library(readr)

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
