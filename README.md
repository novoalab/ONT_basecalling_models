# Benchmarking of modification-aware basecalling models

![](./img/Fig4B_snapshot.png "Figure 4B")

## Project Overview  
This repository contains the code required to reproduce all figures that are part of our manuscript on benchmarking modification-aware basecalling models for direct RNA-sequencing. The following [dorado](https://github.com/nanoporetech/dorado) models were tested:

| Basecalling Models | Compatible<br />Modifications | Modifications<br />Model<br />Version | Data<br />Sampling<br />Frequency |
| :-------- | :------- | :--- | :--- |
| **rna004_130bps_hac@v5.1.0** | m5C<br />m6A_DRACH<br />inosine_m6A<br />pseU | v1<br />v1<br />v1<br />v1 | 4 kHz |
| **rna004_130bps_sup@v5.1.0** | m5C<br />m6A_DRACH<br />inosine_m6A<br />pseU | v1<br />v1<br />v1<br />v1 | 4 kHz |

## Installation

###  **Clone this repository** 

```bash
git clone https://github.com/YOUR-USERNAME/ONT_basecalling_models.git

cd ONT_basecalling_models
```

###  **Set up dependencies** 

Ensure R (> 4.3.2) and required packages are installed.

Install key dependencies manually:

| Package | Version |
| :-------- | :------- |
| ggdist |	3.3.2 |
| rstatix |	0.7.2 |
|scales |	1.3.0 |
|ggpubr |	0.6.0 |
|usethis |	2.2.3 |
|ggrepel |	0.9.5 |
|MetBrewer |	0.2.0 |
|xlsx |	0.6.5 |
|readxl |	1.4.3 |
|janitor |	2.2.0 |
|ggforce |	0.4.2 |
|ggsci |	3.0.2 |
|skimr |	2.1.5 |
|here | 1.0.1 |
|backports |	1.4.1 |
|lubridate |	1.9.3 |
|forcats |	1.0.0 |
|stringr |	1.5.1 |
|dplyr |	1.1.4 |
|purrr |	1.0.4 |
|readr |	2.1.5 |
|tidyr |	1.3.1 |
|tibble |	3.2.1 |
|ggplot2 |	3.5.1 |
|tidyverse |	2.0.0 |
|openintro |	2.4.0 |

## Citation
If you find this work/code useful, please cite our manuscript: 

Diensthuber G, Milenkovic I, Llovera L, Milovanovic A, Pelizzari F and Novoa EM. **Systematic benchmarking of basecalling models for RNA modification detection with highly-multiplexed nanopore sequencing** [bioRxiv, 2025](https://www.biorxiv.org/content/10.1101/2025.07.11.663424v1). 

