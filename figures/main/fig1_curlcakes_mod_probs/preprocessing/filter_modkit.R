###################################################
## Processing modkit output and annotating sites ##
###################################################

# Task 0: Load packages


pkgs <- c("openintro","tidyverse","backports","here","skimr","dplyr", "ggplot2", "ggsci","ggforce",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr","scales")

lapply(pkgs, library, character.only = TRUE)

#######################
# Task 1: Import Data #
#######################

# Task 1.1: Import ref

ref <- read.csv(here("../__ref/epinano_output_plusstrand.csv"), header = T, sep = "\t")

# Task 1.2: Import modkit output

dir_hac_pool1 <- setwd(here("..","__data","curlcakes","raw","hac_0.1_filter-pass", "Pool1"))
dir_hac_pool3 <- setwd(here("..","__data","curlcakes","raw","hac_0.1_filter-pass", "Pool3"))


file_list_hac_pool1 <- list.files(path = dir_hac_pool1)
file_list_hac_pool3 <- list.files(path = dir_hac_pool3)

# Task 1.3: Function for reading in all files from a directory 

# Create a loop to read in every file of the directory and append it to the initialized data.frame plus add a new column that contains the name of the pool
# Comment: Could be sped up with fread and data.tables

read_dir <- function(file_list, work_dir){
  
  setwd(work_dir)
  
  dataset <- data.frame()
  
  for (i in 1:length(file_list)){
    temp_data <- read_tsv(file_list[i]) #each file will be read in, specify which columns you need read in to avoid any errors # specifying col_types is essential to see spike_ins
    temp_data$sample <-  gsub("\\..*\\.union\\.bed", "", file_list[i])#clean the data as needed, in this case I am creating a new column that indicates which file each row of data came from
    dataset <- plyr::rbind.fill(dataset, temp_data) #for each iteration, bind the new data to the building dataset
  }
  
  rm(i)
  rm(temp_data)
  
  return(dataset)
}

hac_raw_pool1 <- read_dir(file_list = file_list_hac_pool1, work_dir = dir_hac_pool1)
hac_raw_pool3 <- read_dir(file_list = file_list_hac_pool3, work_dir = dir_hac_pool3)



####################################
# Task 2: Prefiltering of raw data #
####################################

filter_df <- function(df) {
  
  df_filtered <- df %>%
    # Remove barcodes not present in the actual run
    dplyr::select(-matches(".*pod5---bc_(5[2-9]|[6-8][0-9]|9[0-6]).*")) %>% 
    # Remove failed m5U barcodes from this run
    dplyr::select(-matches(".*pod5---bc_(16|17|18).*")) %>% 
    # Filter out ranges > 1nt position (introduced by beduniong)
    filter(end - start == 1) %>% 
    # Split filename column into appropriate columns
    separate(col = sample, into = c("sample","accuracy_model"), sep = "\\+") %>%
    separate(col = accuracy_model, into = c("accuracy","model"), sep = "@")
  
  return(df_filtered)
  
}


hac_pool1_filtered <- filter_df(hac_raw_pool1)
hac_pool3_filtered <- filter_df(hac_raw_pool3)


rm(hac_raw_pool1)
rm(hac_raw_pool3)

##########################
# Task 3: Data wrangling #
##########################

# Task 3.1: Function to create a column to overlap on

create_overlap <- function(df, col1, col2) {
  
  df_overlap <- df %>% 
  mutate(overlap=paste({{col1}}, {{col2}}, sep = "_"))
  
}


hac_pool1_filtered <- create_overlap(df=hac_pool1_filtered, col1=chrom, col2=end)
hac_pool3_filtered <- create_overlap(df=hac_pool3_filtered, col1=chrom, col2=end)
ref <- create_overlap(df=ref, col1=X.Ref, col2=pos)




overlap_merge <- function(df, df_ref, model){
  
  stoich <- df %>% 
    filter(model == !!model) %>% 
    full_join(df_ref) %>% 
    fill(sample,accuracy,model, .direction = "updown") %>% 
    dplyr::select(X.Ref, pos, base, strand, model, accuracy, sample, starts_with("percent_mod"), overlap) %>% 
    pivot_longer(names_to = "modifications", values_to = "stoichiometry", cols = c(starts_with("percent_mod"))) %>%
    # Change barcode names to sample names
    mutate(modifications = case_when(
      str_detect(modifications, "percent_mod_pod5---bc_1\\b") ~ "UNM_rep1",
      str_detect(modifications, "percent_mod_pod5---bc_2\\b") ~ "UNM_rep2",
      str_detect(modifications, "percent_mod_pod5---bc_3\\b") ~ "UNM_rep3",
      str_detect(modifications, "percent_mod_pod5---bc_4\\b") ~ "ac4C_rep1",
      str_detect(modifications, "percent_mod_pod5---bc_5\\b") ~ "ac4C_rep2",
      str_detect(modifications, "percent_mod_pod5---bc_6\\b") ~ "ac4C_rep3",
      str_detect(modifications, "percent_mod_pod5---bc_7\\b") ~ "hm5C_rep1",
      str_detect(modifications, "percent_mod_pod5---bc_8\\b") ~ "hm5C_rep2",
      str_detect(modifications, "percent_mod_pod5---bc_9\\b") ~ "hm5C_rep3",
      str_detect(modifications, "percent_mod_pod5---bc_10\\b") ~ "m1Y_rep1",
      str_detect(modifications, "percent_mod_pod5---bc_11\\b") ~ "m1Y_rep2",
      str_detect(modifications, "percent_mod_pod5---bc_12\\b") ~ "m1Y_rep3",
      str_detect(modifications, "percent_mod_pod5---bc_13\\b") ~ "m5C_rep1",
      str_detect(modifications, "percent_mod_pod5---bc_14\\b") ~ "m5C_rep2",
      str_detect(modifications, "percent_mod_pod5---bc_15\\b") ~ "m5C_rep3",
      
      str_detect(modifications, "percent_mod_pod5---bc_19\\b") ~ "m6A_12.5_rep1",
      str_detect(modifications, "percent_mod_pod5---bc_20\\b") ~ "m6A_12.5_rep2",
      str_detect(modifications, "percent_mod_pod5---bc_21\\b") ~ "m6A_12.5_rep3",
      str_detect(modifications, "percent_mod_pod5---bc_22\\b") ~ "m6A_25_rep1",
      str_detect(modifications, "percent_mod_pod5---bc_23\\b") ~ "m6A_25_rep2",
      str_detect(modifications, "percent_mod_pod5---bc_24\\b") ~ "m6A_25_rep3",
      str_detect(modifications, "percent_mod_pod5---bc_25\\b") ~ "m6A_50_rep1",
      str_detect(modifications, "percent_mod_pod5---bc_26\\b") ~ "m6A_50_rep2",
      str_detect(modifications, "percent_mod_pod5---bc_27\\b") ~ "m6A_50_rep3",
      str_detect(modifications, "percent_mod_pod5---bc_28\\b") ~ "m6A_75_rep1",
      str_detect(modifications, "percent_mod_pod5---bc_29\\b") ~ "m6A_75_rep2",
      str_detect(modifications, "percent_mod_pod5---bc_30\\b") ~ "m6A_75_rep3",
      str_detect(modifications, "percent_mod_pod5---bc_31\\b") ~ "m6A_100_rep1",
      str_detect(modifications, "percent_mod_pod5---bc_32\\b") ~ "m6A_100_rep2",
      str_detect(modifications, "percent_mod_pod5---bc_33\\b") ~ "m6A_100_rep3",
      str_detect(modifications, "percent_mod_pod5---bc_34\\b") ~ "Y_12.5_rep1",
      str_detect(modifications, "percent_mod_pod5---bc_35\\b") ~ "Y_12.5_rep2",
      str_detect(modifications, "percent_mod_pod5---bc_36\\b") ~ "Y_12.5_rep3",
      str_detect(modifications, "percent_mod_pod5---bc_37\\b") ~ "Y_25_rep1",
      str_detect(modifications, "percent_mod_pod5---bc_38\\b") ~ "Y_25_rep2",
      str_detect(modifications, "percent_mod_pod5---bc_39\\b") ~ "Y_25_rep3",
      str_detect(modifications, "percent_mod_pod5---bc_40\\b") ~ "Y_50_rep1",
      str_detect(modifications, "percent_mod_pod5---bc_41\\b") ~ "Y_50_rep2",
      str_detect(modifications, "percent_mod_pod5---bc_42\\b") ~ "Y_50_rep3",
      str_detect(modifications, "percent_mod_pod5---bc_43\\b") ~ "Y_75_rep1",
      str_detect(modifications, "percent_mod_pod5---bc_44\\b") ~ "Y_75_rep2",
      str_detect(modifications, "percent_mod_pod5---bc_45\\b") ~ "Y_75_rep3",
      str_detect(modifications, "percent_mod_pod5---bc_46\\b") ~ "Y_100_rep1",
      str_detect(modifications, "percent_mod_pod5---bc_47\\b") ~ "Y_100_rep2",
      str_detect(modifications, "percent_mod_pod5---bc_48\\b") ~ "Y_100_rep3",
      str_detect(modifications, "percent_mod_pod5---bc_49\\b") ~ "m5U_100_rep1",
      str_detect(modifications, "percent_mod_pod5---bc_50\\b") ~ "m5U_100_rep2",
      str_detect(modifications, "percent_mod_pod5---bc_51\\b") ~ "m5U_100_rep3",
      
      # Barcodes not present in run RNA241003_Pool1 are filtered out
      FALSE ~ NA
    )) %>%
    #separate(col = modifications, into = c("modification","replicate"), sep = "_(?=[^_]*$)", remove = T) %>% 
    unite(overlap_2, c("overlap","modifications"), remove = F)
  
  
  
  cov <- df %>%  
    filter(model == !!model) %>% 
    full_join(df_ref) %>% 
    fill(sample,accuracy,model, .direction = "updown") %>% 
    ### Get coverage in long format
    dplyr::select(starts_with("cov"), overlap) %>% 
    pivot_longer(names_to = "modifications_2", values_to = "coverage", cols = c(starts_with("cov"))) %>%
    # Change barcode names to sample names
    mutate(modifications_2 = case_when(
      str_detect(modifications_2, "cov_pod5---bc_1\\b") ~ "UNM_rep1",
      str_detect(modifications_2, "cov_pod5---bc_2\\b") ~ "UNM_rep2",
      str_detect(modifications_2, "cov_pod5---bc_3\\b") ~ "UNM_rep3",
      str_detect(modifications_2, "cov_pod5---bc_4\\b") ~ "ac4C_rep1",
      str_detect(modifications_2, "cov_pod5---bc_5\\b") ~ "ac4C_rep2",
      str_detect(modifications_2, "cov_pod5---bc_6\\b") ~ "ac4C_rep3",
      str_detect(modifications_2, "cov_pod5---bc_7\\b") ~ "hm5C_rep1",
      str_detect(modifications_2, "cov_pod5---bc_8\\b") ~ "hm5C_rep2",
      str_detect(modifications_2, "cov_pod5---bc_9\\b") ~ "hm5C_rep3",
      str_detect(modifications_2, "cov_pod5---bc_10\\b") ~ "m1Y_rep1",
      str_detect(modifications_2, "cov_pod5---bc_11\\b") ~ "m1Y_rep2",
      str_detect(modifications_2, "cov_pod5---bc_12\\b") ~ "m1Y_rep3",
      str_detect(modifications_2, "cov_pod5---bc_13\\b") ~ "m5C_rep1",
      str_detect(modifications_2, "cov_pod5---bc_14\\b") ~ "m5C_rep2",
      str_detect(modifications_2, "cov_pod5---bc_15\\b") ~ "m5C_rep3",
      
      str_detect(modifications_2, "cov_pod5---bc_19\\b") ~ "m6A_12.5_rep1",
      str_detect(modifications_2, "cov_pod5---bc_20\\b") ~ "m6A_12.5_rep2",
      str_detect(modifications_2, "cov_pod5---bc_21\\b") ~ "m6A_12.5_rep3",
      str_detect(modifications_2, "cov_pod5---bc_22\\b") ~ "m6A_25_rep1",
      str_detect(modifications_2, "cov_pod5---bc_23\\b") ~ "m6A_25_rep2",
      str_detect(modifications_2, "cov_pod5---bc_24\\b") ~ "m6A_25_rep3",
      str_detect(modifications_2, "cov_pod5---bc_25\\b") ~ "m6A_50_rep1",
      str_detect(modifications_2, "cov_pod5---bc_26\\b") ~ "m6A_50_rep2",
      str_detect(modifications_2, "cov_pod5---bc_27\\b") ~ "m6A_50_rep3",
      str_detect(modifications_2, "cov_pod5---bc_28\\b") ~ "m6A_75_rep1",
      str_detect(modifications_2, "cov_pod5---bc_29\\b") ~ "m6A_75_rep2",
      str_detect(modifications_2, "cov_pod5---bc_30\\b") ~ "m6A_75_rep3",
      str_detect(modifications_2, "cov_pod5---bc_31\\b") ~ "m6A_100_rep1",
      str_detect(modifications_2, "cov_pod5---bc_32\\b") ~ "m6A_100_rep2",
      str_detect(modifications_2, "cov_pod5---bc_33\\b") ~ "m6A_100_rep3",
      str_detect(modifications_2, "cov_pod5---bc_34\\b") ~ "Y_12.5_rep1",
      str_detect(modifications_2, "cov_pod5---bc_35\\b") ~ "Y_12.5_rep2",
      str_detect(modifications_2, "cov_pod5---bc_36\\b") ~ "Y_12.5_rep3",
      str_detect(modifications_2, "cov_pod5---bc_37\\b") ~ "Y_25_rep1",
      str_detect(modifications_2, "cov_pod5---bc_38\\b") ~ "Y_25_rep2",
      str_detect(modifications_2, "cov_pod5---bc_39\\b") ~ "Y_25_rep3",
      str_detect(modifications_2, "cov_pod5---bc_40\\b") ~ "Y_50_rep1",
      str_detect(modifications_2, "cov_pod5---bc_41\\b") ~ "Y_50_rep2",
      str_detect(modifications_2, "cov_pod5---bc_42\\b") ~ "Y_50_rep3",
      str_detect(modifications_2, "cov_pod5---bc_43\\b") ~ "Y_75_rep1",
      str_detect(modifications_2, "cov_pod5---bc_44\\b") ~ "Y_75_rep2",
      str_detect(modifications_2, "cov_pod5---bc_45\\b") ~ "Y_75_rep3",
      str_detect(modifications_2, "cov_pod5---bc_46\\b") ~ "Y_100_rep1",
      str_detect(modifications_2, "cov_pod5---bc_47\\b") ~ "Y_100_rep2",
      str_detect(modifications_2, "cov_pod5---bc_48\\b") ~ "Y_100_rep3",
      str_detect(modifications_2, "cov_pod5---bc_49\\b") ~ "m5U_100_rep1",
      str_detect(modifications_2, "cov_pod5---bc_50\\b") ~ "m5U_100_rep2",
      str_detect(modifications_2, "cov_pod5---bc_51\\b") ~ "m5U_100_rep3",
      
      # Barcodes not present in run RNA241003_Pool1 are filtered out
      FALSE ~ NA
    )) %>% 
    unite(overlap_2, c("overlap","modifications_2"), remove = F) %>% 
    dplyr::select(overlap_2, coverage)
  
  
  final_df <- full_join(stoich, cov ,by = join_by(overlap_2)) %>% 
    separate(col = modifications, into = c("modification","replicate"), sep = "_(?=[^_]*$)", remove = T) %>% 
    select(-overlap_2)
  
  return(final_df)
  
}


pool1_pseU <- overlap_merge(df = hac_pool1_filtered, df_ref = ref, model = "pseU") %>%
  # Remove m5U_100 from pool1
  filter(!((modification == "m5U_100")))

pool1_m6A_DRACH <- overlap_merge(df = hac_pool1_filtered, df_ref = ref, model = "m6A_DRACH") %>%
  # Remove m5U_100 from pool1
  filter(!((modification == "m5U_100")))

pool1_inosine_m6A_17596 <- overlap_merge(df = hac_pool1_filtered, df_ref = ref, model = "inosine_m6A-17596") %>%
  # Remove m5U_100 from pool1
  filter(!((modification == "m5U_100")))

pool1_inosine_m6A_a <- overlap_merge(df = hac_pool1_filtered, df_ref = ref, model = "inosine_m6A-a") %>%
  # Remove m5U_100 from pool1
  filter(!((modification == "m5U_100")))

pool1_m5C <- overlap_merge(df = hac_pool1_filtered, df_ref = ref, model = "m5C") %>%
  # Remove m5U_100 from pool1
  filter(!((modification == "m5U_100")))



pool3_pseU <- overlap_merge(df = hac_pool3_filtered, df_ref = ref, model = "pseU") %>% 
  # Keep only m5U_100 from pool3
  filter(modification == "m5U_100")

pool3_m6A_DRACH <- overlap_merge(df = hac_pool3_filtered, df_ref = ref, model = "m6A_DRACH") %>% 
  # Keep only m5U_100 from pool3
  filter(modification == "m5U_100")

pool3_inosine_m6A_17596 <- overlap_merge(df = hac_pool3_filtered, df_ref = ref, model = "inosine_m6A-17596") %>% 
  # Keep only m5U_100 from pool3
  filter(modification == "m5U_100")

pool3_inosine_m6A_a <- overlap_merge(df = hac_pool3_filtered, df_ref = ref, model = "inosine_m6A-a") %>% 
  # Keep only m5U_100 from pool3
  filter(modification == "m5U_100")

pool3_m5C <- overlap_merge(df = hac_pool3_filtered, df_ref = ref, model = "m5C") %>% 
  # Keep only m5U_100 from pool3
  filter(modification == "m5U_100")


master_df <- rbind(pool1_pseU, pool1_m6A_DRACH, pool1_m5C, pool1_inosine_m6A_a, pool1_inosine_m6A_17596,
                   pool3_pseU, pool3_m6A_DRACH, pool3_m5C, pool3_inosine_m6A_a, pool3_inosine_m6A_17596)

master_df <- rbind(pool1_pseU, pool1_m6A_DRACH, pool1_m5C, pool1_inosine_m6A_a, pool1_inosine_m6A_17596,
                   pool3_pseU)

write.table(master_df, file=here("..","__data","curlcakes","processed","hac_0.1_filter-pass", "master_table.tsv"), sep="\t", quote=FALSE,row.names=FALSE)





