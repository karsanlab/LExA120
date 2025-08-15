###############################################################################
# Linear Prediction Scores (LPSs) and classification calls for RHL30 task
#
# Author         | Dr. Shujun Huang 
# Maintainer     | Joshua Bridgers (jbridgers@bcgsc.ca)
# PI (Lab)       | Dr. Aly Karsan (Karsan Lab)
#                |   https://www.bcgsc.ca/labs/karsan-lab
# 
# Date Created   | 2021-12-12
# Description    | Script to generate LExA cHL risk classification LPS for all rounds
# Details        
# - Import samples information and their historical calls
# - Data pre-processing:
#   - Read the NanoString data from the RCC files
#   - Construct an expression matrix (rows are for the LExA probes, columns are for the samples, 
#     cells are raw count expression data)
# 
# - The analyses include the following steps:
#   - Step 1: calculate the RHL30 model risk scores using RHL30 package and determine risk groups 
#             associate with Post-Autologous Stem Cell Transplantation (post-ASCT) outcomes for classical 
#             Hodgkin Lymphoma (cHL) samples from the LExA NanoString expression counts;
#   - Step 2: merge the LExA calls with the previous RHL30 calls;
#   - Step 3: compare the LExA results with the known RHL30 results
#   - Step 4: Output the results for further visualization
#                
# License        | MIT License (see https://opensource.org/licenses/MIT)
#
###############################################################################



# set up --------------------------------------------------------------------------------------

# Load required libraries
library(nanostringr)
library(readxl)
library(RHL30)
library(caret)
library(e1071)
library(ggrepel)
library(tidyverse)
library(here)

# Source defined helpers
source("src/utils.R")






# import sample library ids -------------------------------------------------------------------

# Read in library IDs of samples in different rounds 
IDs <- read.csv(here::here("data", "library_IDs_of_62_samples_20210322.csv"), 
                header = TRUE, stringsAsFactors = FALSE)

IDs_more <- read.csv(here::here("data", "library_IDs_of_240_samples_20210513.csv"), 
                     header = TRUE, stringsAsFactors = FALSE)

IDs_uhn <- read.csv(here::here("data", "nanostring", "RHL30_LexA_comparison_samples", "UHN_HL_sample_IDs.csv"), 
                    header = TRUE, stringsAsFactors = FALSE) %>% 
  separate(Sample_ID, into = c("AMDL_ID", "PMH_ID"), sep = "_")




# import historical lps and calls -------------------------------------------------------------

# Read in the historical panel RHL30 scores for reference
rhl30_calls <- read_csv(here::here("data", "dlbcl90_and_rhl30", "historical_scores_of_62_samples_20210322.csv")) %>% 
  dplyr::rename(score = risk_score)
colnames(rhl30_calls) <- paste0(colnames(rhl30_calls), ".rhl30")

rhl30_calls_more <- read_csv(here::here("data", "dlbcl90_and_rhl30", "historical_scores_of_240_samples_20210513.csv")) %>% 
  dplyr::rename(score = risk_score)
colnames(rhl30_calls_more) <- paste0(colnames(rhl30_calls_more), ".rhl30")

rhl30_calls_uhn <- read_excel(here::here("data", "dlbcl90_and_rhl30", "historical_scores_of_UHN_additional_cHL_samples.xlsx"))
colnames(rhl30_calls_uhn) <- paste0(colnames(rhl30_calls_uhn), ".rhl30")




# retrieve RHL30 model info -------------------------------------------------------------------

# RHL30 model coefficients and genes from RHL30 package
rhl30_model_df <- get_rhl30_model_coef_df()
hk_genes <- filter(rhl30_model_df, gene_type == "housekeeper") %>% pull("gene_name")




# define analysis configuration of all rounds -------------------------------------------------

round_configs <- list(
  
  # BCC Rounds - using existing helper functions
  BCC_Round_1 = list(
    round_name = "BCC_Round1",
    rcc_paths = here::here("data", "nanostring", "round1_20200501_D48053_to_D48064_RCC"),
    ids_data = IDs,
    reference_calls = rhl30_calls,
    round_col = "LibraryID_Round1",
    sample_naming_function = function(x) sapply(x, function(y) unlist(strsplit(y, "\\."))[3]),
    exclude_samples = "AK1133"
  ),
  
  BCC_Round_2 = list(
    round_name = "BCC_Round2",
    rcc_paths = here::here("data", "nanostring", "round2_20200514_D70132_to_D70143_RCC"),
    ids_data = IDs,
    reference_calls = rhl30_calls,
    round_col = "LibraryID_Round2",
    sample_naming_function = function(x) sapply(x, function(y) unlist(strsplit(y, "\\."))[2]),
    exclude_samples = "AK1133"
  ),
  
  BCC_Round_3 = list(
    round_name = "BCC_Round3",
    rcc_paths = here::here("data", "nanostring", "round3_20200703_AK1122_to_AK1133_RCC"),
    ids_data = IDs,
    reference_calls = rhl30_calls,
    round_col = "LibraryID_Round3",
    exclude_samples = "AK1133"
  ),
  
  BCC_Round_4 = list(
    round_name = "BCC_Round4",
    rcc_paths = here::here("data", "nanostring", "round4_batch3"),
    ids_data = IDs,
    reference_calls = rhl30_calls,
    round_col = "LibraryID_Round4",
    exclude_samples = "AK1133"
  ),
  
  # UHN Rounds
  UHN_Round_1 = list(
    round_name = "UHN_Round1",
    rcc_paths = list(
      here::here("data", "nanostring", "UHN_round1", "UHN_20200717"),
      here::here("data", "nanostring", "UHN_round1", "UHN_20200724")
    ),
    ids_data = IDs,
    reference_calls = rhl30_calls,
    round_col = "UHN_Round1",
    sample_naming_function = list(
      function(x) paste0(rep("UHN_066_01_"), sprintf("%02d", seq(1:12))),
      function(x) paste0(rep("UHN_070_01_"), sprintf("%02d", seq(1:12)))
    ),
    exclude_samples = "AK1133"
  ),
  
  UHN_Round_2 = list(
    round_name = "UHN_Round2",
    rcc_paths = here::here("data", "nanostring", "UHN_rounds2-4", "round2"),
    ids_data = IDs,
    reference_calls = rhl30_calls,
    round_col = "UHN_Round2",
    sample_naming_function = function(x) {
      list.files(path = here::here("data", "nanostring", "UHN_rounds2-4", "round2"), pattern = ".RCC")
    },
    exclude_samples = "AK1133"
  ),
  
  UHN_Round_3 = list(
    round_name = "UHN_Round3",
    rcc_paths = here::here("data", "nanostring", "UHN_rounds2-4", "round3"),
    ids_data = IDs,
    reference_calls = rhl30_calls,
    round_col = "UHN_Round3",
    sample_naming_function = function(x) {
      list.files(path = here::here("data", "nanostring", "UHN_rounds2-4", "round3"), pattern = ".RCC")
    },
    exclude_samples = "AK1133"
  ),
  
  UHN_Round_4 = list(
    round_name = "UHN_Round4",
    rcc_paths = here::here("data", "nanostring", "UHN_rounds2-4", "round4"),
    ids_data = IDs,
    reference_calls = rhl30_calls,
    round_col = "UHN_Round4",
    sample_naming_function = function(x) {
      list.files(path = here::here("data", "nanostring", "UHN_rounds2-4", "round4"), pattern = ".RCC")
    },
    exclude_samples = "AK1133"
  ),
  
  # Special cases
  UHN_additional_cHLs = list(
    round_name = "UHN_Round5",
    rcc_paths = here::here("data", "nanostring", "RHL30_LexA_comparison_samples"),
    ids_data = IDs_uhn,
    reference_calls = rhl30_calls_uhn,
    round_col = "RCC_file_ID",
    sample_naming_function = function(x) {
      list.files(path = here::here("data", "nanostring", "RHL30_LexA_comparison_samples"), pattern = ".RCC")
    },
    special_processing = "uhn_additional"
  ),
  
  BCC_additional_cHLs = list(
    round_name = "BCC_Round5",
    rcc_paths = here::here("data", "nanostring", "additional_cHL_runs"),
    ids_data = IDs_more,
    reference_calls = rhl30_calls_more,
    round_col = "Additional_Round",
    sample_naming_function = function(x) sapply(x, function(y) unlist(strsplit(y, "\\."))[2]),
    special_processing = "bcc_additional",
    exclude_samples = c("F48581", "F48595", "F48596", "F48597")
  )
)




# run all rounds ------------------------------------------------------------------------------

all_results <- list()
all_data <- list()
summary_table <- data.frame()

for (round_name in names(round_configs)) {
  
  config <- round_configs[[round_name]]
  
  tryCatch({
    # Process the round using existing helper functions
    result <- process_single_round(config, hk_genes, rhl30_model_df)
    
    # Store results
    all_results[[round_name]] <- result$results
    all_data[[round_name]] <- result$data
    
    # Add to summary table
    summary_row <- data.frame(
      Round = round_name,
      Size = result$results$res_risk["Size"],
      Pearson = result$results$res_risk["Pearson"],
      R_squared = result$results$res_risk["R_squared"],
      Misclassifications = result$results$res_risk["Misclassification"],
      stringsAsFactors = FALSE
    )
    
    summary_table <- rbind(summary_table, summary_row)
    
  }, error = function(e) {
    all_results[[round_name]] <- NULL
    all_data[[round_name]] <- NULL
  })
}


