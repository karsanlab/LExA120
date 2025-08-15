###############################################################################
# Linear Prediction Scores (LPSs) and classification calls for DLBCL90 tasks
#
# Author         | Dr. Shujun Huang 
# Maintainer     | Joshua Bridgers (jbridgers@bcgsc.ca)
# PI (Lab)       | Dr. Aly Karsan (Karsan Lab)
#                |   https://www.bcgsc.ca/labs/karsan-lab
# 
# Date Created   | 2021-12-12
# Description    | Script to generate DLBCL90 related LPS and classification calls, 
#                | including COO, PMBL, and DHITsig
# Details        
# - Import samples information and their historical calls
# - Data pre-processing:
#   - Read the NanoString data from the RCC files
#   - Construct an expression matrix (rows are for the LExA probes, columns are for the samples, 
#     cells are raw count expression data)
# 
# - The analyses include the following steps:
#   - Step 1: The DLBCL90 r package was applied to the LExA assay expression data for three separate tasks: 
#             - the PMBL LPS scores and calls, 
#             - the COO LPS scores and  calls, 
#             - the DHITsig LPS scores and calls. 
#   - Step 2: Merge the LExA results with known results from the previously separately-run DLBCL90 assay.  
#   - Step 3: Compare the LExA results with the known DLBCL90 results
#   - Step 4: Output the results
#                
# License        | MIT License (see https://opensource.org/licenses/MIT)
#
###############################################################################



# set up --------------------------------------------------------------------------------------

# Load required libraries
library(nanostringr)
library(readxl)
library(DLBCL90)
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




# import historical lps and calls -------------------------------------------------------------

# Read in the historical panel DLBCL90 scores for reference
dlbcl90_calls <- read_csv(here::here("data", "dlbcl90_and_rhl30", 
                                     "historical_scores_of_62_samples_20210322.csv"))
colnames(dlbcl90_calls) <- paste0(colnames(dlbcl90_calls), ".dlbcl90")

dlbcl90_calls_more <- read_csv(here::here("data", "dlbcl90_and_rhl30", 
                                          "historical_scores_of_240_samples_20210513.csv")) %>% 
  dplyr::rename(score = risk_score)
colnames(dlbcl90_calls_more) <- paste0(colnames(dlbcl90_calls_more), ".dlbcl90")




# define analysis configuration of all rounds -------------------------------------------------

rounds_config <- list(
  BCC_Round_1 = list(
    round_name = "BCC_Round1",
    rcc_paths = list(
      here::here("data", "nanostring", "round1_20200417_D48041_to_D48052_RCC"),
      here::here("data", "nanostring", "round1_rerun_20200504_D48049_to_D48050_RCC"),
      here::here("data", "nanostring", "round1_20200501_D48053_to_D48064_RCC")
    ),
    ids_data = IDs,
    dlbcl90_calls = dlbcl90_calls,
    round_col = "LibraryID_Round1",
    sample_naming_function = list(
      function(x) {
        sample_ids <- c("D48041", "D48042", "D48043", "D48044", "D48045", "D48046", 
                        "D48047", "D48048", "D48049", "D48050", "D48051", "D48052")
        # Exclude D48049 and D48050 from first batch as they're in second batch
        sample_ids[!sample_ids %in% c("D48049", "D48050")]
      },
      function(x) c("D48049", "D48050"),
      function(x) c("D48053", "D48054", "D48055", "D48056")
    ),
    prefix_removal = NULL,
    geomean_cut = NULL,
    exclude_samples = NULL,
    id_transformation = NULL
  ),
  
  BCC_Round_2 = list(
    round_name = "BCC_Round2",
    rcc_paths = list(
      here::here("data", "nanostring", "round2_20200513_D70120_to_D70131_RCC"),
      here::here("data", "nanostring", "round2_20200514_D70132_to_D70143_RCC")
    ),
    ids_data = IDs,
    dlbcl90_calls = dlbcl90_calls,
    round_col = "LibraryID_Round2",
    sample_naming_function = list(
      function(x) sapply(x, function(xx) unlist(strsplit(xx, "\\."))[3]),
      function(x) sapply(x, function(xx) unlist(strsplit(xx, "\\."))[2])
    ),
    prefix_removal = NULL,
    geomean_cut = NULL,
    exclude_samples = NULL,
    id_transformation = NULL
  ),
  
  BCC_Round_3 = list(
    round_name = "BCC_Round3",
    rcc_paths = list(
      here::here("data", "nanostring", "round3_20200630_AK1110_to_AK1121_RCC"),
      here::here("data", "nanostring", "round3_20200703_AK1122_to_AK1133_RCC")
    ),
    ids_data = IDs,
    dlbcl90_calls = dlbcl90_calls,
    round_col = "LibraryID_Round3",
    sample_naming_function = NULL,
    prefix_removal = "AK.",
    geomean_cut = NULL,
    exclude_samples = NULL,
    id_transformation = NULL
  ),
  
  BCC_Round_4 = list(
    round_name = "BCC_Round4",
    rcc_paths = list(
      here::here("data", "nanostring", "round4_batch1"),
      here::here("data", "nanostring", "round4_batch2"),
      here::here("data", "nanostring", "round4_batch3"),
      here::here("data", "nanostring", "round4_batch4")
    ),
    ids_data = IDs,
    dlbcl90_calls = dlbcl90_calls,
    round_col = "LibraryID_Round4",
    sample_naming_function = list(
      NULL,  # Use default for batch 1
      NULL,  # Use default for batch 2
      function(x) {
        sample_ids <- c("F04999", "F10327", "F10328", "F10329", "F10330", "F10331")
        sample_ids[1]  # Only first sample for batch 3
      },
      NULL   # Use default for batch 4
    ),
    prefix_removal = "AK.",
    geomean_cut = NULL,
    exclude_samples = NULL,
    id_transformation = NULL
  ),
  
  UHN_Round_1 = list(
    round_name = "UHN_Round1",
    rcc_paths = list(
      here::here("data", "nanostring", "UHN_round1", "UHN_20200717"),
      here::here("data", "nanostring", "UHN_round1", "UHN_20200724")
    ),
    ids_data = IDs,
    dlbcl90_calls = dlbcl90_calls,
    round_col = "UHN_Round1",
    sample_naming_function = list(
      function(x) paste0(rep("UHN_066_01_"), sprintf("%02d", seq(1:12))),
      function(x) paste0(rep("UHN_070_01_"), sprintf("%02d", seq(1:12)))
    ),
    prefix_removal = NULL,
    geomean_cut = NULL,
    exclude_samples = NULL,
    id_transformation = NULL
  ),
  
  UHN_Round_2 = list(
    round_name = "UHN_Round2",
    rcc_paths = list(here::here("data", "nanostring", "UHN_rounds2-4", "round2")),
    ids_data = IDs,
    dlbcl90_calls = dlbcl90_calls,
    round_col = "UHN_Round2",
    sample_naming_function = function(path) {
      list.files(path = path, pattern = ".RCC")
    },
    prefix_removal = NULL,
    geomean_cut = NULL,
    exclude_samples = NULL,
    id_transformation = function(x) gsub("-", ".", x)
  ),
  
  UHN_Round_3 = list(
    round_name = "UHN_Round3",
    rcc_paths = list(here::here("data", "nanostring", "UHN_rounds2-4", "round3")),
    ids_data = IDs,
    dlbcl90_calls = dlbcl90_calls,
    round_col = "UHN_Round3",
    sample_naming_function = function(path) {
      list.files(path = path, pattern = ".RCC")
    },
    prefix_removal = NULL,
    geomean_cut = NULL,
    exclude_samples = NULL,
    id_transformation = function(x) gsub("-", ".", x)
  ),
  
  UHN_Round_4 = list(
    round_name = "UHN_Round4",
    rcc_paths = list(here::here("data", "nanostring", "UHN_rounds2-4", "round4")),
    ids_data = IDs,
    dlbcl90_calls = dlbcl90_calls,
    round_col = "UHN_Round4",
    sample_naming_function = function(path) {
      list.files(path = path, pattern = ".RCC")
    },
    prefix_removal = NULL,
    geomean_cut = 1,
    exclude_samples = NULL,
    id_transformation = function(x) gsub("-", ".", x)
  ),
  
  BCC_Round_5 = list(
    round_name = "BCC_Round5",
    rcc_paths = list(here::here("data", "nanostring", "additional_DLBCL_runs")),
    ids_data = IDs_more,
    dlbcl90_calls = dlbcl90_calls_more,
    round_col = "Additional_Round",
    sample_naming_function = NULL,
    prefix_removal = "AK.",
    geomean_cut = 1,
    exclude_samples = NULL,
    id_transformation = NULL
  )
)



# run all rounds ------------------------------------------------------------------------------

all_results <- list()
all_data <- list()

for (round_id in names(rounds_config)) {

  # Get configuration for current round
  config <- rounds_config[[round_id]]
  
  # Extract configuration parameters
  round_name <- config$round_name
  rcc_paths <- config$rcc_paths
  ids_data <- config$ids_data
  dlbcl90_calls <- config$dlbcl90_calls
  round_col <- config$round_col
  sample_naming_function <- config$sample_naming_function
  prefix_removal <- config$prefix_removal
  geomean_cut <- config$geomean_cut
  exclude_samples <- config$exclude_samples
  id_transformation <- config$id_transformation
  
  # Run analysis with error handling
  tryCatch({
    
    results <- run_dlbcl90_analysis(
      round_name = round_name,
      rcc_paths = rcc_paths,
      ids_data = ids_data,
      dlbcl90_calls = dlbcl90_calls,
      round_col = round_col,
      sample_naming_function = sample_naming_function,
      prefix_removal = prefix_removal,
      geomean_cut = geomean_cut,
      exclude_samples = exclude_samples,
      id_transformation = id_transformation
    )
    
    # Store results
    all_results[[round_id]] <- results$results
    all_data[[round_id]] <- results$data
    
  }, error = function(e) {
    all_results[[round_id]] <- NULL
    all_data[[round_id]] <- NULL
  })
}

