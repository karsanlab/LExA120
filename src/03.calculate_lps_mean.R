###############################################################################
#
# Author         | Dr. Shujun Huang
# Maintainer     | Joshua Bridgers (jbridgers@bcgsc.ca)
# PI (Lab)       | Dr. Aly Karsan (Karsan Lab)
#                |   https://www.bcgsc.ca/labs/karsan-lab
#
# Date Created   | 2021-12-12
# Description    | calculate mean lps of Rounds 1-3 for BCC and UHN
#
# License        | MIT License (see https://opensource.org/licenses/MIT)
#
###############################################################################


# set up --------------------------------------------------------------------------------------

# Load required libraries
library(here)
library(tidyverse)

# Source defined helpers
source("src/utils.R")




# calculate mean and re-assign classes --------------------------------------------------------

# Process all combinations of sites and analysis types
sites <- c("BCC", "UHN")
analysis_types <- c("dlbcl90", "rhl30")

results <- list()

for (site in sites) {
  for (analysis_type in analysis_types) {
    
    tryCatch({
      results[[site]][[analysis_type]] <- process_site_analysis(site, analysis_type)
      
    }, error = function(e) {
      warning(sprintf("Error in %s %s analysis: %s", site, toupper(analysis_type), e$message))
    })
    message("")  # Add spacing between analyses
  }
}


# Summary of results
for (site in sites) {
  for (analysis_type in analysis_types) {
    if (!is.null(results[[site]][[analysis_type]])) {
      file_name <- sprintf("mydat_%s_mean_for_%s_calls.rda", site, analysis_type)
      n_samples <- nrow(results[[site]][[analysis_type]])
      message(sprintf("- %s (%d samples)", file_name, n_samples))
    }
  }
}