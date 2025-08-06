###############################################################################
#
# Author         | Dr. Shujun Huang
# Maintainer     | Joshua Bridgers (jbridgers@bcgsc.ca)
# PI (Lab)       | Dr. Aly Karsan (Karsan Lab)
#                |   https://www.bcgsc.ca/labs/karsan-lab
#
# Date Created   | 2021-12-12
# Description    | Derive the calibration factors between LExA with the two
#                |   historical panels (DLBCL90 and RHL30) using BCC Rounds 4-5
#                |   and then apply the calibrators to BCC Rounds 1-3.
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




# load the uncorrected lexa lps ---------------------------------------------------------------

rounds <- paste0("BCC_round", 1:5)
panel_types <- c("dlbcl90", "rhl30")

# Load all data files dynamically
for (round in rounds) {
  for (panel in panel_types) {
    file_path <- here::here("data", "generated_lps", sprintf("analysis_for_%s_to_%s.rda", round, panel))
    load(file_path)
  }
}


# derive calibration factors using BCC rounds 4-5 ---------------------------------------------

# Derive all calibration factors
dlbcl90_factors <- derive_dlbcl90_calibration_factors() # COO, PMBL, and DHIT correction factors
rhl30_factors <- derive_rhl30_calibration_factors()     # cHL risk correction factors

# Extract individual models for backward compatibility
lm_coo <- dlbcl90_factors$coo
lm_pmbl <- dlbcl90_factors$pmbl
lm_dhit <- dlbcl90_factors$dhit
lm_risk <- rhl30_factors

# Save calibration factors
save(lm_coo, lm_pmbl, lm_dhit, lm_risk, 
     file = here::here("data", "lps_correction_factors.rda"))




# apply calibration to rounds 1-3 for DLBCL90 -------------------------------------------------

rounds_to_correct <- 1:3

# uncorrected lps for DLBCL90 tasks
dlbcl90_round_data <- list(
  mydat_BCC_round1_for_dlbcl90_calls,
  mydat_BCC_round2_for_dlbcl90_calls,
  mydat_BCC_round3_for_dlbcl90_calls
)

# Process DLBCL90 (COO, PMBL, and DHIT) corrections
for (i in rounds_to_correct) {
  corrected_results <- apply_dlbcl90_calibration(dlbcl90_round_data[[i]], dlbcl90_factors)
  process_and_save_corrected_round(i, "dlbcl90", corrected_results)
}



# apply calibration to rounds 1-3 for RHL30 ---------------------------------------------------

rounds_to_correct <- 1:3

# Uncorrected lps for RHL30 task
rhl30_round_data <- list(
  mydat_BCC_round1_for_rhl30_calls,
  mydat_BCC_round2_for_rhl30_calls,
  mydat_BCC_round3_for_rhl30_calls
)

# Process RHL30 corrections
for (i in rounds_to_correct) {
  corrected_results <- apply_rhl30_calibration(rhl30_round_data[[i]], rhl30_factors)
  process_and_save_corrected_round(i, "rhl30", corrected_results)
}

