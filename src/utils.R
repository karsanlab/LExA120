###############################################################################
# LEXA120 Analysis Helper Functions
#
# Author         | Dr. Shujun Huang 
# Maintainer     | Joshua Bridgers (jbridgers@bcgsc.ca)
# PI (Lab)       | Dr. Aly Karsan (Karsan Lab)
#                |   https://www.bcgsc.ca/labs/karsan-lab
# 
# Date Created   | 2021-12-12
# Description    | Helper functions for LEXA120 analysis tasks
#
# License        | MIT License (see https://opensource.org/licenses/MIT)
#
###############################################################################


# data processing helper functions ------------------------------------------------------------

#' Read and process RCC files from multiple paths
#' @param rcc_paths List of paths containing RCC files
#' @param sample_naming_function Function to generate sample names from file paths
#' @param prefix_removal Character string to remove from column names (e.g., "AK.")
#' @return Matrix with processed expression data
process_rcc_files <- function(rcc_paths, sample_naming_function = NULL, prefix_removal = NULL) {
  
  # Handle single path vs multiple paths
  if (length(rcc_paths) == 1) {
    filesIn <- nanostringr::read_rcc(path = rcc_paths[[1]])
    csvfile <- data.frame(filesIn$raw)
    rownames(csvfile) <- csvfile$Name
    csvfile <- csvfile[, -c(1:3)]
    
    if (!is.null(sample_naming_function)) {
      colnames(csvfile) <- sample_naming_function(colnames(csvfile))
    } else if (!is.null(prefix_removal)) {
      colnames(csvfile) <- gsub(prefix_removal, "", colnames(csvfile))
    }
    
  } else {
    # Handle multiple paths (batches)
    csvfile_list <- list()
    
    for (i in seq_along(rcc_paths)) {
      filesIn <- nanostringr::read_rcc(path = rcc_paths[[i]])
      temp_csvfile <- data.frame(filesIn$raw)
      rownames(temp_csvfile) <- temp_csvfile$Name
      temp_csvfile <- temp_csvfile[, -c(1:3)]
      
      # Apply naming function if provided
      if (!is.null(sample_naming_function)) {
        if (is.list(sample_naming_function)) {
          colnames(temp_csvfile) <- sample_naming_function[[i]](colnames(temp_csvfile))
        } else {
          colnames(temp_csvfile) <- sample_naming_function(colnames(temp_csvfile))
        }
      } else if (!is.null(prefix_removal)) {
        colnames(temp_csvfile) <- gsub(prefix_removal, "", colnames(temp_csvfile))
      }
      
      csvfile_list[[i]] <- temp_csvfile
    }
    
    # Merge all csvfiles
    csvfile <- csvfile_list[[1]]
    if (length(csvfile_list) > 1) {
      for (i in 2:length(csvfile_list)) {
        if (identical(rownames(csvfile), rownames(csvfile_list[[i]]))) {
          csvfile <- bind_cols(csvfile, csvfile_list[[i]])
        } else {
          csvfile <- merge(csvfile, csvfile_list[[i]], by = "row.names", all = FALSE)
          rownames(csvfile) <- csvfile$Row.names
          csvfile <- csvfile[, -1]
        }
      }
    }
  }
  
  # Replace values < 1 with 1 to avoid log(0) = inf
  tmp1 <- which(csvfile < 1, arr.ind = TRUE)
  if (nrow(tmp1) > 0) {
    csvfile[tmp1] <- 1
  }
  
  return(csvfile)
}

#' Perform DLBCL90 classification analysis
#' @param csvfile Expression matrix
#' @param round_name Name of the current round
#' @param geomean_cut Geometric mean cutoff for DLBCL90calls
#' @return DLBCL90 classification results
perform_dlbcl90_analysis <- function(csvfile, round_name, geomean_cut = NULL) {
  
  # Write expression data
  output_file <- here::here("cache", sprintf("lexa_expression_counts_for_DLBCL90_tasks_%s.csv", round_name))
  write.csv(csvfile, output_file)
  
  # Perform DLBCL90 calls
  calls_output_file <- here::here("cache", sprintf("lexa_calls_for_DLBCL90_tasks_%s.csv", round_name))
  
  if (!is.null(geomean_cut)) {
    lexa120_calls <- DLBCL90calls(output_file,
                                  nHeading = 1,
                                  outfileName = calls_output_file,
                                  geomeanCut = geomean_cut)
  } else {
    lexa120_calls <- DLBCL90calls(output_file,
                                  nHeading = 1,
                                  outfileName = calls_output_file)
  }
  
  lexa120_calls <- as_tibble(lexa120_calls)
  colnames(lexa120_calls) <- paste0(colnames(lexa120_calls), ".lexa120")
  lexa120_calls$sampleName <- lexa120_calls$sampleName.lexa120
  
  return(lexa120_calls)
}

#' Perform RHL30 classification analysis
#' @param csvfile Expression matrix
#' @param hk_genes Housekeeping genes for normalization
#' @param rhl30_model_df RHL30 model data frame
#' @return RHL30 classification results
perform_rhl30_analysis <- function(csvfile, hk_genes, rhl30_model_df) {
  
  csvfile <- as.matrix(csvfile)
  
  # Normalize the expression matrix
  sample_normalizer_values <- get_sample_normalizer_value(csvfile, hk_genes)
  
  # Check for samples below threshold
  low_samples <- names(sample_normalizer_values[sample_normalizer_values <= 35])
  if (length(low_samples) > 0) {
    message("Samples below normalization threshold (35): ", paste(low_samples, collapse = ", "))
  }
  
  exprs_mat_norm <- normalize_exprs_mat(csvfile, sample_normalizer_values)
  
  # Calculate RHL30 scores
  rhl30_df <- get_rhl30_scores_df(exprs_mat_norm, rhl30_model_df)
  
  lexa120_calls <- inner_join(rhl30_df,
                              data.frame(sample_id = names(sample_normalizer_values), 
                                         normalizer = sample_normalizer_values),
                              by = "sample_id") %>%
    mutate(risk_class = ifelse(score > 10.4, "High", "Low")) %>%
    dplyr::rename(sampleName = sample_id) %>%
    mutate(sampleName = gsub("X", "", sampleName))
  
  colnames(lexa120_calls) <- paste0(colnames(lexa120_calls), ".lexa120")
  
  return(lexa120_calls)
}

#' Merge analysis results with sample IDs
#' @param ids_data Sample ID mapping data
#' @param reference_calls Reference method results (DLBCL90 or RHL30)
#' @param lexa_calls LExA method results
#' @param round_col Column name for the round in ids_data
#' @param reference_id_col Column name for reference ID matching
#' @param sample_type_filter Sample type to include (optional)
#' @param id_transformation Function to transform IDs if needed
#' @return Merged dataset
merge_analysis_results <- function(ids_data, reference_calls, lexa_calls, 
                                   round_col, reference_id_col, 
                                   sample_type_filter = NULL, 
                                   id_transformation = NULL) {
  
  # Filter IDs data for current round
  ids_filtered <- ids_data %>% filter(!is.na(.data[[round_col]]))
  
  # Apply sample type filter if specified
  if (!is.null(sample_type_filter)) {
    ids_filtered <- ids_filtered %>% filter(SampleType %in% sample_type_filter)
  }
  
  # Merge with reference calls
  mydat <- dplyr::inner_join(ids_filtered, 
                             reference_calls, 
                             by = c("ResID" = reference_id_col))
  
  # Apply ID transformation if needed
  if (!is.null(id_transformation)) {
    mydat <- mydat %>% mutate(!!round_col := id_transformation(.data[[round_col]]))
  }
  
  # Merge with LExA calls
  mydat <- dplyr::inner_join(mydat, lexa_calls, 
                             by = c(round_col = "sampleName.lexa120"))
  
  return(mydat)
}

#' Calculate classification performance metrics
#' @param data Merged dataset
#' @param score_cols List with reference and lexa score column names
#' @param class_cols List with reference and lexa class column names
#' @param set_name Name of the dataset
#' @return Performance metrics vector
calculate_performance_metrics <- function(data, score_cols, class_cols, set_name) {
  
  # Calculate correlation and linear model
  corr <- cor(data[[score_cols$lexa]], data[[score_cols$reference]])
  lm_model <- lm(data[[score_cols$reference]] ~ data[[score_cols$lexa]])
  ci <- confint(lm_model, names(lm_model$coefficients)[2], level = 0.95)
  
  # Create results vector
  results <- c(
    Set = set_name,
    Size = nrow(data),
    Pearson = round(corr, 3),
    Intercept = round(summary(lm_model)$coefficients[1, 1], 3),
    Slope = round(summary(lm_model)$coefficients[2, 1], 3),
    CI = paste0(round(ci[1], 3), "~", round(ci[2], 3)),
    R_squared = round(summary(lm_model)$adj.r.squared, 3),
    Misclassification = sum(data[[class_cols$reference]] != data[[class_cols$lexa]])
  )
  
  return(results)
}

#' Identify and report misclassified samples
#' @param data Merged dataset
#' @param class_cols List with reference and lexa class column names
#' @param sample_id_col Column name for sample identifiers
#' @return Vector of misclassified sample names
report_misclassifications <- function(data, class_cols, sample_id_col = "OriginalSourceName") {
  
  misclassified <- data[[sample_id_col]][data[[class_cols$reference]] != data[[class_cols$lexa]]]
  
  if (length(misclassified) > 0) {
    message("Misclassified samples: ", paste(misclassified, collapse = ", "))
  } else {
    message("No misclassified samples")
  }
  
  return(misclassified)
}

#' Save analysis results
#' @param data Merged dataset
#' @param results_list List of performance metrics
#' @param round_name Name of the round
#' @param analysis_type Type of analysis ("dlbcl90" or "rhl30")
save_analysis_results <- function(data, results_list, round_name, analysis_type) {
  
  # Create variable names
  data_var_name <- paste0("mydat_", gsub("_", "_", round_name), "_for_", analysis_type, "_calls")
  
  # Create environment for saving
  save_env <- new.env()
  assign(data_var_name, data, envir = save_env)
  
  # Assign results to environment
  for (i in seq_along(results_list)) {
    result_var_name <- paste0(names(results_list)[i], "_", gsub("_", "_", round_name))
    assign(result_var_name, results_list[[i]], envir = save_env)
  }
  
  # Save file
  output_file <- here::here("data", "generated_lps", 
                            sprintf("analysis_for_%s_to_%s.rda", round_name, analysis_type))
  
  save(list = ls(envir = save_env), file = output_file, envir = save_env)
  
  # Also save CSV
  csv_file <- here::here("results", sprintf("Merged_LExA_and_%s_results_%s.csv", 
                                            toupper(analysis_type), round_name))
  write_csv(data, csv_file)
}



# mean calculation helper functions -----------------------------------------------------------

# Helper function to load multiple analysis results
load_analysis_results <- function(site, rounds, analysis_type) {
  results <- list()
  for (round in rounds) {
    file_path <- here::here("data", "generated_lps", 
                            sprintf("analysis_for_%s_round%d_to_%s.rda", site, round, analysis_type))
    load(file_path, envir = .GlobalEnv)
    
    # Store the loaded data object name for reference
    data_var_name <- sprintf("mydat_%s_round%d_for_%s_calls", site, round, analysis_type)
    results[[paste0("round", round)]] <- get(data_var_name)
  }
  return(results)
}

# Generic function to calculate means across rounds for any analysis type
calculate_rounds_means <- function(round_data_list, analysis_type, score_columns, class_thresholds) {
  
  # Extract ResID from first round for merging
  base_data <- round_data_list[[1]] %>% 
    select(ResID, OriginalSourceName, SampleType, 
           LibraryID_Round1, LibraryID_Round2, LibraryID_Round3, 
           UHN_Round1, UHN_Round2, UHN_Round3,
           ends_with(paste0(".", analysis_type)))
  
  # Combine score data from all rounds
  combined_scores <- bind_cols(
    round_data_list[[1]] %>% 
      select(ResID, ends_with(".lexa120")) %>% 
      column_to_rownames(var = "ResID") %>% 
      setNames(paste0(names(.), '.round1')),
    
    round_data_list[[2]] %>% 
      select(ResID, ends_with(".lexa120")) %>% 
      column_to_rownames(var = "ResID") %>% 
      setNames(paste0(names(.), '.round2')),
    
    round_data_list[[3]] %>% 
      select(ResID, ends_with(".lexa120")) %>% 
      column_to_rownames(var = "ResID") %>% 
      setNames(paste0(names(.), '.round3'))
  ) %>% 
    rownames_to_column("ResID")
  
  # Calculate means for each score column
  rounds_means <- combined_scores %>% select(ResID)
  
  for (score_col in score_columns) {
    mean_col_name <- paste0(score_col, ".lexa120_mean")
    rounds_means[[mean_col_name]] <- rowMeans(
      select(combined_scores, starts_with(paste0(score_col, ".lexa120")))
    )
  }
  
  # Apply classification thresholds
  if (analysis_type == "dlbcl90") {
    rounds_means <- rounds_means %>%
      mutate(
        DLBCLcall.lexa120_mean = case_when(
          DLBCLp.lexa120_mean > 0.9 ~ "ABC",
          DLBCLp.lexa120_mean < 0.1 ~ "GCB",
          TRUE ~ "UNCLASS"
        ),
        PMBLcall.lexa120_mean = case_when(
          PMBLp.lexa120_mean > 0.9 ~ "PMBL",
          PMBLp.lexa120_mean < 0.1 ~ "DLBCL",
          TRUE ~ "Unclear"
        ),
        DHITsig_class.lexa120_mean = case_when(
          DHITsig_prob_POS.lexa120_mean > 0.8 ~ "POS",
          DHITsig_prob_POS.lexa120_mean < 0.2 ~ "NEG",
          TRUE ~ "UNCLASS"
        )
      )
  } else if (analysis_type == "rhl30") {
    rounds_means <- rounds_means %>%
      mutate(
        risk_class.lexa120_mean = ifelse(score.lexa120_mean > 10.4, "High", "Low")
      )
  }
  
  # Merge with base data and clean up column names
  final_data <- inner_join(base_data, rounds_means, by = "ResID") %>%
    select(!matches("round1|round2|round3"))
  
  colnames(final_data) <- gsub("_mean", "", colnames(final_data))
  
  return(final_data)
}

# Function to process and save results for a site and analysis type
process_site_analysis <- function(site, analysis_type) {
  
  message(sprintf("Processing %s analysis for %s", toupper(analysis_type), toupper(site)))
  
  # Load analysis results for rounds 1-3
  round_data_list <- load_analysis_results(site, 1:3, analysis_type)
  
  # Define score columns and thresholds based on analysis type
  if (analysis_type == "dlbcl90") {
    score_columns <- c("DLBCLval", "DLBCLp", "PMBLval", "PMBLp", 
                       "DHITsig_score", "DHITsig_prob_POS", "DHITsig_prob_NEG")
    class_thresholds <- list(
      dlbcl = c(abc_threshold = 0.9, gcb_threshold = 0.1),
      pmbl = c(pmbl_threshold = 0.9, dlbcl_threshold = 0.1),
      dhit = c(pos_threshold = 0.8, neg_threshold = 0.2)
    )
  } else if (analysis_type == "rhl30") {
    score_columns <- c("score")
    class_thresholds <- list(risk = c(high_threshold = 10.4))
  }
  
  # Calculate means across rounds
  mean_results <- calculate_rounds_means(round_data_list, analysis_type, 
                                         score_columns, class_thresholds)
  
  # Display classification results
  if (analysis_type == "dlbcl90") {
    message("DLBCL COO calls:")
    print(select(mean_results, starts_with("DLBCLcall")))
    message("PMBL calls:")
    print(select(mean_results, starts_with("PMBLcall")))
    message("DHIT calls:")
    print(select(mean_results, starts_with("DHITsig_class")))
  } else if (analysis_type == "rhl30") {
    message("Risk classification calls:")
    print(select(mean_results, contains("score")))
    print(select(mean_results, starts_with("risk_class")))
  }
  
  # Save results using consistent naming
  output_var_name <- sprintf("mydat_%s_mean_for_%s_calls", site, analysis_type)
  output_file <- here::here("data", "generated_lps", paste0(output_var_name, ".rda"))
  
  # Create environment for saving to avoid conflicts
  save_env <- new.env()
  assign(output_var_name, mean_results, envir = save_env)
  
  save(list = output_var_name, file = output_file, envir = save_env)

  return(mean_results)
}


# calibration helper functions ----------------------------------------------------------------

#' Derive calibration factors for DLBCL90 classifications
#' @return List of linear models for COO, PMBL, and DHIT
derive_dlbcl90_calibration_factors <- function() {
  
  # Combine data from rounds 4-5 for calibration
  mydat <- bind_rows(
    mydat_BCC_round4_for_dlbcl90_calls,
    mydat_BCC_round5_for_dlbcl90_calls
  )
  
  # Derive correction factors for each classification
  lm_coo <- lm(DLBCLval.dlbcl90 ~ DLBCLval.lexa120, data = mydat)
  lm_pmbl <- lm(PMBLval.dlbcl90 ~ PMBLval.lexa120, data = mydat)
  lm_dhit <- lm(DHITsig_score.dlbcl90 ~ DHITsig_score.lexa120, 
                data = mydat[!is.na(mydat$DHITsig_score.dlbcl90), ])
  
  return(list(coo = lm_coo, pmbl = lm_pmbl, dhit = lm_dhit))
}

#' Derive calibration factors for RHL30 risk classification
#' @return Linear model for risk classification
derive_rhl30_calibration_factors <- function() {
  
  # Use round 5 data (no cHL samples in BCC round 4)
  mydat <- mydat_BCC_round5_for_rhl30_calls
  
  # Get correction factor for cHL risk classification
  lm_risk <- lm(score.rhl30 ~ score.lexa120, data = mydat)
  
  return(lm_risk)
}

# function to get DHITsig call from LPS score
get_DHITcall <- function(LPS_score) {
  # parameters from extdata of package DLBCL90
  testLPSmean <- 2.470818
  refLPSmean <- -23.051925
  testLPSsd <- 10.133290
  refLPSsd <- 8.483427
  
  classProbCut <- 0.8
  
  testGroup <- "POS"
  refGroup <- "NEG"
  
  LPS_prob_test <- getProb(LPS_score, groupMeans = c(testLPSmean, refLPSmean), groupSds = c(testLPSsd, refLPSsd))
  LPS_prob_ref <- getProb(LPS_score, groupMeans = c(refLPSmean, testLPSmean), groupSds = c(refLPSsd, testLPSsd))
  
  LPS_class <- rep("UNCLASS", length(LPS_score))
  LPS_class[which(LPS_prob_test >= classProbCut)] <- testGroup
  LPS_class[which(LPS_prob_ref >= classProbCut)] <- refGroup
  
  return(LPS_class)
}


# function to get COO call from LPS score
get_DLBCLcall <- function(DLBCLval) {
  load(system.file("extdata", "model.table.rda", package = "DLBCL90"))
  model.table <- model.table
  paramx <- model.table[65:75, 2]
  
  index.ABC.mn <- 6
  index.GCB.mn <- 7
  index.ABC.sd <- 8
  index.GCB.sd <- 9
  
  p1 <- dnorm(DLBCLval, paramx[index.ABC.mn], paramx[index.ABC.sd])
  p2 <- dnorm(DLBCLval, paramx[index.GCB.mn], paramx[index.GCB.sd])
  DLBCLp <- p1 / (p1 + p2)
  DLBCLcall <- 1 + (DLBCLp > 0.1) + (DLBCLp > 0.9)
  DLBCLcall <- c("GCB", "UNCLASS", "ABC")[DLBCLcall]
  
  return(DLBCLcall)
}


# function to get PMBL call from LPS score
get_PMBLcall <- function(PMBLval) {
  load(system.file("extdata", "model.table.rda", package = "DLBCL90"))
  model.table <- model.table
  paramx <- model.table[65:75, 2]
  
  index.DLBCL.mn <- 2
  index.PMBL.mn <- 3
  index.DLBCL.sd <- 4
  index.PMBL.sd <- 5
  
  p1 <- dnorm(PMBLval, paramx[index.DLBCL.mn], paramx[index.DLBCL.sd])
  p2 <- dnorm(PMBLval, paramx[index.PMBL.mn], paramx[index.PMBL.sd])
  PMBLp <- p2 / (p1 + p2)
  
  PMBLcall <- 1 + (PMBLp > 0.1) + (PMBLp > 0.9)
  PMBLcall <- c("DLBCL", "Unclear", "PMBL")[PMBLcall]
  
  return(PMBLcall)
}

#' Apply DLBCL90 calibration to a single round
#' @param round_data Data frame containing uncorrected LExA scores
#' @param calibration_models List of linear models for calibration
#' @return List containing corrected data and performance metrics
apply_dlbcl90_calibration <- function(round_data, calibration_models) {
  
  # Apply calibration using model coefficients
  mydat <- round_data %>%
    mutate(
      DLBCLval.lexa120 = DLBCLval.lexa120 * summary(calibration_models$coo)$coefficients[2, 1] + 
        summary(calibration_models$coo)$coefficients[1, 1],
      PMBLval.lexa120 = PMBLval.lexa120 * summary(calibration_models$pmbl)$coefficients[2, 1] + 
        summary(calibration_models$pmbl)$coefficients[1, 1],
      DHITsig_score.lexa120 = DHITsig_score.lexa120 * summary(calibration_models$dhit)$coefficients[2, 1] + 
        summary(calibration_models$dhit)$coefficients[1, 1]
    ) %>%
    mutate(
      DLBCLcall.lexa120 = get_DLBCLcall(DLBCLval.lexa120),
      PMBLcall.lexa120 = get_PMBLcall(PMBLval.lexa120),
      DHITsig_class.lexa120 = get_DHITcall(DHITsig_score.lexa120)
    )
  
  # Calculate performance metrics using existing helper function
  res_coo <- calculate_performance_metrics(mydat,
                                           list(lexa = "DLBCLval.lexa120", reference = "DLBCLval.dlbcl90"),
                                           list(lexa = "DLBCLcall.lexa120", reference = "DLBCLcall.dlbcl90"),
                                           "DLBCL")
  
  res_pmbl <- calculate_performance_metrics(mydat,
                                            list(lexa = "PMBLval.lexa120", reference = "PMBLval.dlbcl90"),
                                            list(lexa = "PMBLcall.lexa120", reference = "PMBLcall.dlbcl90"),
                                            "DLBCL")
  
  res_dhit <- calculate_performance_metrics(mydat,
                                            list(lexa = "DHITsig_score.lexa120", reference = "DHITsig_score.dlbcl90"),
                                            list(lexa = "DHITsig_class.lexa120", reference = "DHITsig_class.dlbcl90"),
                                            "DLBCL")
  
  return(list(
    data = mydat,
    results = list(coo = res_coo, pmbl = res_pmbl, dhit = res_dhit)
  ))
}

#' Apply RHL30 calibration to a single round
#' @param round_data Data frame containing uncorrected LExA scores
#' @param calibration_model Linear model for risk calibration
#' @return List containing corrected data and performance metrics
apply_rhl30_calibration <- function(round_data, calibration_model) {
  
  # Apply calibration using model coefficients
  mydat <- round_data %>%
    mutate(
      score.lexa120 = score.lexa120 * summary(calibration_model)$coefficients[2, 1] + 
        summary(calibration_model)$coefficients[1, 1]
    ) %>%
    mutate(risk_class.lexa120 = ifelse(score.lexa120 > 10.4, "High", "Low"))
  
  # Calculate performance metrics using existing helper function
  res_risk <- calculate_performance_metrics(mydat,
                                            list(lexa = "score.lexa120", reference = "score.rhl30"),
                                            list(lexa = "risk_class.lexa120", reference = "risk_class.rhl30"),
                                            "cHL")
  
  return(list(
    data = mydat,
    results = list(risk = res_risk)
  ))
}

#' Process and save corrected results for a single round
#' @param round_number Round number (1, 2, or 3)
#' @param analysis_type Either "dlbcl90" or "rhl30"
#' @param corrected_results List containing corrected data and results
process_and_save_corrected_round <- function(round_number, analysis_type, corrected_results) {
  
  round_name <- paste0("BCC_round", round_number)
  
  if (analysis_type == "dlbcl90") {
    # Create variable names for DLBCL90 results
    data_var_name <- paste0("mydat_", round_name, "_for_dlbcl90_calls_corrected")
    coo_var_name <- paste0("res_coo_", round_name, "_corrected")
    pmbl_var_name <- paste0("res_pmbl_", round_name, "_corrected")
    dhit_var_name <- paste0("res_dhit_", round_name, "_corrected")
    
    # Save results
    assign(data_var_name, corrected_results$data)
    assign(coo_var_name, corrected_results$results$coo)
    assign(pmbl_var_name, corrected_results$results$pmbl)
    assign(dhit_var_name, corrected_results$results$dhit)
    
    save(list = c(data_var_name, coo_var_name, pmbl_var_name, dhit_var_name),
         file = here::here("data", "generated_lps", 
                           sprintf("analysis_for_%s_to_dlbcl90_corrected.rda", round_name)))
    
  } else if (analysis_type == "rhl30") {
    # Create variable names for RHL30 results
    data_var_name <- paste0("mydat_", round_name, "_for_rhl30_calls_corrected")
    risk_var_name <- paste0("res_risk_", round_name, "_corrected")
    
    # Save results
    assign(data_var_name, corrected_results$data)
    assign(risk_var_name, corrected_results$results$risk)
    
    save(list = c(data_var_name, risk_var_name),
         file = here::here("data", "generated_lps", 
                           sprintf("analysis_for_%s_to_rhl30_corrected.rda", round_name)))
  }
}


# main DLBCL90 analysis wrapper functions -----------------------------------------------------

#' Complete DLBCL90 analysis workflow
#' @param round_name Name of the analysis round
#' @param rcc_paths List of paths containing RCC files
#' @param ids_data Sample ID mapping data
#' @param dlbcl90_calls Reference DLBCL90 calls
#' @param round_col Column name for the round in ids_data
#' @param sample_naming_function Function for sample naming (optional)
#' @param prefix_removal Prefix to remove from sample names (optional)
#' @param geomean_cut Geometric mean cutoff (optional)
#' @param exclude_samples Samples to exclude (optional)
#' @param id_transformation Function to transform sample IDs (optional)
run_dlbcl90_analysis <- function(round_name, rcc_paths, ids_data, dlbcl90_calls, 
                                 round_col, sample_naming_function = NULL, 
                                 prefix_removal = "AK.", geomean_cut = NULL,
                                 exclude_samples = NULL, id_transformation = NULL) {
  
  message("Starting DLBCL90 analysis for round: ", round_name)
  
  # Process RCC files
  csvfile <- process_rcc_files(rcc_paths, sample_naming_function, prefix_removal)
  
  # Perform DLBCL90 analysis
  lexa120_calls <- perform_dlbcl90_analysis(csvfile, round_name, geomean_cut)
  
  # Merge results
  mydat <- merge_analysis_results(ids_data, dlbcl90_calls, lexa120_calls, 
                                  round_col, "ResID.dlbcl90", 
                                  sample_type_filter = c("DLBCL", "PMBL"),
                                  id_transformation = id_transformation)
  
  # Exclude samples if specified
  if (!is.null(exclude_samples)) {
    mydat <- mydat %>% filter(!OriginalSourceName %in% exclude_samples)
  }
  
  # Calculate performance metrics for each classification type
  coo_results <- calculate_performance_metrics(mydat, 
                                               list(lexa = "DLBCLval.lexa120", reference = "DLBCLval.dlbcl90"),
                                               list(lexa = "DLBCLcall.lexa120", reference = "DLBCLcall.dlbcl90"),
                                               "DLBCL")
  
  pmbl_results <- calculate_performance_metrics(mydat,
                                                list(lexa = "PMBLval.lexa120", reference = "PMBLval.dlbcl90"),
                                                list(lexa = "PMBLcall.lexa120", reference = "PMBLcall.dlbcl90"),
                                                "DLBCL")
  
  dhit_results <- calculate_performance_metrics(mydat,
                                                list(lexa = "DHITsig_score.lexa120", reference = "DHITsig_score.dlbcl90"),
                                                list(lexa = "DHITsig_class.lexa120", reference = "DHITsig_class.dlbcl90"),
                                                "DLBCL")
  
  # Report misclassifications
  message("COO misclassifications:")
  report_misclassifications(mydat, list(lexa = "DLBCLcall.lexa120", reference = "DLBCLcall.dlbcl90"))
  
  message("PMBL misclassifications:")
  report_misclassifications(mydat, list(lexa = "PMBLcall.lexa120", reference = "PMBLcall.dlbcl90"))
  
  message("DHIT misclassifications:")
  report_misclassifications(mydat, list(lexa = "DHITsig_class.lexa120", reference = "DHITsig_class.dlbcl90"))
  
  # Save results
  results_list <- list(res_coo = coo_results, res_pmbl = pmbl_results, res_dhit = dhit_results)
  save_analysis_results(mydat, results_list, round_name, "dlbcl90")
  
  return(list(data = mydat, results = results_list))
}




# main RHL30 analysis wrapper functions -------------------------------------------------------

# Handle UHN additional samples processing
process_uhn_additional <- function(config, hk_genes, rhl30_model_df) {
  
  # Read RCC files manually since this case doesn't fit standard helper
  filesIn <- nanostringr::read_rcc(path = config$rcc_paths)
  csvfile <- data.frame(filesIn$raw)
  rownames(csvfile) <- csvfile$Name
  csvfile <- csvfile[, -c(1:3)]
  colnames(csvfile) <- list.files(path = config$rcc_paths, pattern = ".RCC")
  
  # Replace values < 1 with 1
  tmp <- which(csvfile < 1, arr.ind = TRUE) 
  if(nrow(tmp) > 0) {csvfile[tmp] <- 1}
  csvfile <- as.matrix(csvfile)   
  
  # Use existing RHL30 analysis helper
  lexa120_calls <- perform_rhl30_analysis(csvfile, hk_genes, rhl30_model_df)
  
  # Special merging for UHN additional
  mydat <- dplyr::inner_join(config$ids_data, 
                             config$reference_calls, 
                             by = c("PMH_ID" = "PMHID.rhl30")) 
  
  mydat <- dplyr::inner_join(mydat, lexa120_calls, 
                             by = c("RCC_file_ID" = "sampleName.lexa120")) %>% 
    dplyr::rename(score.rhl30 = `RHL30 Score.rhl30`) %>% 
    mutate(risk_class.rhl30 = ifelse(score.rhl30 > 10.4, "High", "Low")) %>% 
    select(UHN_Additional_Round = RCC_file_ID, ResID = AMDL_ID, 
           External_ID = PMH_ID, OriginalSourceName = PMH_ID, 
           score.rhl30, risk_class.rhl30, score.lexa120, risk_class.lexa120)
  
  # Calculate performance and save results
  results <- calculate_performance_and_save(mydat, config$round_name)
  
  return(list(data = mydat, results = results))
}

# Handle BCC additional samples processing  
process_bcc_additional <- function(config, hk_genes, rhl30_model_df) {
  
  # Read RCC files manually
  filesIn <- nanostringr::read_rcc(path = config$rcc_paths)
  csvfile <- data.frame(filesIn$raw)
  rownames(csvfile) <- csvfile$Name
  csvfile <- csvfile[, -c(1:3)]
  colnames(csvfile) <- sapply(colnames(csvfile), function(xx)(unlist(strsplit(xx, "\\."))[2]))
  
  # Replace values < 1 with 1
  tmp <- which(csvfile < 1, arr.ind = TRUE) 
  if(nrow(tmp) > 0) {csvfile[tmp] <- 1}
  csvfile <- as.matrix(csvfile)   
  
  # Use existing RHL30 analysis helper
  lexa120_calls <- perform_rhl30_analysis(csvfile, hk_genes, rhl30_model_df)
  
  # Special merging for BCC additional
  mydat <- dplyr::inner_join(config$ids_data %>% filter(!is.na(Additional_Round)), 
                             config$reference_calls, 
                             by = c("ResID" = "ResID.rhl30"))
  
  mydat <- dplyr::inner_join(mydat, lexa120_calls, 
                             by = c("Additional_Round" = "sampleName.lexa120")) %>% 
    filter(!is.na(score.rhl30)) %>% 
    filter(!Additional_Round %in% config$exclude_samples)
  
  # Calculate performance and save results
  results <- calculate_performance_and_save(mydat, config$round_name)
  
  return(list(data = mydat, results = results))
}

process_single_round <- function(config, hk_genes, rhl30_model_df) {
  
  round_name <- config$round_name
  
  # Handle special processing cases that don't fit the standard helper pattern
  if (!is.null(config$special_processing)) {
    
    if (config$special_processing == "uhn_additional") {
      # Process UHN additional samples (24 cHL samples)
      result <- process_uhn_additional(config, hk_genes, rhl30_model_df)
      
    } else if (config$special_processing == "bcc_additional") {
      # Process BCC additional samples
      result <- process_bcc_additional(config, hk_genes, rhl30_model_df)
      
    }
  } else {
    # Use the existing run_rhl30_analysis helper function for standard cases
    result <- run_rhl30_analysis(
      round_name = round_name,
      rcc_paths = config$rcc_paths,
      ids_data = config$ids_data,
      rhl30_calls = config$reference_calls,
      hk_genes = hk_genes,
      rhl30_model_df = rhl30_model_df,
      round_col = config$round_col,
      sample_naming_function = config$sample_naming_function,
      exclude_samples = config$exclude_samples
    )
  }
  
  return(result)
}

# Calculate performance metrics and save results (using existing helper functions)
calculate_performance_and_save <- function(mydat, round_name) {
  
  # Use existing helper function for performance calculation
  res_risk <- calculate_performance_metrics(mydat,
                                            list(lexa = "score.lexa120", reference = "score.rhl30"),
                                            list(lexa = "risk_class.lexa120", reference = "risk_class.rhl30"),
                                            "cHL")
  
  # Report mis-classifications using existing helper
  report_misclassifications(mydat, 
                            list(lexa = "risk_class.lexa120", reference = "risk_class.rhl30"))
  
  # Save results using existing helper
  results_list <- list(res_risk = res_risk)
  save_analysis_results(mydat, results_list, round_name, "rhl30")
  
  return(results_list)
}

#' Complete RHL30 analysis workflow
#' @param round_name Name of the analysis round
#' @param rcc_paths List of paths containing RCC files
#' @param ids_data Sample ID mapping data
#' @param rhl30_calls Reference RHL30 calls
#' @param hk_genes Housekeeping genes for normalization
#' @param rhl30_model_df RHL30 model data frame
#' @param round_col Column name for the round in ids_data
#' @param sample_naming_function Function for sample naming (optional)
#' @param prefix_removal Prefix to remove from sample names (optional)
#' @param exclude_samples Samples to exclude (optional)
#' @param id_transformation Function to transform sample IDs (optional)
run_rhl30_analysis <- function(round_name, rcc_paths, ids_data, rhl30_calls, 
                               hk_genes, rhl30_model_df, round_col, 
                               sample_naming_function = NULL, prefix_removal = NULL,
                               exclude_samples = "AK1133", id_transformation = NULL) {
  
  message("Starting RHL30 analysis for round: ", round_name)
  
  # Process RCC files
  csvfile <- process_rcc_files(rcc_paths, sample_naming_function, prefix_removal)
  
  # Perform RHL30 analysis
  lexa120_calls <- perform_rhl30_analysis(csvfile, hk_genes, rhl30_model_df)
  
  # Merge results
  mydat <- merge_analysis_results(ids_data, rhl30_calls, lexa120_calls, 
                                  round_col, "ResID.rhl30", 
                                  sample_type_filter = "cHL",
                                  id_transformation = id_transformation)
  
  # Exclude samples if specified
  if (!is.null(exclude_samples)) {
    mydat <- mydat %>% filter(!OriginalSourceName %in% exclude_samples)
  }
  
  # Calculate performance metrics
  risk_results <- calculate_performance_metrics(mydat,
                                                list(lexa = "score.lexa120", reference = "score.rhl30"),
                                                list(lexa = "risk_class.lexa120", reference = "risk_class.rhl30"),
                                                "cHL")
  
  # Report mis-classifications
  message("Risk classification misclassifications:")
  report_misclassifications(mydat, list(lexa = "risk_class.lexa120", reference = "risk_class.rhl30"))
  
  # Save results
  results_list <- list(res_risk = risk_results)
  save_analysis_results(mydat, results_list, round_name, "rhl30")
  
  return(list(data = mydat, results = results_list))
}
