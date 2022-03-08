## This script aims to take in the LExA expression RCC files (after nSolver QC), generate the uncorrected LPSs, 
## correct the LPSs (calibration factors derived from LExA validation BCC Rounds 4-5), and classify the samples 
## The future LExA analysis pipeline can be based on this script
## 2021-12-12

#### Overview of the script tasks
#### For DLBCL samples
# Step 1: Read in the LExA rcc file for a batch of given samples 
# Step 2: Generate the LPS scores (uncorrected)
    # Use the DLBCL90 r package to generate the COO, PMBL, and DHIT LSP scores
# Step 3: Correcte the LPS scores (Derived from BCC Rounds 4-5 data) 
    # for COO (DLBCL90 LPS ~ LExA LPS): y = 1.078x - 5.611
    # for PMBL (DLBCL90 LPS ~ LExA LPS): y = 1.01x + 5.559
    # for DHIT	(DLBCL90 LPS ~ LExA LPS): y = 1.038x + 2.401
# Step 4: Perform the classification using the corrected LPS scores
    # need to calculate the probability
    # For COO: ABC vs GCB
    # For PMBL: PMBLvs DLBCL
    # For DHIT: NEG vs POS


#### For cHL samples
# Step 1: Read in the LExA rcc file for a batch of given samples  (Let's use BCC Round 2 data as an example)
# Step 2: Generate the LPS scores (uncorrected)
        # Use the RHL30 r package
# Step 3: Correct the LPS scores (Derived from BCC Rounds 4-5 data) 
        # for cHL Risk	(RHL30 LPS ~ LExA LPS): y = 1.004x + 0.231
# Step 4: Perform the classification using the corrected LPS scores
        # high risk vs low risk (10.4 as the cutoff from the RHL30 paper)




#### Required R packages and some self-defined functions
# Install and load nanostringr r package for reading NanoString expression rcc files
# install.packages("nanostringr")
library(nanostringr)
# Install and load the DLBCL90 r package (From Aixiang Jiang)
# we will use the DLBCL90calls function from this package to generate LPS scores (related to the DLBCL90 assay)
# install.packages("/projects/karsanlab/shuang/softwares/DLBCL90_1.1.0.tar.gz", repos = NULL, type = "source")
library(DLBCL90)

# Install and load the RHL30 r package by Chan (who developed the RHL30 signature)
# we will use the get_sample_normalizer_value + normalize_exprs_mat + get_rhl30_scores_df functions from this package 
# to generate cHL risk LPS scores (related to the RHL30 assay) 
# https://github.com/tinyheero/RHL30
# devtools::install_github("tinyheero/RHL30")
library(RHL30)

# Load the other packages for general data manipulation
library(tidyverse)
library(here)


#### Load the calibration factores between DLBCL90/RHL30 and LExA
# see CLEAN_calibrate_lps.R for more details
load(here::here("data", "lps_correction_factors.rda"))
lm_coo
lm_pmbl
lm_dhit
lm_risk

#### Let's use the 8 cHL samples and 16 DLBCL samples in BCC Round 2 as an example to perform the analyses step by step
## the current round
round = "BCC_Round_2"



#############################################################################
##  For cHL samples
#############################################################################
#### Step 1: Read in LExA NanoString data of cHL samples ##############################
## Here I used the cHLs we ran in LExA validation Round 2 as an example
## RCC files of the 8 cHL samples (together 4 PMBLs) run on 2020-05-14
filesIn = nanostringr::read_rcc(path = here::here("data", "nanostring", "round2_20200514_D70132_to_D70143_RCC"))
csvfile = data.frame(filesIn$raw)
rownames(csvfile) = csvfile$Name
csvfile = csvfile[, -c(1:3)]
csvfile = csvfile[, c(5:12)] # focusing on only cHLs
colnames(csvfile) = sapply(colnames(csvfile), function(x)(unlist(strsplit(x, "\\."))[2]))

## to avoid log(0) is inf, for all 0, change to 1
tmp = which(csvfile < 1, arr.ind = TRUE) 
if(dim(tmp)[1] > 0) {csvfile[tmp] = 1}
csvfile = as.matrix(csvfile)   


#### Step 2: Calculate the uncorrected LPSs using RHL30 r package ##############################
## RHL30 model coefficients and genes
rhl30_model_df = get_rhl30_model_coef_df()
hk_genes = filter(rhl30_model_df, gene_type == "housekeeper") %>% pull("gene_name")

## Normalize the expression matrix
sample_normalizer_values = get_sample_normalizer_value(csvfile, hk_genes)
## In the paper, a threshold of 35 was set to exclude poor quality samples. 
# This was done because very low normalizer values often lead to very high normalized expression values. 
# here we are not applying this threshold since we have our nSolver QC check before we this LPS calculation task
# in the QC check, the low-quality samples are already removed
# but we can still take a look at any samples passing the QC check but still failing this 35 normalizer threshold
names(sample_normalizer_values[sample_normalizer_values <= 35]) 
exprs_mat_norm = normalize_exprs_mat(csvfile, sample_normalizer_values)

## Calcuate the LPSs and assign the classes (uncorrected)
risk_score = get_rhl30_scores_df(exprs_mat_norm, rhl30_model_df)

lexa_calls_for_rhl30 = inner_join(risk_score, 
                                  data.frame(sample_id = names(sample_normalizer_values), normalizer = sample_normalizer_values),
                                  by = "sample_id") %>% 
    mutate(risk_class = ifelse(score > 10.4, "High", "Low")) %>% 
    dplyr::rename(sampleName = sample_id) %>% 
    mutate(sampleName = gsub("X", "", sampleName))

colnames(lexa_calls_for_rhl30) = paste0(colnames(lexa_calls_for_rhl30), ".lexa")


#### Steps 3-4: Correct the LPS scores and re-assgine the class ##############################
## Simply apply the slope and intercept (from the function of RHL30 LPSs ~ LExA LPSs) to our obtained LExA LPSs
## Re-assign the high- vs low-risk classes using the corrected LPSs
corrected_lexa_calls_for_rhl30 = lexa_calls_for_rhl30 %>%
    mutate(score.lexa = score.lexa * summary(lm_risk)$coefficients[2, 1] + summary(lm_risk)$coefficients[1, 1]) %>% 
    mutate(risk_class.lexa = ifelse(score.lexa > 10.4, "High", "Low"))


## Output the LPSs + classes for the cHLs
write.csv(corrected_lexa_calls_for_rhl30, 
          file = here::here("results",sprintf("analysis_for_%s_to_rhl30_corrected.csv", round)))





#############################################################################
##  For DLBCL samples
#############################################################################
#### Step 1: Read in LExA NanoString data of DLBCL samples ##############################
## Here I used the DLBCLs we ran in LExA validation Round 2 as an example
#### RCC files for the 12 DLBCLs run on 2020-05-13
filesIn = nanostringr::read_rcc(path = here::here("data", "nanostring", "round2_20200513_D70120_to_D70131_RCC"))
csvfile = data.frame(filesIn$raw)
rownames(csvfile) = csvfile$Name
csvfile = csvfile[, -c(1:3)]
colnames(csvfile) = sapply(colnames(csvfile), function(xx)(unlist(strsplit(xx, "\\."))[3]))

#### RCC files for two of the 4 PMBLs run on 2020-05-14
filesIn1 = nanostringr::read_rcc(path = here::here("data", "nanostring", "round2_20200514_D70132_to_D70143_RCC"))
csvfile1 = data.frame(filesIn1$raw)
rownames(csvfile1) = csvfile1$Name
csvfile1 = csvfile1[, -c(1:3)]
csvfile1 = csvfile1[, c(1:4)]
colnames(csvfile1) = sapply(colnames(csvfile1), function(xx)(unlist(strsplit(xx, "\\."))[2]))

#### Merge the two csv files
csvfile_all = merge(csvfile, csvfile1, by = "row.names", all = F)
rownames(csvfile_all) = csvfile_all$Row.names
csvfile_all = csvfile_all[, -1]

#### to avoid log(0) is inf, for all 0, change to 1
tmp1 = which(csvfile_all < 1, arr.ind = TRUE) 
if(dim(tmp1)[1] > 0) {csvfile_all [tmp1] = 1}

write.csv(csvfile_all, here::here("cache", sprintf("lexa_expression_counts_for_DLBCL90_tasks_%s.csv", round)))



#### Step 2: Calculate the uncorrected LPSs using DLBCL90 r package ##############################
#### Calculate the COO, PMBL, and DHIT LPSs and assign the classes using the DLBCL90 package (uncorrected)
## How the DLBCL90calls function work: https://www.bcgsc.ca/jira/browse/KARSANBIO-2334
expression_file = here::here("cache", sprintf("lexa_expression_counts_for_DLBCL90_tasks_%s.csv", round))
lexa_calls_for_dlbcl90 = DLBCL90calls(expression_file,
                                      nHeading = 1, 
                                      outfileName = here::here("cache", sprintf("lexa_calls_for_DLBCL90_tasks_%s.csv", round))) 

lexa_calls_for_dlbcl90 = as_tibble(lexa_calls_for_dlbcl90)
colnames(lexa_calls_for_dlbcl90) = paste0(colnames(lexa_calls_for_dlbcl90), ".lexa")



#### Steps 3-4: Correct the LPS scores and re-assgine the class ##############################
## Simply apply the slope and intercept (from the function of DLBCL90 LPSs ~ LExA LPSs) to our obtained LExA LPSs
## Re-assign the COO, PMBL, and DHIT classes using the corrected LPSs using the three self-defined functions
source('/projects/karsanlab/shuang/all_lexa120_validation_analysis/report/functions/get_DLBCLcall.R') # assign ABC vs GCB
source('/projects/karsanlab/shuang/all_lexa120_validation_analysis/report/functions/get_PMBLcall.R')  # assign PMBL vs DLBCL
source('/projects/karsanlab/shuang/all_lexa120_validation_analysis/report/functions/get_DHITcall.R')  # assign POS vs NEG
corrected_lexa_calls_for_dlbcl90 = lexa_calls_for_dlbcl90 %>% 
    mutate(DLBCLval.lexa = DLBCLval.lexa * summary(lm_coo)$coefficients[2, 1] + summary(lm_coo)$coefficients[1, 1],
           PMBLval.lexa = PMBLval.lexa * summary(lm_pmbl)$coefficients[2, 1] + summary(lm_pmbl)$coefficients[1, 1],
           DHITsig_score.lexa = DHITsig_score.lexa * summary(lm_dhit)$coefficients[2, 1] + summary(lm_dhit)$coefficients[1, 1]) %>% 
    mutate(DLBCLcall.lexa = get_DLBCLcall(DLBCLval.lexa),
           PMBLcall.lexa = get_PMBLcall(PMBLval.lexa),
           DHITsig_class.lexa = get_DHITcall(DHITsig_score.lexa))



## Output the LPSs + classes for the cHLs
write.csv(corrected_lexa_calls_for_dlbcl90, 
          file = here::here("results",sprintf("analysis_for_%s_to_dlbcl90_corrected.csv", round)))



