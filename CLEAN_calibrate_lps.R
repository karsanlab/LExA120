## This script aims to derive the calibration factors between LExA with the two historical panels (DLBCL90 and RHL30) using BCC Rounds 4-5
## and then apply the calibrators to BCC Rounds 1-3
rm(list = ls())
library(tidyverse)
library(here)
library(miceadds)
source.all(here::here("report", "functions"))


##############################################################################################
###### Load the uncorrected LExA LPSs
##############################################################################################
## BCC round 1 KARSANBIO-2417
load(here::here("data", "generated_lps", "analysis_for_BCC_round1_to_dlbcl90.rda"))
load(here::here("data", "generated_lps", "analysis_for_BCC_round1_to_rhl30.rda"))
## BCC round 2 KARSANBIO-2357
load(here::here("data", "generated_lps", "analysis_for_BCC_round2_to_dlbcl90.rda"))
load(here::here("data", "generated_lps", "analysis_for_BCC_round2_to_rhl30.rda"))
## BCC round 3 KARSANBIO-2438
load(here::here("data", "generated_lps", "analysis_for_BCC_round3_to_dlbcl90.rda"))
load(here::here("data", "generated_lps", "analysis_for_BCC_round3_to_rhl30.rda"))
## Round 4 KARSANBIO-2545
load(here::here("data", "generated_lps", "analysis_for_BCC_round4_to_dlbcl90.rda"))
load(here::here("data", "generated_lps", "analysis_for_BCC_round4_to_rhl30.rda"))
# Round 5
load(here::here("data", "generated_lps", "analysis_for_BCC_additional_round_to_dlbcl90.rda"))
load(here::here("data", "generated_lps", "analysis_for_BCC_additional_round_to_rhl30.rda"))




##############################################################################################
###### Derive COO, PMBL, and DHIT correction factors
##############################################################################################
## The data used to determine the LExA/DLBCL calibration factors
mydat = bind_rows(mydat_BCC_round4_for_dlbcl90_calls, mydat_BCC_additional_round_for_dlbcl90_calls) 

## Get the correction factors for each classification
lm_coo = lm(DLBCLval.dlbcl90 ~ DLBCLval.lexa120, data = mydat)
lm_equation_coo = format_lm_equation(lm_coo)
# [1] "y = 1.078x - 5.611"    BCC R4 + BCC addtional ** Final calibration factors

lm_pmbl = lm(PMBLval.dlbcl90 ~ PMBLval.lexa120, data = mydat)
lm_equation_pmbl = format_lm_equation(lm_pmbl)
# [1] "y = 1.01x + 5.559"  BCC R4 + BCC addtional ** Final calibration factors

lm_dhit = lm(DHITsig_score.dlbcl90 ~ DHITsig_score.lexa120, data = mydat[!is.na(mydat$DHITsig_score.dlbcl90),])
lm_equation_dhit = format_lm_equation(lm_dhit)
# [1] "y = 1.038x + 2.401" BCC R4 + BCC addtional ** Final calibration factors


##############################################################################################
######  Derive cHL risk correction factors
##############################################################################################
## The data used to determine the LExA/RHL30 calibration factors
mydat1 = mydat_BCC_additional_round_for_rhl30_calls 

### Get the correction factors for cHL risk classification
lm_risk = lm(score.rhl30 ~ score.lexa120, data = mydat1)
lm_equation_risk = format_lm_equation(lm_risk)
#[1] "y = 1.004x + 0.231" BCC addtional ** Final calibration factors

save(lm_coo, lm_pmbl, lm_dhit, lm_risk, file = here::here("data", "lps_correction_factors.rda"))



##############################################################################################
###### Apply the COO, PMBL, and DHIT correction factors
##############################################################################################
##### Apply the correction factors to BCC Round 1 LPSs
# Calibrate BCC Round 1 LPSs
mydat = mydat_BCC_round1_for_dlbcl90_calls %>% 
    mutate(DLBCLval.lexa120 = DLBCLval.lexa120 * summary(lm_coo)$coefficients[2, 1] + summary(lm_coo)$coefficients[1, 1],
           PMBLval.lexa120 = PMBLval.lexa120 * summary(lm_pmbl)$coefficients[2, 1] + summary(lm_pmbl)$coefficients[1, 1],
           DHITsig_score.lexa120 = DHITsig_score.lexa120 * summary(lm_dhit)$coefficients[2, 1] + summary(lm_dhit)$coefficients[1, 1]) %>% 
    mutate(DLBCLcall.lexa120 = get_DLBCLcall(DLBCLval.lexa120),
           PMBLcall.lexa120 = get_PMBLcall(PMBLval.lexa120),
           DHITsig_class.lexa120 = get_DHITcall(DHITsig_score.lexa120))

# Classification performance based on the corrected LPSs
corr = cor(mydat$DLBCLval.lexa120, mydat$DLBCLval.dlbcl90)
lm = lm(DLBCLval.dlbcl90 ~ DLBCLval.lexa120, data = mydat)
ci = confint(lm, 'DLBCLval.lexa120', level = 0.95)
lm_equation = format_lm_equation(lm)
res_coo = c(Set = "DLBCL",
            Size = nrow(mydat),
            Pearson = round(corr, 3), 
            Intercept = round(summary(lm)$coefficients[1, 1], 3),
            Slope = round(summary(lm)$coefficients[2, 1], 3), 
            CI = paste0(round(ci[1], 3), "~", round(ci[2], 3)),
            R_squared = round(summary(lm)$adj.r.squared, 3),
            Misclassification = sum(mydat$DLBCLcall.dlbcl90 != mydat$DLBCLcall.lexa120))


corr = cor(mydat$PMBLval.lexa120, mydat$PMBLval.dlbcl90)
lm = lm(PMBLval.dlbcl90 ~ PMBLval.lexa120, data = mydat)
ci = confint(lm, 'PMBLval.lexa120', level = 0.95)
lm_equation = format_lm_equation(lm)
res_pmbl = c(Set = "DLBCL",
             Size = nrow(mydat),
             Pearson = round(corr, 3), 
             Intercept = round(summary(lm)$coefficients[1, 1], 3),
             Slope = round(summary(lm)$coefficients[2, 1], 3), 
             CI = paste0(round(ci[1], 3), "~", round(ci[2], 3)),
             R_squared = round(summary(lm)$adj.r.squared, 3),
             Misclassification = sum(mydat$PMBLcall.dlbcl90 != mydat$PMBLcall.lexa120))


corr = cor(mydat$DHITsig_score.lexa120, mydat$DHITsig_score.dlbcl90)
lm = lm(DHITsig_score.dlbcl90 ~ DHITsig_score.lexa120, data = mydat)
ci = confint(lm, 'DHITsig_score.lexa120', level = 0.95)
lm_equation = format_lm_equation(lm) 
res_dhit = c(Set = "DLBCL",
             Size = nrow(mydat),
             Pearson = round(corr, 3), 
             Intercept = round(summary(lm)$coefficients[1, 1], 3),
             Slope = round(summary(lm)$coefficients[2, 1], 3), 
             CI = paste0(round(ci[1], 3), "~", round(ci[2], 3)),
             R_squared = round(summary(lm)$adj.r.squared, 3),
             Misclassification = sum(mydat$DHITsig_class.dlbcl90 != mydat$DHITsig_class.lexa120))


# Output the analysis results in Rda
mydat_BCC_round1_for_dlbcl90_calls_corrected  = mydat
res_coo_BCC_round1_corrected = res_coo
res_pmbl_BCC_round1_corrected = res_pmbl 
res_dhit_BCC_round1_corrected = res_dhit

save(mydat_BCC_round1_for_dlbcl90_calls_corrected, 
     res_coo_BCC_round1_corrected, 
     res_pmbl_BCC_round1_corrected, 
     res_dhit_BCC_round1_corrected,
     file = here::here("data", "generated_lps", "analysis_for_BCC_round1_to_dlbcl90_corrected.rda"))


##### Apply the correction factors to BCC Round 2 LPSs
# Calibrate BCC Round 2 LPSs
mydat = mydat_BCC_round2_for_dlbcl90_calls %>% 
    mutate(DLBCLval.lexa120 = DLBCLval.lexa120 * summary(lm_coo)$coefficients[2, 1] + summary(lm_coo)$coefficients[1, 1],
           PMBLval.lexa120 = PMBLval.lexa120 * summary(lm_pmbl)$coefficients[2, 1] + summary(lm_pmbl)$coefficients[1, 1],
           DHITsig_score.lexa120 = DHITsig_score.lexa120 * summary(lm_dhit)$coefficients[2, 1] + summary(lm_dhit)$coefficients[1, 1]) %>% 
    mutate(DLBCLcall.lexa120 = get_DLBCLcall(DLBCLval.lexa120),
           PMBLcall.lexa120 = get_PMBLcall(PMBLval.lexa120),
           DHITsig_class.lexa120 = get_DHITcall(DHITsig_score.lexa120))

## Classification performance based on the corrected LPSs
corr = cor(mydat$DLBCLval.lexa120, mydat$DLBCLval.dlbcl90)
lm = lm(DLBCLval.dlbcl90 ~ DLBCLval.lexa120, data = mydat)
ci = confint(lm, 'DLBCLval.lexa120', level = 0.95)
lm_equation = format_lm_equation(lm)
res_coo = c(Set = "DLBCL",
            Size = nrow(mydat),
            Pearson = round(corr, 3), 
            Intercept = round(summary(lm)$coefficients[1, 1], 3),
            Slope = round(summary(lm)$coefficients[2, 1], 3), 
            CI = paste0(round(ci[1], 3), "~", round(ci[2], 3)),
            R_squared = round(summary(lm)$adj.r.squared, 3),
            Misclassification = sum(mydat$DLBCLcall.dlbcl90 != mydat$DLBCLcall.lexa120))


corr = cor(mydat$PMBLval.lexa120, mydat$PMBLval.dlbcl90)
lm = lm(PMBLval.dlbcl90 ~ PMBLval.lexa120, data = mydat)
ci = confint(lm, 'PMBLval.lexa120', level = 0.95)
lm_equation = format_lm_equation(lm)
res_pmbl = c(Set = "DLBCL",
             Size = nrow(mydat),
             Pearson = round(corr, 3), 
             Intercept = round(summary(lm)$coefficients[1, 1], 3),
             Slope = round(summary(lm)$coefficients[2, 1], 3), 
             CI = paste0(round(ci[1], 3), "~", round(ci[2], 3)),
             R_squared = round(summary(lm)$adj.r.squared, 3),
             Misclassification = sum(mydat$PMBLcall.dlbcl90 != mydat$PMBLcall.lexa120))


corr = cor(mydat$DHITsig_score.lexa120, mydat$DHITsig_score.dlbcl90)
lm = lm(DHITsig_score.dlbcl90 ~ DHITsig_score.lexa120, data = mydat)
ci = confint(lm, 'DHITsig_score.lexa120', level = 0.95)
lm_equation = format_lm_equation(lm) 
res_dhit = c(Set = "DLBCL",
             Size = nrow(mydat),
             Pearson = round(corr, 3), 
             Intercept = round(summary(lm)$coefficients[1, 1], 3),
             Slope = round(summary(lm)$coefficients[2, 1], 3), 
             CI = paste0(round(ci[1], 3), "~", round(ci[2], 3)),
             R_squared = round(summary(lm)$adj.r.squared, 3),
             Misclassification = sum(mydat$DHITsig_class.dlbcl90 != mydat$DHITsig_class.lexa120))


## Output the analysis results in Rda
mydat_BCC_round2_for_dlbcl90_calls_corrected  = mydat
res_coo_BCC_round2_corrected = res_coo
res_pmbl_BCC_round2_corrected = res_pmbl 
res_dhit_BCC_round2_corrected = res_dhit

save(mydat_BCC_round2_for_dlbcl90_calls_corrected, 
     res_coo_BCC_round2_corrected, 
     res_pmbl_BCC_round2_corrected, 
     res_dhit_BCC_round2_corrected,
     file = here::here("data", "generated_lps", "analysis_for_BCC_round2_to_dlbcl90_corrected.rda"))


##### Apply the correction factors to BCC Round 3 LPSs
## Calibrate BCC Round 2 LPSs
mydat = mydat_BCC_round3_for_dlbcl90_calls %>% 
    mutate(DLBCLval.lexa120 = DLBCLval.lexa120 * summary(lm_coo)$coefficients[2, 1] + summary(lm_coo)$coefficients[1, 1],
           PMBLval.lexa120 = PMBLval.lexa120 * summary(lm_pmbl)$coefficients[2, 1] + (summary(lm_pmbl)$coefficients[1, 1]),
           DHITsig_score.lexa120 = DHITsig_score.lexa120 * summary(lm_dhit)$coefficients[2, 1] + summary(lm_dhit)$coefficients[1, 1]) %>% 
    mutate(DLBCLcall.lexa120 = get_DLBCLcall(DLBCLval.lexa120),
           PMBLcall.lexa120 = get_PMBLcall(PMBLval.lexa120),
           DHITsig_class.lexa120 = get_DHITcall(DHITsig_score.lexa120))

## Classification performance based on the corrected LPSs
corr = cor(mydat$DLBCLval.lexa120, mydat$DLBCLval.dlbcl90)
lm = lm(DLBCLval.dlbcl90 ~ DLBCLval.lexa120, data = mydat)
ci = confint(lm, 'DLBCLval.lexa120', level = 0.95)
lm_equation = format_lm_equation(lm)
res_coo = c(Set = "DLBCL",
            Size = nrow(mydat),
            Pearson = round(corr, 3), 
            Intercept = round(summary(lm)$coefficients[1, 1], 3),
            Slope = round(summary(lm)$coefficients[2, 1], 3), 
            CI = paste0(round(ci[1], 3), "~", round(ci[2], 3)),
            R_squared = round(summary(lm)$adj.r.squared, 3),
            Misclassification = sum(mydat$DLBCLcall.dlbcl90 != mydat$DLBCLcall.lexa120))


corr = cor(mydat$PMBLval.lexa120, mydat$PMBLval.dlbcl90)
lm = lm(PMBLval.dlbcl90 ~ PMBLval.lexa120, data = mydat)
ci = confint(lm, 'PMBLval.lexa120', level = 0.95)
lm_equation = format_lm_equation(lm)
res_pmbl = c(Set = "DLBCL",
             Size = nrow(mydat),
             Pearson = round(corr, 3), 
             Intercept = round(summary(lm)$coefficients[1, 1], 3),
             Slope = round(summary(lm)$coefficients[2, 1], 3), 
             CI = paste0(round(ci[1], 3), "~", round(ci[2], 3)),
             R_squared = round(summary(lm)$adj.r.squared, 3),
             Misclassification = sum(mydat$PMBLcall.dlbcl90 != mydat$PMBLcall.lexa120))


corr = cor(mydat$DHITsig_score.lexa120, mydat$DHITsig_score.dlbcl90)
lm = lm(DHITsig_score.dlbcl90 ~ DHITsig_score.lexa120, data = mydat)
ci = confint(lm, 'DHITsig_score.lexa120', level = 0.95)
lm_equation = format_lm_equation(lm) 
res_dhit = c(Set = "DLBCL",
             Size = nrow(mydat),
             Pearson = round(corr, 3), 
             Intercept = round(summary(lm)$coefficients[1, 1], 3),
             Slope = round(summary(lm)$coefficients[2, 1], 3), 
             CI = paste0(round(ci[1], 3), "~", round(ci[2], 3)),
             R_squared = round(summary(lm)$adj.r.squared, 3),
             Misclassification = sum(mydat$DHITsig_class.dlbcl90 != mydat$DHITsig_class.lexa120))


## Output the analysis results in Rda
mydat_BCC_round3_for_dlbcl90_calls_corrected  = mydat
res_coo_BCC_round3_corrected = res_coo
res_pmbl_BCC_round3_corrected = res_pmbl 
res_dhit_BCC_round3_corrected = res_dhit

save(mydat_BCC_round3_for_dlbcl90_calls_corrected, 
     res_coo_BCC_round3_corrected, 
     res_pmbl_BCC_round3_corrected, 
     res_dhit_BCC_round3_corrected,
     file = here::here("data", "generated_lps", "analysis_for_BCC_round3_to_dlbcl90_corrected.rda"))



##############################################################################################
###### Apply the cHL risk correction factors
##############################################################################################
##### Apply the correction factors to BCC Round 1 LPSs
## Calibrate BCC Round 1 LPSs
mydat = mydat_BCC_round1_for_rhl30_calls %>%
    mutate(score.lexa120 = score.lexa120 * summary(lm_risk)$coefficients[2, 1] + summary(lm_risk)$coefficients[1, 1]) %>% 
    mutate(risk_class.lexa120 = ifelse(score.lexa120 > 10.4, "High", "Low"))

## Classification performance based on the corrected LPSs
corr = cor(mydat$score.lexa120, mydat$score.rhl30)
lm = lm(score.rhl30 ~ score.lexa120, data = mydat)
ci = confint(lm, 'score.lexa120', level = 0.95)
lm_equation = format_lm_equation(lm)
res_risk = c(Set = "cHL",
             Size = nrow(mydat),
             Pearson = round(corr, 3), 
             Intercept = round(summary(lm)$coefficients[1, 1], 3),
             Slope = round(summary(lm)$coefficients[2, 1], 3), 
             CI = paste0(round(ci[1], 3), "~", round(ci[2], 3)),
             R_squared = round(summary(lm)$adj.r.squared, 3),
             Misclassification = sum(mydat$risk_class.rhl30 != mydat$risk_class.lexa120))


## Output the analysis results in Rda
mydat_BCC_round1_for_rhl30_calls_corrected = mydat
res_risk_BCC_round1_corrected = res_risk

save(mydat_BCC_round1_for_rhl30_calls_corrected, 
     res_risk_BCC_round1_corrected, 
     file = here::here("data", "generated_lps", "analysis_for_BCC_round1_to_rhl30_corrected.rda"))


##### Apply the correction factors to BCC Round 2 LPSs
## Calibrate BCC Round 2 LPSs
mydat = mydat_BCC_round2_for_rhl30_calls %>%
    mutate(score.lexa120 = score.lexa120 * summary(lm_risk)$coefficients[2, 1] + summary(lm_risk)$coefficients[1, 1]) %>% 
    mutate(risk_class.lexa120 = ifelse(score.lexa120 > 10.4, "High", "Low"))

## Classification performance based on the corrected LPSs
corr = cor(mydat$score.lexa120, mydat$score.rhl30)
lm = lm(score.rhl30 ~ score.lexa120, data = mydat)
ci = confint(lm, 'score.lexa120', level = 0.95)
lm_equation = format_lm_equation(lm)
res_risk = c(Set = "cHL",
             Size = nrow(mydat),
             Pearson = round(corr, 3), 
             Intercept = round(summary(lm)$coefficients[1, 1], 3),
             Slope = round(summary(lm)$coefficients[2, 1], 3), 
             CI = paste0(round(ci[1], 3), "~", round(ci[2], 3)),
             R_squared = round(summary(lm)$adj.r.squared, 3),
             Misclassification = sum(mydat$risk_class.rhl30 != mydat$risk_class.lexa120))


## Output the analysis results in Rda
mydat_BCC_round2_for_rhl30_calls_corrected = mydat
res_risk_BCC_round2_corrected = res_risk

save(mydat_BCC_round2_for_rhl30_calls_corrected, 
     res_risk_BCC_round2_corrected, 
     file = here::here("data", "generated_lps", "analysis_for_BCC_round2_to_rhl30_corrected.rda"))


##### Apply the correction factors to BCC Round 3 LPSs
## Calibrate BCC Round 3 LPSs
mydat = mydat_BCC_round3_for_rhl30_calls %>%
    mutate(score.lexa120 = score.lexa120 * summary(lm_risk)$coefficients[2, 1] + summary(lm_risk)$coefficients[1, 1]) %>% 
    mutate(risk_class.lexa120 = ifelse(score.lexa120 > 10.4, "High", "Low"))

## Classification performance based on the corrected LPSs
corr = cor(mydat$score.lexa120, mydat$score.rhl30)
lm = lm(score.rhl30 ~ score.lexa120, data = mydat)
ci = confint(lm, 'score.lexa120', level = 0.95)
lm_equation = format_lm_equation(lm)
res_risk = c(Set = "cHL",
             Size = nrow(mydat),
             Pearson = round(corr, 3), 
             Intercept = round(summary(lm)$coefficients[1, 1], 3),
             Slope = round(summary(lm)$coefficients[2, 1], 3), 
             CI = paste0(round(ci[1], 3), "~", round(ci[2], 3)),
             R_squared = round(summary(lm)$adj.r.squared, 3),
             Misclassification = sum(mydat$risk_class.rhl30 != mydat$risk_class.lexa120))


## Output the analysis results in Rda
mydat_BCC_round3_for_rhl30_calls_corrected = mydat
res_risk_BCC_round3_corrected = res_risk

save(mydat_BCC_round3_for_rhl30_calls_corrected, 
     res_risk_BCC_round3_corrected, 
     file = here::here("data", "generated_lps", "analysis_for_BCC_round3_to_rhl30_corrected.rda"))


