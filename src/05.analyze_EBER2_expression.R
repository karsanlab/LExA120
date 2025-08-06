###############################################################################
#
# Author         | Dr. Shujun Huang 
# Maintainer     | Joshua Bridgers (jbridgers@bcgsc.ca)
# PI (Lab)       | Dr. Aly Karsan (Karsan Lab)
#                |   https://www.bcgsc.ca/labs/karsan-lab
# 
# Date Created   | 2021-12-12
# Description    | Create ROC curves on EBV status. Plot EBER2 expression and 
#                | EBER-ISH status
#
# License        | MIT License (see https://opensource.org/licenses/MIT)
#
###############################################################################



# set up --------------------------------------------------------------------------------------

library(here)
library(readxl)
library(plotROC)
library(pROC)
library(ROCR)
library(caret)
library(tidyverse)



# load expression data  -----------------------------------------------------------------------

# raw and log2 transformed expression dataframes of bcc and uhn rounds 1-4
load(here::here("data", "lexa_relevant_Expr.rda"))
libraries <- read.csv(here::here("data", "libraryIDs.csv"), header = T, stringsAsFactors = F)



# load EBV phenotype data ---------------------------------------------------------------------

# The 8 Round 1-3 cHL cases
original_ebv_status <- read_excel(here::here("data", "ebv_and_bcl6", "HL_EBER_status_8_cases_Oct_2020.xlsx"), sheet = 1)
colnames(original_ebv_status) <- c("SourceName", "External_ID", "SampleType", "SampleSource", "TMA", "HRS", "Background")

# The 48 Round 5 cHL cases
lexa_ebv_status <- read_excel(here::here("data", "ebv_and_bcl6", "20210330_48rHL_RNA_submission.xlsx"), sheet = 1)



# process expression data ---------------------------------------------------------------------

# Eight cHL cases from Round 1-3 replicates
hls_expr <- bind_rows(bcc_round1_logExpr, bcc_round2_logExpr, bcc_round3_logExpr) %>%
  dplyr::select(SampleName, SampleType, UniqID, ResID, RoundID, EBER2) %>%
  filter(SampleType == "cHL")

hls_expr_mean <- hls_expr %>%
  group_by(ResID) %>%
  summarise(EBER2_mean = mean(EBER2)) %>%
  ungroup()

original_hls_expr <- hls_expr %>%
  left_join(hls_expr_mean, by = "ResID") %>%
  filter(RoundID == "BCC-R1") %>%
  mutate(RoundID = "BCC_R1-3") %>%
  dplyr::select(-EBER2) %>%
  dplyr::rename(EBER2 = EBER2_mean)

# Round 5 cHL cases
lexa_hl_expr <- additional_hl_logExpr %>%
  dplyr::select(SampleName, SampleType, UniqID, ResID, RoundID, EBER2)



# merge EBER2 expression + EBV status for analysis --------------------------------------------

# Round 1-3 cHL cases
original_data <- inner_join(
  original_hls_expr,
  original_ebv_status %>% dplyr::select(External_ID, SourceName, HRS, Background),
  by = c("UniqID" = "SourceName")
) %>%
  mutate(Background = ifelse(Background != "NEG", "POS", "NEG")) %>%
  mutate(Group = case_when(
    HRS == "NEG" & Background == "NEG" ~ "3",
    HRS == "NEG" & Background == "POS" ~ "2",
    HRS == "POS" & Background == "NEG" ~ "1"
  )) %>%
  arrange(EBER2)

# Round 5 cHL cases
lexa_data <- inner_join(
  lexa_hl_expr,
  lexa_ebv_status %>%
    dplyr::select(
      External_ID = "HL#",
      Patient_ID,
      HRS = "EBER+ HRS cells",
      Background = "EBER+ Background cells",
      Group
    ),
  by = c("ResID" = "Patient_ID")
) %>%
  arrange(EBER2)

merged_data <- bind_rows(original_data, lexa_data)
merged_data %>%
  group_by(HRS, Background) %>%
  summarise(count = n())



# annotate and cleanup ------------------------------------------------------------------------

# Format the data
EBER2_df <- merged_data %>%
  mutate(Cells = case_when(
    HRS == "POS" | (HRS == "NEG" & Background == "NEG") ~ "Hodgkin-Reed Sternberg cell",
    HRS == "NEG" & Background == "POS" ~ "Background non-tumour cell"
  )) %>%
  mutate(EBER_ish_status = ifelse(HRS == "POS", "Positive", "Negative")) %>%
  mutate(Label = ifelse(EBER_ish_status == "Positive", 1, 0)) %>%
  mutate(
    Group = case_when(Group == "1 or 2" ~ "2", TRUE ~ Group),
    UniqID = gsub("_", "", UniqID)
  ) %>%
  mutate(Housekeeper_normalization_factor = ifelse(UniqID %in% c("AK1133", "AK1252", "AK1266", "AK1267", "AK1268"), "fail", "pass")) %>%
  mutate(Group_label = case_when(Group == 1 ~ "HRS Pos", Group == 2 ~ "Background Pos", Group == 3 ~ "HRS Neg")) %>%
  mutate(Group_label = factor(Group_label, levels = c("HRS Neg", "Background Pos", "HRS Pos"))) %>%
  mutate(Group = factor(Group, levels = c("3", "2", "1")))

# Filter out samples with potential degradation from Nanostring QC
EBER2_df <- EBER2_df %>%
  filter(Housekeeper_normalization_factor == "pass")

write.csv(EBER2_df, here::here("results", "EBER2_expression_51_cHL.csv"), row.names = TRUE)



# compare EBER2 expression between ISH status (three classes) ---------------------------------

# pROC to generate the ROC plot
multiclass.roc(EBER2_df$Group, EBER2_df$EBER2)

# Youden's J statistic on ROC comparisons
coords(
  roc(EBER2_df$Group, EBER2_df$EBER2, levels = c("1", "3"), plot = TRUE, col = "gold1"),
  x = "best",
  best.method = "youden"
)

coords(
  roc(EBER2_df$Group, EBER2_df$EBER2, levels = c("1", "2"), plot = TRUE, add = TRUE, col = "darkblue"),
  x = "best",
  best.method = "youden"
)


coords(
  roc(EBER2_df$Group, EBER2_df$EBER2, levels = c("2", "3"), plot = TRUE, add = TRUE, col = "darkred"),
  x = "best",
  best.method = "youden"
)

# Initial thresholds determined by Youden's J statistic and then rounded to the closest EBER2 expression 
thresholds <- c(11.56, 9.0)

legend(
  x = "bottomright",
  lty = 1,
  cex = 0.8,
  legend = c("Pos vs Neg", "Pos vs Background Pos", "Neg vs Background Pos"),
  col = c("gold1", "darkblue", "darkred")
)
text(
  x = 0.0,
  y = 0.13,
  cex = 0.8,
  pos = 3,
  labels = paste0("AUC = ", round(multiclass.roc(EBER2_df$Group, EBER2_df$EBER2)$auc, 3))
)

# ROCR to generate the ROC plot
# Pos vs Intermediate (1 vs 2)
temp1 <- EBER2_df %>%
  filter(Group != 3) %>%
  droplevels()
pred1 <- prediction(temp1$EBER2, temp1$Group, label.ordering = c(2, 1))
perf1 <- performance(pred1, "tpr", "fpr")
sapply(slotNames(pred1), function(x) slot(pred1, x))

# Pos vs Neg (1 vs 3)
temp2 <- EBER2_df %>%
  filter(Group != 2) %>%
  droplevels()
pred2 <- prediction(temp2$EBER2, temp2$Group, label.ordering = c(3, 1))
perf2 <- performance(pred2, "tpr", "fpr")

# Intermediate vs Neg (2 vs 3)
temp3 <- EBER2_df %>%
  filter(Group != 1) %>%
  droplevels()
pred3 <- prediction(temp3$EBER2, temp3$Group, label.ordering = c(3, 2))
perf3 <- performance(pred3, "tpr", "fpr")

# aucs 
auc.perf1 <- performance(pred1, measure = "auc")
auc.perf2 <- performance(pred2, measure = "auc")
auc.perf3 <- performance(pred3, measure = "auc")

aucs <- c(
  unlist(auc.perf1@y.values),
  unlist(auc.perf2@y.values),
  unlist(auc.perf3@y.values)
)
names(aucs) <- c("pos_vs_bg", "neg_vs_pos", "neg_vs_bg")

# Roc plot
png(
  here::here(
    "plots",
    str_glue("roc_for_three_groups.png")
  ),
  width = 180,
  height = 150,
  units = "mm",
  res = 200
)
plot(
  perf1,
  colorize = F,
  lwd = 1,
  col = "darkblue",
  xlab = "False positive rate (1-Specificity)",
  ylab = "True positive rate (Sensitivity)",
  main = "ROC Curve Assessment of EBER2 Expression and EBERish Status"
)
plot(perf2, add = TRUE, colorize = F, lwd = 1, col = "gold1")
plot(perf3, add = TRUE, colorize = F, lwd = 1, col = "darkred")
text(x = 0.80, y = 0.1, cex = 0.8, pos = 3, labels = paste0("AUC = ", round(mean(aucs), 3)))
abline(a = 0, b = 1, col = "gray", lty = "dashed")
legend(
  x = "bottomright",
  lty = 1,
  cex = 0.6,
  legend = c("Pos vs Neg (AUC = 0.997)", "Pos vs Background Pos (AUC = 0.952)", "Neg vs Background Pos (AUC = 0.898)"),
  col = c("gold1", "darkblue", "darkred")
)
dev.off()

plot(
  x = NA,
  y = NA,
  xlim = c(0, 1),
  ylim = c(0, 1),
  ylab = "True Positive Rate",
  xlab = "False Positive Rate",
  bty = "n"
)
lines(perf1@x.values[[1]], perf1@y.values[[1]], col = 1)
lines(perf2@x.values[[1]], perf2@y.values[[1]], col = 2)
lines(perf3@x.values[[1]], perf3@y.values[[1]], col = 3)

# cutoffs --> the “optimal” cut point -- > the point closest to the TPR of (1) and FPR of (0)
cost.perf1 <- performance(pred1, "cost")
pred1@cutoffs[[1]][which.min(cost.perf1@y.values[[1]])]
cost.perf2 <- performance(pred2, "cost")
pred2@cutoffs[[1]][which.min(cost.perf2@y.values[[1]])]
cost.perf3 <- performance(pred3, "cost")
pred3@cutoffs[[1]][which.min(cost.perf3@y.values[[1]])]

cutoffs <- c(
  pred1@cutoffs[[1]][which.min(cost.perf1@y.values[[1]])],
  pred2@cutoffs[[1]][which.min(cost.perf2@y.values[[1]])],
  pred3@cutoffs[[1]][which.min(cost.perf3@y.values[[1]])]
)
names(cutoffs) <- c("pos_vs_bg", "neg_vs_pos", "neg_vs_bg")

# 2 (control) vs 1 (case): we do not want 2 to be classified as 1, so no FPs
cost.perf1.1 <- performance(pred1, "cost", cost.fp = 4, cost.fn = 1)
pred1@cutoffs[[1]][which.min(cost.perf1.1@y.values[[1]])]

# 3 (control) vs 2 (case): we do not want 2 to be classified as 3, so no FNs
cost.perf3.1 <- performance(pred3, "cost", cost.fp = 1, cost.fn = 12)
pred3@cutoffs[[1]][which.min(cost.perf3.1@y.values[[1]])]

cost.perf3.1 <- performance(pred3, "cost")
pred3@cutoffs[[1]][which.min(cost.perf3.1@y.values[[1]])]

cost.perf3.1 <- performance(pred2, "cost")
pred2@cutoffs[[1]][which.min(cost.perf3.1@y.values[[1]])]

cutoffs <- c(
  pred1@cutoffs[[1]][which.min(cost.perf1.1@y.values[[1]])],
  NA,
  pred2@cutoffs[[1]][which.min(cost.perf3.1@y.values[[1]])]
)
names(cutoffs) <- c("pos_vs_bg", "neg_vs_pos", "neg_vs_bg")

