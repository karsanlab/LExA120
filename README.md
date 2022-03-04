# LExA120
The LExA120 Nanostring panel is a combination of the RHL30 and DLBCL90 NanoString assays. 

This script aims to take in the LExA expression RCC files (after nSolver QC), generate the uncorrected LPSs, 
correct the LPSs (calibration factors derived from LExA validation BCC Rounds 4-5), and classify the samples 
The future LExA analysis pipeline can be based on this script

Author: Shujun Huang
Date: 2021-12-12

## Overview of the script tasks
#### For DLBCL samples
- Step 1: Read in the LExA rcc file for a batch of given samples 
- Step 2: Generate the LPS scores (uncorrected)
  - Use the DLBCL90 r package to generate the COO, PMBL, and DHIT LSP scores
- Step 3: Correcte the LPS scores (Derived from BCC Rounds 4-5 data) 
  - for COO (DLBCL90 LPS ~ LExA LPS): y = 1.078x - 5.611
  - for PMBL (DLBCL90 LPS ~ LExA LPS): y = 1.01x + 5.559
  - for DHIT  (DLBCL90 LPS ~ LExA LPS): y = 1.038x + 2.401
- Step 4: Perform the classification using the corrected LPS scores
  - need to calculate the probability
  - For COO: ABC vs GCB
  - For PMBL: PMBLvs DLBCL
  - For DHIT: NEG vs POS

#### For cHL samples
- Step 1: Read in the LExA rcc file for a batch of given samples  (Let's use BCC Round 2 data as an example)
- Step 2: Generate the LPS scores (uncorrected)
  - Use the RHL30 r package
- Step 3: Correct the LPS scores (Derived from BCC Rounds 4-5 data) 
  - for cHL Risk  (RHL30 LPS ~ LExA LPS): y = 1.004x + 0.231
- Step 4: Perform the classification using the corrected LPS scores
  - high risk vs low risk (10.4 as the cutoff from the RHL30 paper)
