# LExA120
The LExA120 NanoString combined panel is an optimization of the RHL30 and DLBCL90 NanoString assays. 

## Features
* Reads in NanoString rcc files
* Generates uncorrected linear prediction scores (LPS)
* Calibrates LPS derived from LExA validation rounds 4-5 at BC Cancer site
* Classifies COO, PMBCL and DZsig subtype. Additionally classifies cHL risk
* Generates ROC curve for _EBER2_ expresion
  
## Information
* **Author**: Dr. Shujun Huang 
* **Maintainer**: Joshua Bridgers 
* **Prinicple Invesitigator (Lab)**: Dr. Aly Karsan ([Karsan Lab](https://www.bcgsc.ca/labs/karsan-lab))


## Project Structure
```
|
├─ cache # Temporary files generated my repo
├─ data # Data generated my repo
├─ plots # Plots generated my repo
├─ README.md
├─ reports # Contains R and R markdown scripts
├─ resources # Contains DLBCL90 R package
└─ results # Results generated my repo
```

## Citing LExA120
Peter J.B. Sabatini, Josh Bridgers, Shujun Huang, Tong Zhang, Clare Sheen, Tracy Stockley, Robert Kridel, Ian Bosdet, Marco A. Marra, Christian Steidl, David W. Scott, Aly Karsan. Validation of a Modular Gene Expression Assay for Risk Stratification and Subtyping Lymphomas. (In revision).

## License
MIT License (see https://opensource.org/licenses/MIT)

## Contact
Please direct questions to: jbridgers [at] bcgsc [dot] ca
