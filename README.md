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
├─ cache # Temporary files generated the repo
├─ data # Data generated the repo
├─ plots # Plots generated the repo
├─ README.md
├─ reports # Contains R and R markdown scripts
└─ results # Results generated the repo
```

## Citing LExA120
Sabatini PJB, Bridgers J, Huang S, Zhang T, Sheen C, Stockley T, Kridel R, Bosdet I, Marra MA, Steidl C, Scott DW, Karsan A. Validation of a Modular Gene Expression Assay for Risk Stratification and Subtyping Lymphomas. J Mol Diagn. 2025 Oct 10:S1525-1578(25)00239-9. doi: 10.1016/j.jmoldx.2025.09.004. PMID: 41077189.

## License
MIT License (see https://opensource.org/licenses/MIT)

## Contact
Please direct questions to: jbridgers [at] bcgsc [dot] ca
