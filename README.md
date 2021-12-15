# CCS-filter

This repository includes the scripts used for the CCS filtering workflow described in the accompanying manuscript *"Improving Confidence in Lipidomic Annotations by Incorporating Empirical Ion Mobility Regression Analysis and Chemical Class Prediction" (Rose, et. al.)*

### Included folders:
#### Function scripts: 
* Functions built to assist in data clean up and analysis:
  - **append.R** appends CCS data to the feature list
  - **classy.R** includes functions to help with molecular classification
  - **blanksub.R** subtracts features that also appear in the blank feature list
* Compendium regression model scripts: both require a copy of the Compendium csv to run (included in folder)
  - **models_classes.R** 
  - **models_subclasses.R**
* **Filtering functions**: all filtering functions for both the broad (class specific) and fine (feature specific) filters are included in the **filter.R** script. This requires output from the Compendium regression model scripts to run.

#### Proof of Concept Murine Tissue Study: 
* **Analysis.Rmd** full analysis workbook for the data analysis reported in the manuscript, including data cleaning, molecular classification, filtering, and output steps
* **Raw data output files**: 
  - Progenesis output files for features with tentative identifications (**HD_Pos_IDs.csv** and **HD_Neg_IDs.csv**)
  - Progenesis output files for features in the "blank" data files (**blank_pos.csv** and **blank_neg.csv**)
  - Mass Profiler output for IM features with CCS values (**agil_pos.xls** and **agil_neg.xls**)
  - Complete cleaned dataset with classified IDs (**classified_ids.csv**)
* **Sample filter output**:
  - Analysis markdown report html file (**Analysis.html**)
  - Sample filter output data with higher confidence identifications (**sampleoutput_highconfids.csv**)


### To run analysis with included test data:
* Clone the repository
* Set the local repository folder as the working directory in R
* Knit the document or run in individual chunks to work step by step
* Output will appear in main folder