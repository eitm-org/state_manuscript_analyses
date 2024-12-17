# STATE_analyses
Tooling and notebooks for analyzing STATE data


Scripts to make manuscript plots are organized like this:

![STATE figs](https://github.com/user-attachments/assets/fa17e95d-353f-4f96-a7ea-ce51a9c78546)


Script                                     | Input                                       | Output
------------------------------------------ | ------------------------------------------- | --------------
scripts/state_manuscript_data_analysis.ipy |                                             |
scripts/all_paper_figs_data_cleaning.R     | output from state_manuscript_data_analysis, |csvs for plotting
scripts/all_paper_figs.Rmd                 | output from all_paper_figs_data_cleaning    | all figures for manuscript


# ClairS workflows

Nextflow modules and workflows that implements and adapts [ClairS](https://github.com/HKU-BAL/ClairS) (Zheng et. al., 2023) for extracting sSNVs with no normal tissue reference. 
