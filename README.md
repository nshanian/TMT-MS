# TMT-MS

## Exploratory and Differential Abundnance Analysis of TMT-MS data

This repository contains tools for exploratory and differential abundance analysis of TMT-MS data using a t-test and Linear and Nonlinear Mixed Effects Models (LNMM) using the `nlme` R package.

[Click here](https://htmlpreview.github.io/?https://github.com/nshanian/Documents/blob/main/TMT-MS-Proteomics.html) for the HTML version of this workflow with the output and the embedded plots.

Differential Abundance Analysis of TMT-MS data

The overall design is a genetically targeted proximity labeling experiment with dCas9 as the targeting moiety. The goal is to explore this dataset, to initially identify any proteins that are showing differences in abundance across the various conditions. Then lastly, to perform differential abundance analysis between different conditions. 

There are 5 conditions (p1, p2, p3, NoGp, NoGn) each with 3 biological replicates.There are 15 columns representing the different conditions and their replicates. Each column represents a TMT channel. The values listed are abundances or intensities. 

Targeted positive controls:

-   The three positive conditions are p1, p2, and p3, with 3 replicates each. 
-   Each “p” group represents separate single guide sgRNAs, targeting the promoter of FOXP2 just upstream of the transcription start site.

Untargeted negative controls:

-   NoG is the dCas9 on its own with no guide RNA as a negative control.
-   NoGp is a negative control that includes full labeling conditions.
-   NoGn excluded the biotin phenol needed to actually do the labeling and represents background endogenously biotinylated proteins.

First, a normal distribution of data is assumed, and log intensities are compared across conditions, with z scores determined for proteins showing log2 > 1 fold change. Then the data is modeled to a normal distribution and a t-test statistic is calculated with FDR-adjusted p-value to determine the significance of the observed fold changes. 

For differential analysis, the log2 ratios between the promoter targeting cases are compared, that is, averaged p1-p3 controls vs either the averaged negative control with labeling (NoGp) or no labeling (NoGn). Linear and Nonlinear Mixed Effects model is used to extract fold change and p-value for each protein using averaged control vs NoGp or NoGn conditions.  


