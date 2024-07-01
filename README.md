# TMT-MS

## Exploratory and Differential Protein Abundnance Analysis of TMT-MS Data 

This repository contains tools for exploratory and differential abundance analysis of TMT-MS data using a t-test and Linear and Nonlinear Mixed Effects Models (LNMM) using the `nlme` R package.

[Click here](https://htmlpreview.github.io/?https://github.com/nshanian/Documents/blob/main/TMT-MS-Proteomics.html) for the HTML version of this workflow with the output and the embedded plots.

The dataset found in this gDrive folder (named ‘proteomicsSampleData.csv’) contains log2 peptide intensities (column "Log_Intensity") from an isobaric proteomics experiment designed to study changes over time after human cells are subjected to a particular type of stress.  There are 2 cell lines (Stress and Control), and between 1-3 replicates were sampled from each cell line at 6 timepoints.  Multiplexed samples were processed in four batches.  

The data file is organized into the following columns:

●	Gene: A human-readable gene name associated with each protein

●	ProteinID: A unique protein identifier associated with each peptide

●	Peptide: A string corresponding to the different peptides

●	Batch:  An identifier with four levels representing different batches. Each batch was analyzed on a separate, non-consecutive day.

●	Time: Time point in the study in which the sample was extracted

●	CellType: A categorical variable describing weather Stress or Control cells were used

●	Replicate: A categorical variable indicating replicate measurements of the same peptide within a single mass spec run

●	Log_Intensity: The log2-transformed intensity measured by the mass spec

The task is to explore this data set and identify any proteins that might be changing in abundance through time and as a result of the stress induction.

This workflow will import the dataset, perform preprocessing and differential abundance analysis. First, a normal or normal-like distribution is assumed and a t-test is performed. Proteins are filtered based on an FDR-adjusted p-value threshold (< 0.05) and a log2 fold change threshold (> 1). Then Linear and Nonlinear Mixed Effects Models (LNMM) approach is taken to model the data and extract fold change information and the corresponding p-values. Initial results are reported and each model is evaluated by running diagnostics and comparing the residuals, or the difference between the predicted value and the actual value (i.e. the 'error' in the predicted value).

