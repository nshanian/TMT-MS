---
title: "Protein Differential Abundance Analysis"
output: html_document
date: "2024-06-11"
---

```{css, echo = FALSE}
pre, code {white-space:pre !important; overflow-x:auto}
```

```{r setup, include=FALSE}
# when you "knit" this file, do you want the resulting PDF to print the code in each chunk (TRUE = yes)?
knitr::opts_chunk$set(echo = TRUE)

################################################################################
# set your working directory
####CHANGE THIS TO THE APPROPRIATE PATH
knitr::opts_knit$set(root.dir = '~/Desktop/Comp_MS/')
################################################################################
# note that outside of an Rmd code chunk, use `setwd()` to set the working directory in R
```

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

As in a regular R script in RStudio, a single line of code can be run with Command-Enter (Mac OS) or Ctrl-Enter (Windows). Whole chunks of code can be run with Command(/Ctrl) + Shift + Enter **or** by clicking the green "\>" button in the top-right corner of the chunk. Alternatively, these options can also be implemented by selecting lines of code and choosing the desired option from the drop-down menu in the "Run" tab, in the top-right corner of the of Source section of RStudio.

This workflow will import the dataset, perform preprocessing and differential analysis using Linear and Nonlinear Mixed Effects Models approach. Given the structure of the data and the objectives of the analysis, using LNME for differential abundance analysis seems appropriate. `nlme` package allows for flexible modeling of both fixed effects. This flexibility can be beneficial for capturing the hierarchical structure of the data and adjusting for nested experimental designs.

To comment or uncomment a series of lines, highlight the lines and use Command(/Ctrl) + Shift + C.

R packages `preprocessCore`, `ggplot2`, `dplyr`, `reshape`, `multcomp` and `nlme` and `Enhanced Volcano` are required for this workflow.
They can be installed using the following commands:

```{r install packages, eval=F}
# install.packages("preprocessCore")
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("multcomp")
# install.packages("nlme")
# install.packages("reshape2")
# install.packages("EnchancedVolcano")

# If install.packages("packagename") command fails in newer versions of R uncomment and run the commands:
# if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("packagename")
```

Once installed, load the packages.

```{r load libraries}
# Load required libraries
library(preprocessCore) # For data preprocessing
library(reshape2)       # For data preprocessing
library(ggplot2)        # For data visualization
library(dplyr)          # For data manipulation
library(nlme)           # For differential analysis
library(multcomp)       # For multiple comparisons of k groups in general linear models 
```

```{r read counts}
# Load the dataset
data <- read.csv("~/Desktop/Comp_MS/proteomicsData.csv")

head(data)

# Display the structure of the dataset
str(data)
```

```{r}
# Load necessary libraries
library(preprocessCore)
library(reshape2)

# Step 1: Exclude the ProteinID column for transformation and normalization
# and the Pool column for initial processing
intensity_data <- data[, -1]
pool_column <- intensity_data$Pool
intensity_data <- intensity_data[, -ncol(intensity_data)]

# Step 2: Normalize each intensity value by the corresponding value in the "Pool" column
normalized_data <- sweep(intensity_data, 1, pool_column, FUN = "/")

# Step 3: Apply log2 transformation
normalized_data <- log2(normalized_data + 1)  # Adding 1 to avoid log2(0)

# Step 4: Add back the ProteinID column
normalized_data <- cbind(ProteinID = data$ProteinID, normalized_data)
colnames(normalized_data) <- colnames(data)[-ncol(data)]  # Ensure column names are consistent

# Display the first few rows of the normalized data
head(normalized_data)
```

```{r}
# Calculate average log2 intensities for conditions
average_log_intensity_control <- rowMeans(normalized_data[, grep("^p[1-3]_", names(normalized_data))], na.rm = TRUE)
average_log_intensity_NoGp <- rowMeans(normalized_data[, grep("^NoGp", names(normalized_data))], na.rm = TRUE)
average_log_intensity_NoGn <- rowMeans(normalized_data[, grep("^NoGn", names(normalized_data))], na.rm = TRUE)

# Calculate log2 fold changes
log2_fold_change_NoGp <- average_log_intensity_NoGp - average_log_intensity_control
log2_fold_change_NoGn <- average_log_intensity_NoGn - average_log_intensity_control

# Calculate Z-scores for each fold change
z_score_NoGp <- (log2_fold_change_NoGp - mean(log2_fold_change_NoGp, na.rm = TRUE)) / sd(log2_fold_change_NoGp, na.rm = TRUE)
z_score_NoGn <- (log2_fold_change_NoGn - mean(log2_fold_change_NoGn, na.rm = TRUE)) / sd(log2_fold_change_NoGn, na.rm = TRUE)

# Create data frame for results
fold_change_data <- data.frame(
  ProteinID = normalized_data$ProteinID,
  log2_fold_change_NoGp = log2_fold_change_NoGp,
  log2_fold_change_NoGn = log2_fold_change_NoGn,
  z_score_NoGp = z_score_NoGp,
  z_score_NoGn = z_score_NoGn
)

# Filter significant proteins with log2 fold change > 1
significant_proteins_NoGp <- fold_change_data[abs(fold_change_data$log2_fold_change_NoGp) > 1, ]
significant_proteins_NoGn <- fold_change_data[abs(fold_change_data$log2_fold_change_NoGn) > 1, ]

# Save the results to a CSV file
# write.csv(fold_change_data, 'differential_analysis_foldchange_zscores.csv', row.names = FALSE)
# write.csv(significant_proteins_NoGp, 'significant_proteins_NoGp_zscores.csv', row.names = FALSE)
# write.csv(significant_proteins_NoGn, 'significant_proteins_NoGn_zscores.csv', row.names = FALSE)

# Print significant proteins
head(significant_proteins_NoGp)
head(significant_proteins_NoGn)
```

```{r visualize processed data}
# Melt the data for visualization
melted_data <- melt(normalized_data, id.vars = "ProteinID")

# Plot the distribution of log2 intensities
ggplot(melted_data, aes(x = variable, y = value)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Distribution of Log2 Intensities", x = "Condition", y = "Log2 Intensity")
```

```{r}
library(dplyr)
library(tidyr)

# Reshape the data into a longer format
long_data <- normalized_data %>%
  pivot_longer(cols = -ProteinID, names_to = "Condition") %>%
  separate(Condition, into = c("Condition", "Replicate"), sep = "_")

# Remove the numerical suffixes from the column names
long_data$Condition <- gsub("\\_\\d+", "", long_data$Condition)

# Calculate the average for each condition
averaged_data <- long_data %>%
  group_by(ProteinID, Condition) %>%
  summarise(Average = mean(value))

# Reshape the data back to its original wide format
averaged_data <- pivot_wider(averaged_data, names_from = Condition, values_from = Average)

# View the first few rows of the averaged data
head(averaged_data)
```

```{r}
# Melt the data for visualization
melted_data <- melt(averaged_data, id.vars = "ProteinID")

# Add a column for the condition type (control or negative control)
melted_data$ConditionType <- ifelse(melted_data$variable %in% c("p1", "p2", "p3"), "Positive Control", "Negative Control")

# Plot the distribution of log2 intensities using boxplot
ggplot(melted_data, aes(x = variable, y = value, fill = ConditionType)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Positive Control" = "blue", "Negative Control" = "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Distribution of Log2 Intensities (Boxplot)", x = "Condition", y = "Log2 Intensity")

# Plot the distribution of log2 intensities using violin plot
ggplot(melted_data, aes(x = variable, y = value, fill = ConditionType)) +
  geom_violin() +
  scale_fill_manual(values = c("Positive Control" = "blue", "Negative Control" = "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Distribution of Log2 Intensities (Violin Plot)", x = "Condition", y = "Log2 Intensity")
```

```{r}
# Let's calculate fold change and z-scores for each protein's log2 fold change. The z-score will help identify how many 
# standard deviations away a protein's log2 fold change is from the mean log2 fold change of the entire dataset.
# We will calculate z-scores for both Control vs NoGp and Control vs NoGn comparisons.

# Load required libraries
library(data.table)  # for fread and fwrite

# Calculate mean log2 intensity for each protein
normalized_data$mean_intensity <- rowMeans(normalized_data[ , -1])

# Calculate the average log2 intensity for the control group (p1, p2, p3)
control_group_means <- rowMeans(normalized_data[, grep("^p[1-3]_", names(normalized_data))], na.rm = TRUE)

# Calculate the average log2 intensity for the NoGp group
average_NoGp <- rowMeans(normalized_data[, grep("^NoGp", names(normalized_data))], na.rm = TRUE)

# Calculate the average log2 intensity for the NoGn group
average_NoGn <- rowMeans(normalized_data[, grep("^NoGn", names(normalized_data))], na.rm = TRUE)

# Calculate log2 fold changes
log2_fold_change_NoGp <- log2(control_group_means / average_NoGp)
log2_fold_change_NoGn <- log2(control_group_means / average_NoGn)

# Create a data frame to store the results
fold_change_data <- data.frame(
  ProteinID = normalized_data$ProteinID,
  log2_fold_change_NoGp = log2_fold_change_NoGp,
  log2_fold_change_NoGn = log2_fold_change_NoGn
)

# Calculate Z-scores for each fold change
fold_change_data$z_score_NoGp <- (fold_change_data$log2_fold_change_NoGp - mean(fold_change_data$log2_fold_change_NoGp, na.rm = TRUE)) / sd(fold_change_data$log2_fold_change_NoGp, na.rm = TRUE)
fold_change_data$z_score_NoGn <- (fold_change_data$log2_fold_change_NoGn - mean(fold_change_data$log2_fold_change_NoGn, na.rm = TRUE)) / sd(fold_change_data$log2_fold_change_NoGn, na.rm = TRUE)

# Filter significant proteins with log2 fold change > 2
significant_proteins_NoGp <- fold_change_data[abs(fold_change_data$log2_fold_change_NoGp) > 1, ]
significant_proteins_NoGn <- fold_change_data[abs(fold_change_data$log2_fold_change_NoGn) > 1, ]

# Save the results to a CSV file
# write.csv(fold_change_data, 'foldchange_zscores.csv', row.names = FALSE)
# write.csv(significant_proteins_NoGp, 'NoGp_zscores.csv', row.names = FALSE)
# write.csv(significant_proteins_NoGn, 'NoGn_zscores.csv', row.names = FALSE)

# Print significant proteins
head(significant_proteins_NoGp)
head(significant_proteins_NoGn)
```

```{r}
# Load necessary libraries if not already loaded
library(reshape2)  # For melt function
library(multcomp)  # For p-value adjustment
library(stats)     # For t.test and p.adjust functions

# Melt the normalized data for visualization
model_data <- melt(normalized_data, id.vars = "ProteinID")
colnames(model_data) <- c("ProteinID", "Condition", "Log2_Intensity")

# Calculate the average log2 intensity for each condition
control_group <- c("p1_1", "p1_2", "p1_3", "p2_1", "p2_2", "p2_3", "p3_1", "p3_2", "p3_3")
average_log_intensity_control <- rowMeans(normalized_data[, control_group], na.rm = TRUE)
average_NoGp <- rowMeans(normalized_data[, grep("^NoGp", colnames(normalized_data))], na.rm = TRUE)
average_NoGn <- rowMeans(normalized_data[, grep("^NoGn", colnames(normalized_data))], na.rm = TRUE)

# Calculate log2 fold changes
log2_fold_change_NoGp <- average_NoGp - average_log_intensity_control
log2_fold_change_NoGn <- average_NoGn - average_log_intensity_control

# Perform t-tests to get p-values, handling NAs and ensuring numeric data
p_values_NoGp <- apply(normalized_data, 1, function(row) {
  control_values <- as.numeric(row[control_group])
  nogp_values <- as.numeric(row[grep("^NoGp", names(row))])
  control_values <- control_values[!is.na(control_values)]
  nogp_values <- nogp_values[!is.na(nogp_values)]
  if (length(control_values) > 1 && length(nogp_values) > 1 && !all(control_values == control_values[1]) && !all(nogp_values == nogp_values[1])) {
    t.test(control_values, nogp_values)$p.value
  } else {
    NA
  }
})

p_values_NoGn <- apply(normalized_data, 1, function(row) {
  control_values <- as.numeric(row[control_group])
  nogn_values <- as.numeric(row[grep("^NoGn", names(row))])
  control_values <- control_values[!is.na(control_values)]
  nogn_values <- nogn_values[!is.na(nogn_values)]
  if (length(control_values) > 1 && length(nogn_values) > 1 && !all(control_values == control_values[1]) && !all(nogn_values == nogn_values[1])) {
    t.test(control_values, nogn_values)$p.value
  } else {
    NA
  }
})

# Store results in a data frame
results_NoGp <- data.frame(ProteinID = normalized_data$ProteinID, log2_fold_change = log2_fold_change_NoGp, p_value = p_values_NoGp)
results_NoGn <- data.frame(ProteinID = normalized_data$ProteinID, log2_fold_change = log2_fold_change_NoGn, p_value = p_values_NoGn)

# Filter significant results
significant_proteins_NoGp <- results_NoGp[results_NoGp$log2_fold_change > 1 & results_NoGp$p_value < 0.05, ]
significant_proteins_NoGn <- results_NoGn[results_NoGn$log2_fold_change > 1 & results_NoGn$p_value < 0.05, ]

# Adjust p-values for multiple testing using FDR correction
results_NoGp$p_adjusted <- p.adjust(results_NoGp$p_value, method = "fdr")
results_NoGn$p_adjusted <- p.adjust(results_NoGn$p_value, method = "fdr")

# Filter significant results after FDR adjustment
significant_proteins_NoGp_fdr <- results_NoGp[results_NoGp$log2_fold_change > 1 & results_NoGp$p_adjusted < 0.05, ]
significant_proteins_NoGn_fdr <- results_NoGn[results_NoGn$log2_fold_change > 1 & results_NoGn$p_adjusted < 0.05, ]

# Save results to CSV
# write.csv(significant_proteins_NoGp_fdr, 'results_foldchange_NoGp_fdr.csv', row.names = FALSE)
# write.csv(significant_proteins_NoGn_fdr, 'results_foldchange_NoGn_fdr.csv', row.names = FALSE)

# Print significant proteins and their adjusted p-values
head(significant_proteins_NoGp_fdr)
head(significant_proteins_NoGn_fdr)
```


```{r}
# Visualize with MA plots
plot(average_log_intensity_control, log2_fold_change_NoGp, 
     xlab = "Average Log2 Intensity (Control)", ylab = "Log2 Fold Change (vs NoGp)",
     main = "MA Plot: Control vs NoGp",
     col = ifelse(abs(log2_fold_change_NoGp) >= 1 & p_values_NoGp < 0.05, "red", "black"))

plot(average_log_intensity_control, log2_fold_change_NoGn, 
     xlab = "Average Log2 Intensity (Control)", ylab = "Log2 Fold Change (vs NoGn)",
     main = "MA Plot: Control vs NoGn",
     col = ifelse(abs(log2_fold_change_NoGn) >= 1 & p_values_NoGn < 0.05, "red", "black"))
```

```{r}
library(nlme)
library(multcomp)

# Step 2: Define a function to perform the analysis for each protein
analyze_protein <- function(protein) {
  protein_data <- model_data[model_data$ProteinID == protein, ]
  
  # Ensure no NA, NaN, or Inf values in the Log2_Intensity column
  protein_data <- protein_data[complete.cases(protein_data$Log2_Intensity), ]
  
  # Check if there are enough data points to fit the model
  if (nrow(protein_data) < 2) {
    return(data.frame(ProteinID = protein, log2_fold_change_NoGp = NA, log2_fold_change_NoGn = NA, p_value_NoGp = NA, p_value_NoGn = NA))
  }
  
  tryCatch({
    # Fit linear mixed model
    model <- lme(Log2_Intensity ~ Condition, random = ~1 | ProteinID, data = protein_data)
    
    # Perform hypothesis tests for Condition
    test_result <- glht(model, linfct = mcp(Condition = "Tukey"))
    
    # Extract the p-value for the contrasts of interest (Control vs NoGp and Control vs NoGn)
    p_values <- summary(test_result)$test$pvalues
    log2_fold_change_NoGp <- fixef(model)["ConditionNoGp"]
    log2_fold_change_NoGn <- fixef(model)["ConditionNoGn"]
    
    data.frame(ProteinID = protein, log2_fold_change_NoGp = log2_fold_change_NoGp, log2_fold_change_NoGn = log2_fold_change_NoGn, p_value_NoGp = p_values["ConditionNoGp"], p_value_NoGn = p_values["ConditionNoGn"])
  }, error = function(e) {
    # Return NA values in case of an error
    data.frame(ProteinID = protein, log2_fold_change_NoGp = NA, log2_fold_change_NoGn = NA, p_value_NoGp = NA, p_value_NoGn = NA)
  })
}

# Step 3: Apply the Function to Each Protein

# Get unique proteins
unique_proteins <- unique(model_data$ProteinID)

# Apply the analysis function to each protein
results <- lapply(unique_proteins, analyze_protein)

# Convert results to a data frame
results_df <- do.call(rbind, lapply(results, as.data.frame))

# Display the results
head(results_df)

# Save the results to a CSV file
# write.csv(results_df, 'results_differential_analysis.csv', row.names = FALSE)
```

```{r}
# Adjust p-values using FDR correction
results_df$p_value_NoGp <- p.adjust(results_df$p_value_NoGp, method = "fdr")
results_df$p_value_NoGn <- p.adjust(results_df$p_value_NoGn, method = "fdr")

# Print the first few rows of the fold change data
head(results_df)

# Filter significant hits for NoGp
significant_hits_NoGp <- results_df[abs(results_df$log2_fold_change_NoGp) >= 1 & results_df$p_value_NoGp < 0.05, ]

# Filter significant hits for NoGn
significant_hits_NoGn <- results_df[abs(results_df$log2_fold_change_NoGn) >= 1 & results_df$p_value_NoGn < 0.05, ]

# Print the tables of significant hits

head(significant_hits_NoGp)
head(significant_hits_NoGn)

# write.csv(results_df, 'results_differential_analysis.csv', row.names = FALSE)
```




