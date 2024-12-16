---
title: "Protein Differential Abundance Analysis"
output: html_document
date: "2024-06-19"
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
knitr::opts_knit$set(root.dir = '~/Desktop/')
################################################################################
# note that outside of an Rmd code chunk, use `setwd()` to set the working directory in R
```

As in a regular R script in RStudio, a single line of code can be run with Command-Enter (Mac OS) or Ctrl-Enter (Windows). Whole chunks of code can be run with Command(/Ctrl) + Shift + Enter **or** by clicking the green "\>" button in the top-right corner of the chunk. Alternatively, these options can also be implemented by selecting lines of code and choosing the desired option from the drop-down menu in the "Run" tab, in the top-right corner of the of Source section of RStudio.

This workflow will import the dataset, perform preprocessing and differential abundance analysis. First, a normal or normal-like distribution is assumed and a t-test is performed. Proteins are filtered based on an FDR-adjusted p-value threshold (< 0.05) and a log2 fold change threshold (> 1). Then Linear and Nonlinear Mixed Effects Models (LNMM) approach is taken to model the data and extract fold change information and the corresponding p-values. Initial results are reported and each model is evaluated by running diagnostics and comparing the residuals, or the difference between the predicted value and the actual value (i.e. the 'error' in the predicted value). 

Here is the workflow outline: 

   1. Aggregation: Mean log2 intensity for each protein in each condition (Control and Stress) is calculated across replicates.
   2. Filtering: The data is separated into control_data and stress_data for easier manipulation.
   3. Differential Analysis: 
       - T-test: For each protein, a t-test is performed to compare the mean log2 intensities between Stress and Control conditions.
       - Linear and Nonlinear Mixed Effects Models: A more statistically robust approach is taken to account for deviations from normality.
   4. Multiple Testing Correction: p-values are adjusted using the FDR method to control the false discovery rate.        
   5. Model Diagnostics: Each model used is evaluated to assess the validity of the extracted coefficients like log2 fold changes. 

While initially, using a t-test to determine proteins undergoing changes in abundance is useful, given deviations from normality often seen in proteomics experiments, the structure of the data and the objectives of the analysis, using LNMM for differential abundance analysis seems more appropriate. The `nlme` package allows for flexible modeling of both fixed effects (e.g., experimental conditions like Time and CellType) and random effects (e.g., Batch/Replicate). This flexibility can be beneficial for capturing the hierarchical structure of the data and adjusting for nested experimental designs.

To comment or uncomment a series of lines, highlight the lines and use Command(/Ctrl) + Shift + C.

R packages `preprocessCore`, `ggplot2`, `dplyr`, `reshape2`, `stats`, `multcomp`, `EnhancedVolcano` as well as `nlme` are required for this workflow.They can be installed using the following commands:

```{r install packages, eval=F}
# install.packages("preprocessCore")
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("reshape2")
# install.packages("stats")
# install.packages("multcomp")
# install.packages("EnhancedVolcano")
# install.packages("nlme")

# If install.packages("packagename") command fails in newer versions of R uncomment and run the commands:
# if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install("packagename")
```

Once installed, load the packages.

```{r load libraries}
# Load required libraries
library(preprocessCore)   # For data preprocessing
library(ggplot2)          # For data visualization
library(dplyr)            # For data manipulation
library(reshape2)         # For restructuring and aggregating data
library(stats)            # For statistical analysis
library(multcomp)         # For multiple comparisons of k groups in general linear models 
library(EnhancedVolcano)  # For visualization of differential analysis
library(nlme)             # For differential analysis using Linear and Nonlinear Mixed Effects Models (LNMM)
```

```{r read log intensities}
# Load the dataset
data <- read.csv("~/Desktop/proteomicsSampleData.csv")

# Display the structure of the dataset
str(data)
```

```{r data pre-processing}
# Sort the data by CellType, Batch, and Replicate
data <- data[order(data$CellType, data$Batch, data$Replicate), ]

# Convert Log_Intensity column to a matrix with one column
log_intensity_matrix <- as.matrix(data$Log_Intensity)

# Normalize the matrix using preprocessCore
normalized_matrix <- normalize.quantiles(log_intensity_matrix)

# Convert the normalized matrix back to a vector 
normalized_log_intensity <- as.vector(normalized_matrix)

# Update the Log_Intensity column with the normalized values 
data$Log_Intensity <- normalized_log_intensity

# Verify the changes
head(data)
```

```{r visualize processed data}
# Visualization: Create a plot of log2 intensity vs. time, colored by CellType
ggplot(data, aes(x = Time, y = Log_Intensity, color = CellType)) +
  geom_point(alpha = 0.5) +  # Set transparency for better visualization of overlapping points
  scale_color_manual(values = c("blue", "red")) + 
  labs(x = "Time", y = "Log2 Intensity", color = "Cell Type") +
  ggtitle("Peptide Intensity vs. Time") +
  facet_wrap(~CellType)
```

```{r visualize processed data as box plots}
# Visualization: Create a box plot of peptide intensity by Time and Cell Type
ggplot(data, aes(x = factor(Time), y = Log_Intensity, fill = CellType)) +
  geom_boxplot(position = "dodge", color = "black") +
  scale_fill_manual(values = c("blue", "red")) + 
  labs(x = "Time", y = "Log2 Intensity", fill = "Cell Type") +
  ggtitle("Box Plot of Peptide Intensity by Time and Cell Type") +
  theme_minimal()
```

```{r visualize processed data as violin plots}
# Visualization: Create a violin plot of peptide intensity by Time and Cell Type
ggplot(data, aes(x = factor(Time), y = Log_Intensity, fill = CellType)) +
  geom_violin(position = "dodge") +
  scale_fill_manual(values = c("blue", "red")) + 
  labs(x = "Time", y = "Log2 Intensity", fill = "Cell Type") +
  ggtitle("Violin Plot of Peptide Intensity by Time and Cell Type") +
  theme_minimal()
```

```{r summarization}
# Summarization: Calculate mean log2 intensity for each protein
protein_summary <- data %>%
  group_by(Gene, Time, CellType) %>%
  summarize(mean_log2_intensity = mean(Log_Intensity))
```

```{r filering}
# Filtering: Identify proteins with > 1 fold changes in abundance 
filtered_proteins <- protein_summary %>%
  group_by(Gene, CellType) %>%
  summarize(max_intensity_change = max(mean_log2_intensity) - min(mean_log2_intensity)) %>%
  filter(max_intensity_change > 1)  # Set threshold based on desired fold change or significance level
head(filtered_proteins)
```

```{r sorting}
# Additional steps for sorting the data
sorted_data <- arrange(data, CellType, Batch, Replicate)
```

```{r}
# Load necessary libraries again 
library(dplyr)
library(reshape2)
library(multcomp)
library(stats)

# Step 1: Aggregate log2 peptide intensities at the protein level
protein_data <- sorted_data %>%
  group_by(Gene, ProteinID, CellType, Time, Replicate) %>%
  summarize(mean_Log_Intensity = mean(Log_Intensity, na.rm = TRUE))

# Step 2: Prepare data for t-tests
# Filter data for Control and Stress
control_data <- protein_data %>% filter(CellType == "Control")
stress_data <- protein_data %>% filter(CellType == "Stress")

# Step 3: Perform t-tests for each protein
results <- protein_data %>%
  group_by(Gene, ProteinID) %>%
  summarize(
    log2_fold_change = mean(mean_Log_Intensity[CellType == "Stress"], na.rm = TRUE) -
                       mean(mean_Log_Intensity[CellType == "Control"], na.rm = TRUE),
    p_value = t.test(
      mean_Log_Intensity[CellType == "Stress"],
      mean_Log_Intensity[CellType == "Control"]
    )$p.value
  ) %>%
  ungroup()

# Step 4: Adjust p-values using the FDR method
results <- results %>%
  mutate(FDR_adj_p_value = p.adjust(p_value, method = "fdr"))

# Step 5: Filter significant results based on FDR-adjusted p-value and log2 fold change threshold
significant_proteins <- results %>%
  filter(FDR_adj_p_value < 0.05 & abs(log2_fold_change) > 1)

# Print all proteins
# print(results)

# Print significant proteins
print(significant_proteins)

# Save results to CSV
# write.csv(results, 'results_diff_t-test.csv', row.names = FALSE)
# write.csv(significant_proteins, 'results_sig_FC_t-test.csv', row.names = FALSE)
```

There are 9 proteins showing negative fold change when comparing Stress vs Control conditions that pass the FDR-adjusted p-value threshold (< 0.05) and a log2 fold change threshold ( > 1). Negative fold change indicates that these proteins underwent a decrease in abundance going from the control to the stress condition. 

Let's visualize the results using a volcano plot and `EnhancedVolcano` package.

```{r}
# Load the EnhancedVolcano library
library(EnhancedVolcano)

# Create the volcano plot
EnhancedVolcano(
  results,
  lab = results$Gene,  # Labels for the points
  x = 'log2_fold_change',  # Log2 fold changes
  y = 'FDR_adj_p_value',  # FDR-adjusted p-values
  title = 'Volcano Plot Stress vs Control',
  xlab = bquote(~Log[2]~ 'fold change'),
  ylab = bquote(~-Log[10]~ 'FDR'),
  pCutoff = 0.05,  # P-value cutoff
  FCcutoff = 1,  # Fold change cutoff
  pointSize = 3.0,
  labSize = 4.0,
  col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
  colAlpha = 1,
  legendLabels = NULL,  # Remove legend
  legendPos = "none",   # Remove legend position
  legendLabSize = 0,     # Remove legend label size
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  gridlines.major = FALSE,
  gridlines.minor = FALSE
)
```

Next, let's perform some diagnostic tests on our model to see how the residuals are distributed. If they appear to be normally distributed based on the Q-Q plots and Shapiro-Wilk test, the assumptions of the t-test would be reasonably satisfied. If not, we might need to consider using non-parametric tests or transformations to better meet the assumptions of the t-test. 

```{r t-test model diagnostics}
# Load necessary libraries
library(dplyr)
library(reshape2)
library(ggplot2)

# Step 1: Aggregate log2 peptide intensities at the protein level
protein_data <- sorted_data %>%
  group_by(Gene, ProteinID, CellType, Time, Replicate) %>%
  summarize(mean_Log_Intensity = mean(Log_Intensity, na.rm = TRUE))

# Prepare data for residuals calculation
control_data <- protein_data %>% filter(CellType == "Control")
stress_data <- protein_data %>% filter(CellType == "Stress")

# Calculate residuals for each protein
residuals_data <- protein_data %>%
  group_by(Gene, ProteinID) %>%
  mutate(
    mean_control_intensity = mean(mean_Log_Intensity[CellType == "Control"], na.rm = TRUE),
    residuals = mean_Log_Intensity - mean_control_intensity
  ) %>%
  ungroup()

# Step 2: Generate Q-Q plots and perform Shapiro-Wilk test for normality on residuals
residuals_control <- residuals_data %>% filter(CellType == "Control") %>% pull(residuals)
residuals_stress <- residuals_data %>% filter(CellType == "Stress") %>% pull(residuals)

# Generate Q-Q plot for Control residuals
qqnorm(residuals_control)
qqline(residuals_control, col = "red")

# Generate Q-Q plot for Stress residuals
qqnorm(residuals_stress)
qqline(residuals_stress, col = "red")

# Perform Shapiro-Wilk test for normality on residuals
shapiro_control <- shapiro.test(residuals_control)
shapiro_stress <- shapiro.test(residuals_stress)

# Print Shapiro-Wilk test results
print(shapiro_control)
print(shapiro_stress)

# Alternatively, use ggplot2 for Q-Q plots
ggplot(residuals_data, aes(sample = residuals)) +
  stat_qq() +
  stat_qq_line(colour = "red") +
  facet_wrap(~ CellType, scales = "free") +
  theme_minimal() +
  labs(title = "Q-Q Plots of Residuals for Control and Stress Conditions",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles")
```

In the context of the Shapiro-Wilk test, the null hypothesis is that the data is normally distributed. The alternative hypothesis is that the data is not normally distributed. If the p-value is greater than 0.05, we fail to reject the null hypothesis, suggesting that the data does not significantly deviate from normality (i.e., it is reasonable to assume the data is normally distributed). If the p-value is less than 0.05, we reject the null hypothesis, indicating that the data significantly deviates from normality. 

Our W Values: Around 0.96 for both control and stress residuals. Our P-Values: Very small (5.541e-12 for control and 4.854e-12 for stress).

Despite the W values being relatively close to 1, the very small p-values suggest that the deviations from normality are statistically significant. This means that while the data might appear somewhat normally distributed based on the W statistic, the statistical test strongly suggests that there is some deviation from a perfect normal distribution.

Let's use a non-parametric tests that do not assume normality, such as the Wilcoxon rank-sum test to determine the significance of observed fold changes.

```{r}
# Perform Wilcoxon rank-sum test for each protein
results <- protein_data %>%
  group_by(Gene, ProteinID) %>%
  summarize(
    log2_fold_change = mean(mean_Log_Intensity[CellType == "Stress"], na.rm = TRUE) -
                       mean(mean_Log_Intensity[CellType == "Control"], na.rm = TRUE),
    p_value = wilcox.test(
      mean_Log_Intensity[CellType == "Stress"],
      mean_Log_Intensity[CellType == "Control"]
    )$p.value
  ) %>%
  ungroup()

# Adjust p-values using the FDR method
results <- results %>%
  mutate(fdr_adjusted_p_value = p.adjust(p_value, method = "fdr"))

# Filter significant results based on FDR-adjusted p-value and log2 fold change threshold
significant_proteins <- results %>%
  filter(fdr_adjusted_p_value < 0.05 & abs(log2_fold_change) > 1)

# Print significant proteins
print(significant_proteins)
```

Let's visualize the results using a volcano plot as we did following differential analysis using a t-test.

```{r}
# Create the volcano plot
EnhancedVolcano(
  results,
  lab = results$Gene,  # Labels for the points
  x = 'log2_fold_change',  # Log2 fold changes
  y = 'fdr_adjusted_p_value',  # FDR-adjusted p-values
  title = 'Volcano Plot Stress vs Control',
  xlab = bquote(~Log[2]~ 'fold change'),
  ylab = bquote(~-Log[10]~ 'FDR'),
  pCutoff = 0.05,  # P-value cutoff
  FCcutoff = 1,  # Fold change cutoff
  pointSize = 3.0,
  labSize = 4.0,
  col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
  colAlpha = 1,
  legendLabels = NULL,  # Remove legend
  legendPos = "none",   # Remove legend position
  legendLabSize = 0,     # Remove legend label size
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'black',
  gridlines.major = FALSE,
  gridlines.minor = FALSE
)
```

This is more or less what we saw using the parametric model fitting and the t-test statistic determination. 

Let's perform a more statistically robust differential analysis by fitting our data using Linear and Nonlinear Mixed Effects Models (LNMM) approach on the 6 proteins showing max intensity changes in Stress vs Control to determine whether the observed changes are significant. We can then extend the approach to all for all proteins in Stress vs Control conditions at each time point separately.

```{r differential analysis using LNMM for proteins showing max intensity changes}
# Load required library
library(multcomp)
library(nlme)

# Convert CellType to factor
data$CellType <- as.factor(data$CellType)

# Filter data for the selected proteins
selected_genes <- c("LRIG3", "POLR3H", "FGL2", "CACYBP", "NRGN", "COA7")
filtered_data <- data[data$Gene %in% selected_genes, ]

# Calculate log fold change
fold_change <- tapply(filtered_data$Log_Intensity, filtered_data$Gene, function(x) max(x) - min(x))

# Perform statistical tests
results <- sapply(selected_genes, function(gene) {
  # Subset data for the current gene
  gene_data <- filtered_data[filtered_data$Gene == gene, ]
  
  # Fit linear mixed model
  model <- lme(Log_Intensity ~ CellType, random = ~1 | Batch/Replicate, data = gene_data)
  
  # Perform hypothesis tests for CellType
  test_result <- glht(model, linfct = mcp(CellType = "Tukey"))
  
  # Extract the p-value for the contrast of interest (Stress vs Control)
  p_value <- summary(test_result)$test$pvalues[2]
  
  return(p_value)
})

# Adjust p-values using FDR correction
adjusted_p_values <- p.adjust(results, method = "fdr")

# Combine fold change and adjusted p-values into a data frame
results_df <- data.frame(Gene = names(results), Log_Fold_Change = fold_change, P_Value = adjusted_p_values, row.names = NULL)

# Print results
print(results_df)
```

To calculate log-fold changes and p-values for all proteins in Stress vs Control conditions at each time point separately, we'll need to fit separate models for each time point or conduct appropriate contrasts within a single model. Here's how we can approach this:

-   Fit separate linear mixed models for each time point.
-   Extract coefficients and p-values for the comparison between Stress vs Control conditions at each time point.
-   Calculate log-fold changes based on the coefficients.
-   Adjust p-values for multiple testing if necessary.

```{r differential analysis between two cell types at each time point}
library(dplyr)
library(nlme)

# Create an empty data frame to store results
results_list <- list()

# Fit separate linear mixed models for each time point except 29
for (time_point in setdiff(unique(data$Time), 29)) {
  # Subset the data for the current time point
  subset_data <- filter(data, Time == time_point)
  
  # Check if both levels of CellType are present
  if (all(levels(subset_data$CellType) %in% c("Control", "Stress"))) {
    # Fit the linear mixed model
    model <- lme(Log_Intensity ~ CellType, random = ~1 | Batch/Replicate, data = subset_data)
   
    # Extract coefficients and p-values
    coefficients <- fixef(model)
    summary_info <- summary(model)
    
    # Get p-value for CellTypeStress coefficient
    p_value <- summary_info$tTable["CellTypeStress", "p-value"]
    
    # Check if p_value is not NA
    if (!is.na(p_value)) {
      # Store results in a list
      results_list[[as.character(time_point)]] <- c(coefficients["CellTypeStress"], p_value)
    } else {
      print(paste("Skipping time point", time_point, "due to missing p-value for CellTypeStress"))
    }
  } else {
    # Print a message if one of the levels is missing
    print(paste("Skipping time point", time_point, "due to missing levels of CellType"))
  }
}

# Convert the results to a data frame
results_df <- do.call(rbind.data.frame, results_list)
colnames(results_df) <- c("Log_Fold_Change", "P_Value")

# Add the Time column
results_df$Time <- as.numeric(rownames(results_df))

# Print the results
print(results_df)
```

```{r combine data frames}
library(stats)
library(ggplot2)    

# Join the results dataframe with the original data to obtain gene names
results_with_genes <- inner_join(results_df, data, by = "Time")

# Calculate FDR-adjusted p-values
results_with_genes <- results_with_genes %>%
mutate(FDR = p.adjust(p_value, method = "fdr"))

print(results_with_genes)

# Create a volcano plot
# Define the limits for the x-axis and y-axis
x_limits <- c(-2, 2)  # Log2 fold change from -2 to 2
y_limits <- c(0, 20)  # -log10(FDR) from 0 to 20

# Volcano plot
volcano_plot <- ggplot(results_with_genes, aes(x = log_fold_changes, y = -log10(FDR), color = CellType)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  xlim(x_limits) +  # Set x-axis limits
  ylim(y_limits) +  # Set y-axis limits
  labs(x = "Log2 Fold Change (Stress vs Control)", y = "-log10(FDR-adjusted p-value)", color = "Cell Type", title = "Volcano Plot") +
  theme_minimal()

# Show the plot 
# print(volcano_plot)
```

In summary, while the observed max intensity changes and log2 fold changes in abundance for proteins "LRIG3", "POLR3H", "FGL2", "CACYBP", "NRGN", "COA7" as determined by t-test seem noteworthy, the assumption that the data has a normal or normal-like distribution doesn't really hold. The Shapiro-Wilk normality test indicates the data significantly deviates from normality. After using a more statistically robust model to fit the data the lack of statistical significance suggests caution in interpreting these findings as conclusive evidence of differential abundance between the Control and Stress conditions. 

Further experimentation or validation by correlating protein abundance data with mRNA expression data by RNA-seq may be necessary to confirm these observations.

And last but not least, while using a linear mixed model approach for differential abundance analysis in proteomics experiments is well suited for complex multidimensional data, its appropriateness depends on several factors. Here are some considerations:

- Linearity: LMMs assume a linear relationship between the response variable and the predictor variables. 
- Normality of Residuals: LMMs assume that the residuals (errors) are normally distributed. Violations of this assumption can affect the       validity of statistical inference.
- Homogeneity of Variance: LMMs assume homogeneity of variance, meaning that the variance of the residuals is constant across all levels of    the predictors. Heteroscedasticity can lead to biased parameter estimates and incorrect hypothesis testing.
- Independence of Observations: LMMs assume independence of observations within groups defined by the random effects. 

Model Diagnostics: It's essential to check the assumptions of the linear mixed model, such as normality, homoscedasticity, and independence of residuals. Residual diagnostics, such as residual plots and tests, can help evaluate whether these assumptions hold.

To test the assumptions of the linear mixed model, we can perform various diagnostic tests.

```{r linear mixed model diagnostics}
# Fit the linear mixed effects model
model <- lme(Log_Intensity ~ CellType, random = ~1 | Batch/Replicate, data = subset_data)

# Extract residuals and fitted values
resid <- residuals(model)
fitted <- fitted(model)

# Diagnostic plots
par(mfrow = c(2, 2))  # Set up a 2x2 grid for multiple plots

# Residual vs. Fitted Values Plot
plot(fitted, resid, xlab = "Fitted Values", ylab = "Residuals")
abline(h = 0, col = "red")  # Add a red line at y = 0

# Normal Q-Q Plot
qqnorm(resid)
qqline(resid)

# Scale-Location Plot
plot(sqrt(abs(resid)) ~ fitted, xlab = "Fitted Values", ylab = "Sqrt(|Residuals|)")
abline(h = 0, col = "red")  # Add a red line at y = 0

# Residuals vs. Index Plot (for checking for autocorrelation)
index <- 1:length(resid)
plot(index, resid, xlab = "Index", ylab = "Residuals")
abline(h = 0, col = "red")  # Add a red line at y = 0

# Reset the plotting layout
par(mfrow = c(1, 1))
```
The Q-Q plot suggests that the residuals are approximately normally distributed and that the assumption of normality for the residuals is met. However, the presence of a few outliers falling below the diagonal line in the -3 quantile could indicate some deviations from normality in the tails of the distribution. 

While the Q-Q plot indicates that the residuals are normally distributed, the Residuals vs Fitted Values plot suggests potential non-linearity. The groups of overlapping points forming vertical lines could indicate some systematic patterns or structure in the residuals that are not captured by the linear mixed model.

In such cases, further investigation is warranted to understand the underlying reasons for the observed patterns in the residuals, to revisit experimental and data acquisition considerations, and to determine whether any model adjustments, transformations or alternative models are necessary to improve model fit and validity.

