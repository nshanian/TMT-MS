# This workflow will import the dataset, perform preprocessing and differential abundance analysis of TMT-MS proteomics dataset.

setwd("~/Desktop/TMT-MS")

# Load required libraries
library(preprocessCore)   # For data preprocessing
library(ggplot2)          # For data visualization
library(dplyr)            # For data manipulation
library(reshape2)         # For restructuring and aggregating data
library(stats)            # For statistical analysis
library(multcomp)         # For multiple comparisons of k groups in general linear models 
library(EnhancedVolcano)  # For visualization of differential analysis
library(nlme)             # For differential analysis using Linear and Nonlinear Mixed Effects Models (NLME)


# Load the dataset
data <- read.csv("~/Desktop/proteomicsSampleData.csv")

# Display the structure of the dataset
str(data)

# Data pre-processing
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


# Visualization: Create a plot of log2 intensity vs. time, colored by CellType
ggplot(data, aes(x = Time, y = Log_Intensity, color = CellType)) +
  geom_point(alpha = 0.5) +  # Set transparency for better visualization of overlapping points
  scale_color_manual(values = c("blue", "red")) + 
  labs(x = "Time", y = "Log2 Intensity", color = "Cell Type") +
  ggtitle("Peptide Intensity vs. Time") +
  facet_wrap(~CellType)


# Visualization: Create a boxplot of peptide intensity by Time and Cell Type
ggplot(data, aes(x = factor(Time), y = Log_Intensity, fill = CellType)) +
  geom_boxplot(position = "dodge", color = "black") +
  scale_fill_manual(values = c("blue", "red")) + 
  labs(x = "Time", y = "Log2 Intensity", fill = "Cell Type") +
  ggtitle("Box Plot of Peptide Intensity by Time and Cell Type") +
  theme_minimal()


# Visualization: Create a violin plot of peptide intensity by Time and Cell Type
ggplot(data, aes(x = factor(Time), y = Log_Intensity, fill = CellType)) +
  geom_violin(position = "dodge") +
  scale_fill_manual(values = c("blue", "red")) + 
  labs(x = "Time", y = "Log2 Intensity", fill = "Cell Type") +
  ggtitle("Violin Plot of Peptide Intensity by Time and Cell Type") +
  theme_minimal()


# Summarization: Calculate mean log2 intensity for each protein
protein_summary <- data %>%
  group_by(Gene, Time, CellType) %>%
  summarize(mean_log2_intensity = mean(Log_Intensity))


# Filtering: Identify proteins with > 1 fold changes in abundance 
filtered_proteins <- protein_summary %>%
  group_by(Gene, CellType) %>%
  summarize(max_intensity_change = max(mean_log2_intensity) - min(mean_log2_intensity)) %>%
  filter(max_intensity_change > 1)  # Set threshold based on desired fold change or significance level
head(filtered_proteins)


# Additional steps for sorting the data
sorted_data <- arrange(data, CellType, Batch, Replicate)


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


# t-test model diagnostics
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


# Let's visualize the results using a volcano plot as we did following differential analysis using a t-test.


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


# Let's perform a more statistically robust differential analysis by fitting our data using Linear and Nonlinear Mixed Effects Models


# differential analysis using LNMM for proteins showing max intensity changes}
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


# To calculate log-fold changes and p-values for all proteins in Stress vs Control conditions at each time point separately, we'll need to fit separate models for each time point or conduct appropriate contrasts within a single model. Here's how we can approach this:
  
#   Fit separate linear mixed models for each time point.
#   Extract coefficients and p-values for the comparison between Stress vs Control conditions at each time point.
#   Calculate log-fold changes based on the coefficients.
#   Adjust p-values for multiple testing if necessary.

# differential analysis between two cell types at each time point}
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


# combine data frames
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


# To test the assumptions of the linear mixed model, we can perform various diagnostic tests.

# linear mixed model diagnostics}
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


