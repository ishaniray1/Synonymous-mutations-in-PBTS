library(tidyverse)
library(gridExtra)
library(scales)

# Load data
df <- final_file_with_all_rows_1107  # Your dataframe

# ============================================================
# PART 1 THRESHOLD TABLE
# ============================================================

#  thresholds to test (as decimals)
thresholds <- c(0, 0.001, 0.005, 0.01, 0.015, 0.02, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25)

# Function to calculate metrics at each threshold
calculate_metrics <- function(data, af_threshold) {
  # Filter mutations BELOW the threshold 
  filtered <- data %>%
    filter(is.na(gnomAD_4_AF) | gnomAD_4_AF < af_threshold)
  
  # Count unique mutations (by chr-pos-ref-alt)
  unique_muts <- filtered %>%
    distinct(CHROM, POS, REF, ALT) %>%
    nrow()
  
  metrics <- data.frame(
    threshold = af_threshold,
    threshold_pct = af_threshold * 100,
    n_total_mutations = nrow(filtered),
    n_unique_mutations = unique_muts,
    n_genes = n_distinct(filtered$Hugo_Symbol),
    n_patients = n_distinct(filtered$Tumor_Sample_Barcode)
  )
  
  return(metrics)
}

#  metrics for all thresholds
results <- map_dfr(thresholds, ~calculate_metrics(df, .x))

# Print summary table
library(knitr)
# Summary table
summary_table <- results %>%
  dplyr::select(threshold_pct, n_total_mutations, n_unique_mutations, 
                n_genes, n_patients) %>%
  rename(`AF Threshold (%)` = threshold_pct,
         `Total Mutations` = n_total_mutations,
         `Unique Mutations` = n_unique_mutations,
         `Genes` = n_genes,
         `Patients` = n_patients)

kable(summary_table, format = "pipe", digits = 2)

# rate of change
results <- results %>%
  mutate(
    delta_total = n_total_mutations - lag(n_total_mutations),
    delta_unique = n_unique_mutations - lag(n_unique_mutations),
    pct_change_total = (delta_total / lag(n_total_mutations)) * 100,
    pct_change_unique = (delta_unique / lag(n_unique_mutations)) * 100
  )

# Rate of change table
roc_table <- results %>%
  dplyr::select(threshold_pct, delta_total, delta_unique, 
                pct_change_total, pct_change_unique) %>%
  filter(!is.na(delta_total)) %>%
  rename(`AF Threshold (%)` = threshold_pct,
         `Δ Total Muts` = delta_total,
         `Δ Unique Muts` = delta_unique,
         `% Change Total` = pct_change_total,
         `% Change Unique` = pct_change_unique)

kable(roc_table, format = "pipe", digits = 2)


# ============================================================
# PART 2 SAMPLE PREVALENCE VS ALLELE FREQUENCY CORRELATION
# ============================================================


# how many samples each mutation appears in
mutation_prevalence <- df %>%
  filter(!is.na(gnomAD_4_AF)) %>%  # Only analyze mutations with AF data
  group_by(CHROM, POS, REF, ALT) %>%
  summarize(
    n_samples = n(),  # How many times this mutation appears
    mean_AF = mean(gnomAD_4_AF, na.rm = TRUE),
    max_AF = max(gnomAD_4_AF, na.rm = TRUE),
    gene = first(Hugo_Symbol),
    .groups = "drop"
  )

# Categorize by sample prevalence
mutation_prevalence <- mutation_prevalence %>%
  mutate(
    prevalence_category = case_when(
      n_samples == 1 ~ "1 sample",
      n_samples == 2 ~ "2 samples",
      n_samples >= 3 & n_samples <= 5 ~ "3-5 samples",
      n_samples >= 6 & n_samples <= 10 ~ "6-10 samples",
      n_samples > 10 ~ ">10 samples"
    ),
    prevalence_category = factor(prevalence_category, 
                                 levels = c("1 sample", "2 samples", 
                                            "3-5 samples", "6-10 samples", 
                                            ">10 samples"))
  )

# Summary statistics

prevalence_summary <- mutation_prevalence %>%
  group_by(prevalence_category) %>%
  summarize(
    n_mutations = n(),
    median_AF = median(mean_AF, na.rm = TRUE),
    mean_AF = mean(mean_AF, na.rm = TRUE),
    max_AF = max(mean_AF, na.rm = TRUE),
    pct_with_AF_gt_0.01 = sum(mean_AF > 0.01) / n() * 100,
    pct_with_AF_gt_0.05 = sum(mean_AF > 0.05) / n() * 100,
    .groups = "drop"
  )

print(prevalence_summary)

# Correlation test
cor_test <- cor.test(mutation_prevalence$n_samples, 
                     mutation_prevalence$mean_AF, 
                     method = "spearman")

cat(sprintf("\n\nSpearman correlation between sample prevalence and AF:\n"))
cat(sprintf("  rho = %.4f\n", cor_test$estimate))
cat(sprintf("  p-value = %.2e\n", cor_test$p.value))

if (cor_test$p.value < 0.001) {
  cat("  *** Highly significant correlation\n")
} else if (cor_test$p.value < 0.05) {
  cat("  ** Significant correlation\n")
} else {
  cat("  Not significant\n")
}

# ============================================================
# PART 3: VISUALIZATIONS
# ============================================================


p4 <- ggplot(results, aes(x = threshold_pct, y = mutations_per_patient)) +
  geom_line(size = 1.2, color = "purple") +
  geom_point(size = 3, color = "purple") +
  labs(title = "Mutations per Patient vs gnomAD AF Threshold",
       x = "gnomAD AF Threshold (%)",
       y = "Mutations per Patient") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))
p4

# ============================================================

# Plot 1 Threshold elbow curve
p1 <- ggplot(results, aes(x = threshold_pct, y = n_total_mutations)) +
  geom_line(size = 1.2, color = "steelblue") +
  geom_point(size = 3, color = "steelblue") +
  labs(title = "Total Mutations v gnomAD 4 AF Threshold",
       x = "gnomAD 4 AF Threshold (%)",
       y = "Number of Total Mutations") +
  scale_x_log10()+
  scale_y_continuous(labels = comma) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))
p1

p2 <- ggplot(results, aes(x = threshold_pct, y = n_unique_mutations)) +
  geom_line(size = 1.2, color = "darkgreen") +
  geom_point(size = 3, color = "darkgreen") +
  labs(title = "Unique Mutations v gnomAD 4 AF Threshold",
       x = "gnomAD 4 AF Threshold (%)",
       y = "Number of Unique Mutations") +
  scale_x_log10()+
  scale_y_continuous(labels = comma) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))
p2

# Plot 2 Sample prevalence vs AF
p3 <- ggplot(mutation_prevalence, aes(x = n_samples, y = mean_AF)) +
  geom_point(alpha = 0.3, size = 1.5, color = "coral") +
  geom_smooth(method = "lm", color = "darkred", se = TRUE) +
  scale_y_log10(labels = percent_format(accuracy = 0.01)) +
  scale_x_log10() +
  labs(title = "Sample Prevalence v gnomAD 4 AF",
       subtitle = sprintf("Spearman rho = %.3f, p = %.2e", 
                          cor_test$estimate, cor_test$p.value),
       x = "Number of Samples (log scale)",
       y = "gnomAD 4 AF (log scale)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

# Plot 3 Boxplot by prevalence 
p4 <- ggplot(mutation_prevalence, aes(x = prevalence_category, y = mean_AF, 
                                      fill = prevalence_category)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0.3) +
  scale_y_log10(labels = percent_format(accuracy = 0.01)) +
  labs(title = "AF Distribution by Sample Prevalence",
       x = "Sample Prevalence Category",
       y = "gnomAD 4 AF (log scale)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  scale_fill_brewer(palette = "Set2")

# Plot 4 AF distribution with threshold 
p5 <- ggplot(df %>% filter(!is.na(gnomAD_4_AF) & gnomAD_4_AF > 0), 
             aes(x = gnomAD_4_AF)) +
  geom_histogram(bins = 50, fill = "lightblue", color = "black") +
  scale_x_log10(labels = percent_format(accuracy = 0.01)) +
  labs(title = "Distribution of gnomAD 4 AF in Dataset",
       x = "gnomAD 4 AF (log scale)",
       y = "Count") +
  theme_minimal() +
  geom_vline(xintercept = c(0.001, 0.01, 0.05, 0.10), 
             linetype = "dashed", 
             color = c("red", "blue", "green", "purple"),
             size = 1) +
  annotate("text", x = 0.001, y = Inf, label = "0.1%", 
           vjust = 2, hjust = -0.1, color = "red", fontface = "bold") +
  annotate("text", x = 0.01, y = Inf, label = "1%", 
           vjust = 2, hjust = -0.1, color = "blue", fontface = "bold") +
  annotate("text", x = 0.05, y = Inf, label = "5%", 
           vjust = 2, hjust = -0.1, color = "green", fontface = "bold") +
  annotate("text", x = 0.10, y = Inf, label = "10%", 
           vjust = 2, hjust = -0.1, color = "purple", fontface = "bold") +
  theme(plot.title = element_text(face = "bold"))

#  all plots
grid.arrange(p1, p2, ncol = 2)
# Zoom in on first part of x-axis for Total Mutations
p1_zoom <- p1 +
  coord_cartesian(xlim = c(0, 1))  
p1_zoom
# Zoom in on first part of x-axis for Unique Mutations
p2_zoom <- p2 +
  coord_cartesian(xlim = c(0, 1)) 
p2_zoom
grid.arrange(p3, p4, ncol = 2)
print(p5)

# ============================================================


# Recurrent high-AF mutations 
recurrent_high_af <- mutation_prevalence %>%
  filter(n_samples >= 5 & mean_AF > 0.01) %>%
  arrange(desc(n_samples), desc(mean_AF)) %>%
  dplyr::select(gene, CHROM, POS, REF, ALT, n_samples, mean_AF)

cat(sprintf("Found %d recurrent mutations (≥5 samples) with AF > 1%%:\n", 
            nrow(recurrent_high_af)))
print(head(recurrent_high_af, 20))

# non recurrent mutations with high AF 
private_high_af <- mutation_prevalence %>%
  filter(n_samples == 1 & mean_AF > 0.05) %>%
  arrange(desc(mean_AF)) %>%
  dplyr::select(gene, CHROM, POS, REF, ALT, n_samples, mean_AF)

cat(sprintf("\n\nFound %d private mutations with AF > 5%% (suspicious):\n", 
            nrow(private_high_af)))
print(head(private_high_af, 20))

# Export results??????
write.csv(results, "gnomad4_threshold_analysis.csv", row.names = FALSE)
write.csv(mutation_prevalence, "mutation_prevalence_analysis.csv", row.names = FALSE)
write.csv(recurrent_high_af, "recurrent_high_af_mutations.csv", row.names = FALSE)


