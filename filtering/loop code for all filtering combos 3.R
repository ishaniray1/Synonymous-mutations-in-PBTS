# output directory
dir.create("filter_combinations_output_3", showWarnings = FALSE)

# filter functions
apply_gnomad_filter <- function(df, threshold = 0.005) {
  df %>% filter(as.numeric(gnomAD_4_AF) < threshold)
}

apply_vaf_filter <- function(df) {
  df %>% filter(fracN == 0)
}

apply_m6a_filter <- function(df) {
  df %>% filter(in_m6A == TRUE)
}

apply_merip_filter <- function(df) {
  df %>% filter(in_merip == TRUE)
}

apply_m6a_or_merip_filter <- function(df) {
  df %>% filter(in_m6A == TRUE | in_merip == TRUE)
}

# all filter combinations (2^4 = 16 combinations + 3 additional "OR" combinations)
filter_combinations <- list(
  list(name = "No_Filters", filters = c(), label = "No Filters Applied"),
  list(name = "gnomAD_only", filters = c("gnomad"), label = "gnomAD < 0.005"),
  list(name = "VAF_only", filters = c("vaf"), label = "VAF (fracN = 0)"),
  list(name = "m6A_only", filters = c("m6a"), label = "m6A Overlap"),
  list(name = "MeRIP_only", filters = c("merip"), label = "MeRIP Overlap"),
  list(name = "m6A_OR_MeRIP", filters = c("m6a_or_merip"), label = "m6A OR MeRIP"),
  list(name = "gnomAD_VAF", filters = c("gnomad", "vaf"), label = "gnomAD + VAF"),
  list(name = "gnomAD_m6A", filters = c("gnomad", "m6a"), label = "gnomAD + m6A"),
  list(name = "gnomAD_MeRIP", filters = c("gnomad", "merip"), label = "gnomAD + MeRIP"),
  list(name = "gnomAD_m6A_OR_MeRIP", filters = c("gnomad", "m6a_or_merip"), label = "gnomAD + (m6A OR MeRIP)"),
  list(name = "VAF_m6A", filters = c("vaf", "m6a"), label = "VAF + m6A"),
  list(name = "VAF_MeRIP", filters = c("vaf", "merip"), label = "VAF + MeRIP"),
  list(name = "VAF_m6A_OR_MeRIP", filters = c("vaf", "m6a_or_merip"), label = "VAF + (m6A OR MeRIP)"),
  list(name = "m6A_MeRIP", filters = c("m6a", "merip"), label = "m6A + MeRIP"),
  list(name = "gnomAD_VAF_m6A", filters = c("gnomad", "vaf", "m6a"), label = "gnomAD + VAF + m6A"),
  list(name = "gnomAD_VAF_MeRIP", filters = c("gnomad", "vaf", "merip"), label = "gnomAD + VAF + MeRIP"),
  list(name = "gnomAD_VAF_m6A_OR_MeRIP", filters = c("gnomad", "vaf", "m6a_or_merip"), label = "gnomAD + VAF + (m6A OR MeRIP)"),
  list(name = "gnomAD_m6A_MeRIP", filters = c("gnomad", "m6a", "merip"), label = "gnomAD + m6A + MeRIP"),
  list(name = "VAF_m6A_MeRIP", filters = c("vaf", "m6a", "merip"), label = "VAF + m6A + MeRIP"),
  list(name = "All_Four", filters = c("gnomad", "vaf", "m6a", "merip"), label = "gnomAD + VAF + m6A + MeRIP")
)

# Storage for all results
all_results <- list()

# Loop through each filter combination
for (combo in filter_combinations) {
  cat("\n=== Processing:", combo$label, "===\n")
  
  # start with full dataset
  filtered_data <- mut2
  
  # apply filters based on combination
  if ("gnomad" %in% combo$filters) {
    filtered_data <- apply_gnomad_filter(filtered_data)
    cat("After gnomAD filter:", nrow(filtered_data), "mutations\n")
  }
  
  if ("vaf" %in% combo$filters) {
    filtered_data <- apply_vaf_filter(filtered_data)
    cat("After VAF filter:", nrow(filtered_data), "mutations\n")
  }
  
  if ("m6a" %in% combo$filters) {
    filtered_data <- apply_m6a_filter(filtered_data)
    cat("After m6A filter:", nrow(filtered_data), "mutations\n")
  }
  
  if ("merip" %in% combo$filters) {
    filtered_data <- apply_merip_filter(filtered_data)
    cat("After MeRIP filter:", nrow(filtered_data), "mutations\n")
  }
  
  if ("m6a_or_merip" %in% combo$filters) {
    filtered_data <- apply_m6a_or_merip_filter(filtered_data)
    cat("After m6A OR MeRIP filter:", nrow(filtered_data), "mutations\n")
  }
  
  # Deduplicate to unique patient-gene combinations
  filtered_data <- filtered_data %>%
    filter(!duplicated(mut_id2))
  
  cat("After deduplication:", nrow(filtered_data), "unique patient-gene combos\n")
  
  # gene summary
  gene_summary <- filtered_data %>%
    group_by(Hugo_Symbol) %>%
    summarise(
      n_mutations = n(),
      n_distinct_mutations = n_distinct(paste(Chromosome, Start_Position, End_Position, 
                                              Reference_Allele, Tumor_Seq_Allele2, sep = "_")),
      n_patients = n_distinct(Kids.First.Participant.ID),
      n_tumor_types = n_distinct(harmonized_diagnosis),
      median_fracT = median(fracT, na.rm = TRUE),
      median_fracN = median(fracN, na.rm = TRUE),
      mean_gnomad_AF = mean(gnomad_3_1_1_AF, na.rm = TRUE),
      tumor_types = paste(unique(harmonized_diagnosis), collapse = "; "),
      .groups = "drop"
    ) %>%
    arrange(desc(n_patients), desc(n_mutations))
  
  # Save summary table
  write_csv(gene_summary, 
            paste0("filter_combinations_output_3/gene_summary_", combo$name, ".csv"))
  
  # visualization for top 20 genes
  if (nrow(gene_summary) > 0) {
    top_n_genes <- min(20, nrow(gene_summary))
    
    plot <- gene_summary %>%
      head(top_n_genes) %>%
      ggplot(aes(x = reorder(Hugo_Symbol, n_patients), y = n_patients)) +
      geom_col(aes(fill = n_tumor_types), width = 0.7) +
      geom_text(aes(label = paste0(n_patients, " pts\n", 
                                   n_distinct_mutations, " muts\n",
                                   n_tumor_types, " types")),
                hjust = -0.1, size = 2.5, lineheight = 0.9) +
      coord_flip() +
      scale_fill_gradient(low = "lightblue", high = "darkblue", 
                          name = "# Tumor\nTypes") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
      labs(title = paste0("Top ", top_n_genes, " Genes with A→X Mutations"),
           subtitle = paste0("Filter: ", combo$label, 
                             "\nTotal unique patient-gene combinations: ", nrow(filtered_data)),
           x = "Gene", 
           y = "Number of Patients") +
      theme_minimal(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, color = "gray30"),
        axis.text.y = element_text(size = 9),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()
      )
    
    # Save plot
    ggsave(
      filename = paste0("filter_combinations_output_3/top_genes_", combo$name, ".png"),
      plot = plot,
      width = 10,
      height = 8,
      dpi = 300
    )
    
    # save as PDF
    ggsave(
      filename = paste0("filter_combinations_output_3/top_genes_", combo$name, ".pdf"),
      plot = plot,
      width = 10,
      height = 8
    )
    
    cat("Plot saved for", combo$name, "\n")
  } else {
    cat("No data to plot for", combo$name, "\n")
  }
  
  # Store results for comparison
  all_results[[combo$name]] <- list(
    n_mutations = nrow(filtered_data),
    n_genes = n_distinct(filtered_data$Hugo_Symbol),
    n_patients = n_distinct(filtered_data$Kids.First.Participant.ID),
    n_tumor_types = n_distinct(filtered_data$harmonized_diagnosis),
    gene_summary = gene_summary
  )
}

# Create comparison summary table
comparison_df <- data.frame(
  Filter_Combination = sapply(filter_combinations, function(x) x$name),
  Filter_Label = sapply(filter_combinations, function(x) x$label),
  N_Mutations = sapply(all_results, function(x) x$n_mutations),
  N_Genes = sapply(all_results, function(x) x$n_genes),
  N_Patients = sapply(all_results, function(x) x$n_patients),
  N_Tumor_Types = sapply(all_results, function(x) x$n_tumor_types)
)
comparison_df
write_csv(comparison_df, "filter_combinations_output_3/comparison_summary.csv")
print(comparison_df)

# comparison barplot
comparison_plot <- comparison_df %>%
  pivot_longer(cols = c(N_Mutations, N_Genes, N_Patients), 
               names_to = "Metric", values_to = "Count") %>%
  mutate(Metric = factor(Metric, levels = c("N_Mutations", "N_Genes", "N_Patients"),
                         labels = c("Mutations", "Genes", "Patients"))) %>%
  ggplot(aes(x = reorder(Filter_Label, -Count), y = Count, fill = Metric)) +
  geom_col(position = "dodge") +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Impact of Filter Combinations on Dataset Size",
       subtitle = "Comparing 16 filter combinations across mutations, genes, and patients",
       x = "Filter Combination",
       y = "Count",
       fill = "") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 8),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, color = "gray30"))

ggsave("filter_combinations_output_3/filter_comparison_plot.png", 
       comparison_plot, width = 12, height = 10, dpi = 300)

ggsave("filter_combinations_output_3/filter_comparison_plot.pdf", 
       comparison_plot, width = 12, height = 10)
comparison_plot
