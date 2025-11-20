library(dplyr)
library(tidyr)
library(stringr)
#erasers_subset <- read.csv("C:/Users/isr4005/Downloads/erasers_subset.csv", stringsAsFactors = FALSE)

# Step 1 Identify site-specific vs gene-only rows
# Site-specific = has chr, start, end coordinates
erasers_clean <- erasers_subset %>%
  mutate(
    has_site = !is.na(chr) & !is.na(start) & !is.na(end) & 
      chr != "" & !is.na(as.numeric(start)) & !is.na(as.numeric(end)),
    #  gene-site key for site-specific rows
    gene_site_key = ifelse(has_site, 
                           paste(target_gene, chr, start, end, strand, sep = "_"), 
                           NA_character_)
  )

# Step 2 Split into site-specific and gene-only
site_specific_data <- erasers_clean %>% filter(has_site == TRUE)
gene_only_data <- erasers_clean %>% filter(has_site == FALSE)

# Step 3 For each WER, process the data
process_wer <- function(wer_name, site_df, gene_df) {
  
  # Filter for this WER
  site_wer <- site_df %>% filter(wer == wer_name)
  gene_wer <- gene_df %>% filter(wer == wer_name)
  
  results <- list()
  
  #  all unique genes for this WER (from both site and gene-only)
  all_genes <- unique(c(site_wer$target_gene, gene_wer$target_gene))
  
  for (g in all_genes) {
    # Site-specific rows for this gene
    site_rows <- site_wer %>% filter(target_gene == g)
    # Gene-only rows for this gene
    gene_rows <- gene_wer %>% filter(target_gene == g)
    
    # Process each unique site for this gene
    if (nrow(site_rows) > 0) {
      unique_sites <- unique(site_rows$gene_site_key)
      
      for (site_key in unique_sites) {
        site_data <- site_rows %>% filter(gene_site_key == site_key)
        
        #  site-specific + gene-only evidence for this gene-site row
        combined <- bind_rows(site_data, gene_rows)
        
        # Aggregate into one row
        row_result <- aggregate_row(combined, wer_name, g, site_key, 
                                    source_type = ifelse(nrow(gene_rows) > 0, "both", "gene_site_only"))
        results <- c(results, list(row_result))
      }
    }
    
    #  gene-only row (separate from site-specific)
    if (nrow(gene_rows) > 0) {
      row_result <- aggregate_row(gene_rows, wer_name, g, site_key = NA_character_, 
                                  source_type = "gene_only")
      results <- c(results, list(row_result))
    }
  }
  
  if (length(results) > 0) {
    bind_rows(results)
  } else {
    NULL
  }
}

# Step 4 Aggregation function - combines all evidence for a row
aggregate_row <- function(df, wer_name, gene, site_key, source_type) {
  
  # Combine all evidence types
  all_evidence <- paste(unique(df$evidence[!is.na(df$evidence)]), collapse = "; ")
  
  # Boolean flags based on evidence
  has_validation <- any(str_detect(tolower(df$evidence), "validat"), na.rm = TRUE)
  has_perturbation <- any(str_detect(tolower(df$evidence), "perturb"), na.rm = TRUE)
  has_binding <- any(str_detect(tolower(df$evidence), "binding"), na.rm = TRUE)
  
  # Get site info (if available)
  site_info <- df %>% filter(!is.na(chr) & chr != "") %>% slice(1)
  
  # Aggregate log2fc values take first non-NA for each type
  get_first <- function(x) { x <- x[!is.na(x) & x != ""]; if(length(x) > 0) x[1] else NA }
  
  tibble(
    wer = wer_name,
    wer_type = "Eraser",
    target_gene = gene,
    gene_site_key = site_key,
    chr = get_first(df$chr),
    strand = get_first(df$strand),
    start = get_first(df$start),
    end = get_first(df$end),
    source = source_type,
    evidence_combined = all_evidence,
    has_validation = has_validation,
    has_perturbation = has_perturbation,
    has_binding = has_binding,
    # Log2FC columns
    log2fc_te = get_first(df$log2fc_te),
    log2fc_exp = get_first(df$log2fc_exp),
    log2fc_stab = get_first(df$log2fc_stab),
    log2fc_modrip = get_first(df$log2fc_modrip),
    log2fc_as = get_first(df$log2fc_as),
    # Perturbation columns
    perturbation_te = get_first(df$perturbation_te),
    perturbation_exp = get_first(df$perturbation_exp),
    perturbation_stab = get_first(df$perturbation_stab),
    perturbation_modrip = get_first(df$perturbation_modrip),
    perturbation_as = get_first(df$perturbation_as),
    # Direction columns
    direction_te = get_first(df$direction_te),
    direction_exp = get_first(df$direction_exp),
    direction_stab = get_first(df$direction_stab),
    direction_modrip = get_first(df$direction_modrip),
    # Other
    inc_level_difference = get_first(df$inc_level_difference),
    uniprot_name = get_first(df$uniprot_name)
  )
}

# Step 5 Process all WERs
all_wers <- unique(erasers_clean$wer)
final_results <- list()

for (w in all_wers) {
  cat("Processing WER:", w, "\n")
  result <- process_wer(w, site_specific_data, gene_only_data)
  if (!is.null(result)) {
    final_results <- c(final_results, list(result))
  }
}

# Combine all results
erasers_wide <- bind_rows(final_results)
head(erasers_wide)
str(erasers_wide)
cat("\n--- Summary ---\n")
cat("Total rows:", nrow(erasers_wide), "\n")
cat("Unique WERs:", n_distinct(erasers_wide$wer), "\n")
cat("Unique genes:", n_distinct(erasers_wide$target_gene), "\n")
table(erasers_wide$source)
table(erasers_wide$has_validation)
table(erasers_wide$has_perturbation)
table(erasers_wide$has_binding)

# Save
write.csv(erasers_wide, "C:/Users/isr4005/Downloads/erasers_wide_format.csv", row.names = FALSE)
saveRDS(erasers_wide, "C:/Users/isr4005/Downloads/erasers_wide_format.rds")