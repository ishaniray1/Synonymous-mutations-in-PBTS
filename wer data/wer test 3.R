# Filter for writers (case-insensitive) 
witers_subset <- subset(all_rm2target, tolower(wer_type) == "writer")
write.csv(witers_subset, "C:/Users/isr4005/Downloads/witers_subset.csv", row.names = FALSE)
# Filter for readers (case-insensitive)
readers_subset <- subset(all_rm2target, tolower(wer_type) == "reader")
write.csv(readers_subset, "C:/Users/isr4005/Downloads/readers_subset.csv", row.names = FALSE)
# Filter for erasers (case-insensitive)
erasers_subset <- subset(all_rm2target, tolower(wer_type) == "eraser")
write.csv(erasers_subset, "C:/Users/isr4005/Downloads/erasers_subset.csv", row.names = FALSE)

# identify gene-only entries (no site info)
gene_only_entries_eraser <- erasers_subset %>%
  filter(is.na(chr) | is.na(start) | is.na(end)) %>%
  distinct(wer, wer_type, target_gene, .keep_all = TRUE) %>%
  select(
    wer, wer_type, target_gene, perturbation, direction,
     log2fc_exp, log2fc_stab, log2fc_te, log2fc_modrip,
    inc_level_difference, uniprot_name, evidence
  )

#  distinct gene-site combinations with their writers
gene_site_writers_eraser <- erasers_subset %>%
  filter(!is.na(chr) & !is.na(start) & !is.na(end)) %>%
  distinct(target_gene, chr, strand, start, end, wer, wer_type, .keep_all = TRUE)

# For each gene-site combination, add all the gene-only writers
expanded_gene_site_eraser <- gene_site_writers_eraser %>%
  # Get all distinct gene-site combinations
  distinct(target_gene, chr, strand, start, end) %>%
  # Cross join with gene-only writers for the same gene
  left_join(gene_only_entries_eraser, by = "target_gene", relationship = "many-to-many") %>%
  # Combine with the original site-specific writers
  bind_rows(
    gene_site_writers_eraser %>%
      select(
        target_gene, chr, strand, start, end, wer, wer_type, perturbation,
        direction,
        log2fc_exp, log2fc_stab, log2fc_te, log2fc_modrip,
        inc_level_difference, uniprot_name, evidence
      )
  ) %>%
  # DON'T remove duplicates - add evidence type columns
  mutate(
    validation = ifelse(grepl("validation", evidence, ignore.case = TRUE), TRUE, FALSE),
    binding = ifelse(grepl("binding", evidence, ignore.case = TRUE), TRUE, FALSE),
    perturbation_evidence = ifelse(grepl("perturbation", evidence, ignore.case = TRUE), TRUE, FALSE)
  ) %>%
  # Add source column to track where each entry came from
  group_by(target_gene, chr, strand, start, end, wer, wer_type) %>%
  mutate(
    source = case_when(
      n() > 1 & !is.na(first(na.omit(evidence))) ~ "both",
      !is.na(first(na.omit(evidence))) ~ "site_specific",
      TRUE ~ "gene_general"
    ),
    # Combine evidence from multiple entries
    evidence_combined = ifelse(n() > 1,
                               paste(na.omit(unique(evidence)), collapse = ";"),
                               first(evidence))
  ) %>%
  ungroup()

#  include pure gene-only rows (for genes that only have general)
pure_gene_only_eraser <- gene_only_entries_eraser %>%
  anti_join(gene_site_writers_eraser, by = "target_gene") %>%
  mutate(
    chr = NA_character_,
    strand = NA_character_,
    start = NA_integer_,
    end = NA_integer_,
    validation = ifelse(grepl("validation", evidence, ignore.case = TRUE), TRUE, FALSE),
    binding = ifelse(grepl("binding", evidence, ignore.case = TRUE), TRUE, FALSE),
    perturbation_evidence = ifelse(grepl("perturbation", evidence, ignore.case = TRUE), TRUE, FALSE),
    source = "gene_only",
    evidence_combined = evidence
  )

# Combine everything
final_long_eraser <- bind_rows(expanded_gene_site_eraser, pure_gene_only_eraser) %>%
  arrange(target_gene, chr, start, end, wer, wer_type) %>%
  # For the final output with numbered columns
  group_by(target_gene, chr, strand, start, end) %>%
  mutate(eraser_number = row_number()) %>%
  ungroup() %>%
  
  # Create the numbered columns
  select(
    target_gene, chr, strand, start, end, eraser_number,
    wer, wer_type, perturbation, direction,
    log2fc_exp, log2fc_stab, log2fc_te, log2fc_modrip,
    inc_level_difference, uniprot_name,
    evidence_combined, validation, binding,
    perturbation_evidence, source
  ) %>%
  
  # Pivot to wider format
  pivot_wider(
    id_cols = c(target_gene, chr, strand, start, end),
    names_from = eraser_number,
    values_from = c(  wer, wer_type, perturbation, direction,
                      log2fc_exp, log2fc_stab, log2fc_te, log2fc_modrip,
                      inc_level_difference, uniprot_name,
                      evidence_combined, validation, binding,
                      perturbation_evidence, source
    ),
    names_glue = "{.value}_{eraser_number}"
  )