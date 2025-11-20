# Check for duplicates in annotated_eg
annotated_eg %>%
  group_by(CHROM, POS, REF, ALT) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1) %>%
  arrange(desc(n))

# Check for duplicates in file
file %>%
  group_by(Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1) %>%
  arrange(desc(n))

# Check what's actually happening in the join
file %>%
  group_by(Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2) %>%
  summarise(file_count = n(), .groups = "drop") %>%
  left_join(
    annotated_eg %>%
      group_by(CHROM, POS, REF, ALT) %>%
      summarise(annotated_count = n(), .groups = "drop"),
    by = c("Chromosome" = "CHROM", 
           "Start_Position" = "POS",
           "Reference_Allele" = "REF",
           "Tumor_Seq_Allele2" = "ALT")
  ) %>%
  filter(is.na(annotated_count))  # These are missing from annotated_eg

# Remove duplicates from annotated data (keep first occurrence)
annotated_eg_dedup <- annotated_eg %>%
  distinct(CHROM, POS, REF, ALT, .keep_all = TRUE)

# Now join - should preserve all 64k rows
final_file_again <- file %>%
  left_join(
    annotated_eg_dedup,
    by = c("Chromosome" = "CHROM", 
           "Start_Position" = "POS",
           "Reference_Allele" = "REF",
           "Tumor_Seq_Allele2" = "ALT"),
    suffix = c("", ".from_df1")
  )

# Verify
nrow(final_file_again)  # Should be 64,936