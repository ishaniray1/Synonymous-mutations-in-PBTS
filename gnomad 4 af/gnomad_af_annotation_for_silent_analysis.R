colnames(silent_variants_df)
# Load file
file <- read_csv("silent_variants_df.csv")

# Create copy with renamed columns for annotation
file_renamed <- file %>%
  mutate(
    CHROM = Chromosome,
    POS = Start_Position,
    REF = Reference_Allele,
    ALT = Tumor_Seq_Allele2
  )
# Verify required columns 
required_cols <- c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2")
if (!all(required_cols %in% colnames(file))) {
  stop("Missing required columns!")
}
# Run annotation 
annotated <- annotate_variants_batch(file_renamed)
# Merge back to original 
final_file <- file %>%
  left_join(
    annotated,
    by = c("Chromosome" = "CHROM", 
           "Start_Position" = "POS",
           "Reference_Allele" = "REF",
           "Tumor_Seq_Allele2" = "ALT")
  )

# Save
write.table(final_file, "silent_variants_df_with_gnomad_v4.csv", 
            row.names = FALSE)
# OR  as tab-delimited
write.table(final_file, "silent_variants_with_gnomad_v4.txt", 
             sep = "\t", quote = FALSE, row.names = FALSE)