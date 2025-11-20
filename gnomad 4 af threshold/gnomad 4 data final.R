# FUNCTION: Merge updated missing variants back into original dataframe
merge_updated_variants <- function(original_df, updated_missing_df) {
  
  library(dplyr)
  
  cat("=== MERGING UPDATED VARIANTS ===\n\n")
  
  # Initial counts
  cat(sprintf("Original dataframe: %d variants\n", nrow(original_df)))
  cat(sprintf("Updated missing variants: %d variants\n", nrow(updated_missing_df)))
  
  original_missing <- sum(original_df$gnomAD_4_AF == 0 | is.na(original_df$gnomAD_4_AF), 
                          na.rm = TRUE)
  cat(sprintf("Original missing (AF=0 or NA): %d variants\n\n", original_missing))
  
  # Identify rows to replace (where gnomAD_4_AF was 0 or NA)
  rows_to_replace <- original_df %>%
    filter(gnomAD_4_AF == 0 | is.na(gnomAD_4_AF))
  
  cat(sprintf("Rows to replace: %d\n", nrow(rows_to_replace)))
  
  # Keep rows that were already annotated (AF > 0)
  rows_to_keep <- original_df %>%
    filter(gnomAD_4_AF > 0 & !is.na(gnomAD_4_AF))
  
  cat(sprintf("Rows to keep (already annotated): %d\n\n", nrow(rows_to_keep)))
  
  # Merge: Keep good rows + add all updated rows
  # This replaces old missing rows with updated data
  merged_df <- bind_rows(rows_to_keep, updated_missing_df) %>%
    arrange(CHROM, POS)
  
  # Summary statistics
  cat("=== MERGE COMPLETE ===\n\n")
  cat(sprintf("Final dataframe: %d variants\n", nrow(merged_df)))
  cat(sprintf("Now annotated (AF>0): %d variants\n", 
              sum(merged_df$gnomAD_4_AF > 0, na.rm = TRUE)))
  cat(sprintf("Rare but in gnomAD (AF=0, AN>0): %d variants\n",
              sum(merged_df$gnomAD_4_AF == 0 & merged_df$gnomAD_4_AN > 0, na.rm = TRUE)))
  cat(sprintf("Still missing (AF=0/NA, AN=0/NA): %d variants\n",
              sum((merged_df$gnomAD_4_AF == 0 | is.na(merged_df$gnomAD_4_AF)) & 
                    (merged_df$gnomAD_4_AN == 0 | is.na(merged_df$gnomAD_4_AN)), 
                  na.rm = TRUE)))
  
  cat(sprintf("\nImprovement: %d variants newly annotated\n", 
              original_missing - sum((merged_df$gnomAD_4_AF == 0 | is.na(merged_df$gnomAD_4_AF)) & 
                                       (merged_df$gnomAD_4_AN == 0 | is.na(merged_df$gnomAD_4_AN)), 
                                     na.rm = TRUE)))
  
  return(merged_df)
}


# ============================================================
# USAGE EXAMPLE
# ============================================================

# Load your dataframes
annotated_v2_final_p1 <- read.csv("annotated_v2_final_p1.csv", 
                                  header = TRUE, 
                                  stringsAsFactors = FALSE)

missing_only <- read.csv("annotated_v2_final.csv", 
                         header = TRUE, 
                         stringsAsFactors = FALSE)

# Merge the updated variants back in
final_merged_1107 <- merge_updated_variants(
  original_df = annotated_v2_final_p1,
  updated_missing_df = missing_only
)

# Save the final merged result
write.csv(final_merged_1107, "annotated_1107.csv", row.names = FALSE)
cat("\nSaved to: annotated_1107.csv\n")

#annotation into deduplicated file
final_merged_1_1107 <- file_renamed2 %>%
  left_join(
    final_merged_1107,
    by = c("Chromosome" = "CHROM", 
           "Start_Position" = "POS",
           "Reference_Allele" = "REF",
           "Tumor_Seq_Allele2" = "ALT"),
    relationship = "many-to-many"
  )
write.csv(final_merged_1_1107, "annotated_1_1107.csv", row.names = FALSE)

#annotation into final file
final_file_with_all_rows_1107 <- file_renamed %>%
  left_join(
    final_merged_1_1107,
    by = c("Chromosome" = "CHROM", 
           "Start_Position" = "POS",
           "Reference_Allele" = "REF",
           "Tumor_Seq_Allele2" = "ALT"),
    relationship = "many-to-many"
  )

write.csv(final_file_with_all_rows_1107, "annotated_1_complete_1107.csv", row.names = FALSE)

colnames(final_file_with_all_rows_1107)

# Count how many rows have NA in gnomAD_4_AF
na_gnomad4af_count <- final_file_with_all_rows_1107 %>%
  filter(is.na(gnomAD_4_AF)) %>%
  nrow()

# Count how many rows have gnomAD_4_in_gnomAD == FALSE, NA, or blank
bad_gnomad4in_count <- final_file_with_all_rows_1107 %>%
  filter(is.na(gnomAD_4_in_gnomAD) |
           gnomAD_4_in_gnomAD == FALSE |
           gnomAD_4_in_gnomAD == "" ) %>%
  nrow()

# Print results
cat("Rows with gnomAD_4_AF == NA:", na_gnomad4af_count, "\n")
cat("Rows with gnomAD_4_in_gnomAD == FALSE / NA / blank:", bad_gnomad4in_count, "\n")

#---------
# Create a new dataframe with only rows where gnomAD_4_AF is NA
# and include only the renamed columns
gnomad4af_na_subset <- final_file_with_all_rows_1107 %>%
  filter(is.na(gnomAD_4_AF)) %>%
  mutate(
    CHROM = Chromosome,
    POS = Start_Position,
    REF = Reference_Allele,
    ALT = Tumor_Seq_Allele2
  ) %>%
  dplyr::select(CHROM, POS, REF, ALT)

# Run retry again with conservative settings
annotated_na <- annotate_variants_batch(gnomad4af_na_subset)
#----------

library(dplyr)

final_file_with_all_rows_1107 <- final_file_with_all_rows_1107 %>%
  mutate(gnomAD_4_AF = ifelse(is.na(gnomAD_4_AF), 0, gnomAD_4_AF))
write.csv(final_file_with_all_rows_1107, "annotated_1_complete_1107.csv", row.names = FALSE)
