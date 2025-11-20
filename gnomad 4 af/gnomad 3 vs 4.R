# Merge the dataframes based on matching chromosome and position columns
merged_data <- merge(
  silents_variants_df,
  missing_only,
  by.x = c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2"),
  by.y = c("CHROM", "POS", "REF", "ALT"),
  all.x = FALSE,  # Only keep matching rows
  suffixes = c("_silent", "_missing")
)

# Check how many rows matched
cat("Total matching rows:", nrow(merged_data), "\n\n")

# Analyze gnomAD_4_in_gnomAD status
cat("gnomAD_4_in_gnomAD status:\n")
print(table(merged_data$gnomAD_4_in_gnomAD, useNA = "ifany"))
cat("\n")

# Analyze dbSNP_RS status (novel vs has RS ID)
# Novel means the cell contains "novel"
merged_data$dbSNP_novel <- tolower(trimws(merged_data$dbSNP_RS)) == "novel"

cat("dbSNP_RS status:\n")
cat("Novel:", sum(merged_data$dbSNP_novel), "\n")
cat("Has RS ID:", sum(!merged_data$dbSNP_novel), "\n\n")

# Cross-tabulation: gnomAD FALSE vs dbSNP novel status
cat("Cross-tabulation of gnomAD_4_in_gnomAD and dbSNP novel status:\n")
cross_tab <- table(
  gnomAD_in_db = merged_data$gnomAD_4_in_gnomAD,
  dbSNP_novel = merged_data$dbSNP_novel
)
print(cross_tab)
cat("\n")

# Specifically look at variants where gnomAD_4_in_gnomAD is FALSE
gnomad_false <- merged_data[merged_data$gnomAD_4_in_gnomAD == FALSE, ]
cat("Among variants with gnomAD_4_in_gnomAD = FALSE:\n")
cat("Total:", nrow(gnomad_false), "\n")
cat("Novel in dbSNP:", sum(gnomad_false$dbSNP_novel), "\n")
cat("Not novel in dbSNP:", sum(!gnomad_false$dbSNP_novel), "\n\n")

# Summary statistics
cat("Summary:\n")
cat("Variants in gnomAD (TRUE) and novel in dbSNP:", 
    sum(merged_data$gnomAD_4_in_gnomAD == TRUE & merged_data$dbSNP_novel), "\n")
cat("Variants in gnomAD (TRUE) and not novel in dbSNP:", 
    sum(merged_data$gnomAD_4_in_gnomAD == TRUE & !merged_data$dbSNP_novel), "\n")
cat("Variants NOT in gnomAD (FALSE) and novel in dbSNP:", 
    sum(merged_data$gnomAD_4_in_gnomAD == FALSE & merged_data$dbSNP_novel), "\n")
cat("Variants NOT in gnomAD (FALSE) and not novel in dbSNP:", 
    sum(merged_data$gnomAD_4_in_gnomAD == FALSE & !merged_data$dbSNP_novel), "\n")


# Filter for variants where gnomAD_4_in_gnomAD is FALSE and dbSNP_RS is not "novel"
not_in_gnomad_has_rs <- merged_data[
  merged_data$gnomAD_4_in_gnomAD == FALSE & !merged_data$dbSNP_novel,
]

# Get the list of RS IDs
rs_ids <- not_in_gnomad_has_rs$dbSNP_RS

# Display the results
cat("Number of variants NOT in gnomAD but HAS RS ID:", length(rs_ids), "\n\n")

cat("RS IDs:\n")
print(rs_ids)

# If you want unique RS IDs (in case there are duplicates)
cat("\n\nUnique RS IDs:\n")
unique_rs_ids <- unique(rs_ids)
print(unique_rs_ids)
cat("\nTotal unique RS IDs:", length(unique_rs_ids), "\n")

# Optional: Save to a file
write.table(
  data.frame(RS_ID = rs_ids),
  "gnomad_false_with_rs_ids.txt",
  row.names = FALSE,
  quote = FALSE
)

# Optional: View the full data for these variants
View(gnomad_false_with_rs)
