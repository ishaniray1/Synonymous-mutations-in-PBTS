# ========================================
# SUPPORTED INPUT FILE FORMATS
# ========================================

# FORMAT 1: Simple CSV (MINIMUM INFO NEEDED)
# -----------------------------------------
# File: variants.csv
#
# CHROM,POS,REF,ALT
# 1,55516888,G,A
# 7,140453136,A,T
# 17,7577538,C,T
# X,123456789,G,C

# FORMAT 2: CSV with Gene Information (RECOMMENDED)
# -----------------------------------------
# File: variants_with_genes.csv
#
# Gene,CHROM,POS,REF,ALT,Type
# EGFR,7,55249071,G,A,Missense
# TP53,17,7577538,C,T,Nonsense
# BRCA1,17,43094692,G,A,Silent

# FORMAT 3: VCF Format (Standard)
# -----------------------------------------
# File: variants.vcf
#
# #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
# 1       55516888 .      G       A       .       .       .
# 7       140453136 .     A       T       .       .       .
# 17      7577538  rs28934576 C   T       .       .       .

# FORMAT 4: Tab-delimited (TSV)
# -----------------------------------------
# File: variants.tsv
#
# CHROM  POS         REF  ALT  GENE    MUTATION_TYPE
# 1      55516888    G    A    PCSK9   Missense
# 7      55249071    G    A    EGFR    Missense


# ========================================
# FLEXIBLE INPUT READER
# ========================================

read_variant_file <- function(file_path) {
  
  # Detect file format
  if (grepl("\\.vcf$|\\.vcf\\.gz$", file_path)) {
    cat("Detected VCF format\n")
    variants <- read.table(file_path, comment.char = "#", 
                           stringsAsFactors = FALSE)
    
    # Standard VCF columns
    if (ncol(variants) >= 5) {
      colnames(variants)[1:5] <- c("CHROM", "POS", "ID", "REF", "ALT")
    }
    
  } else if (grepl("\\.csv$", file_path)) {
    cat("Detected CSV format\n")
    variants <- read.csv(file_path, stringsAsFactors = FALSE)
    
  } else if (grepl("\\.tsv$|\\.txt$", file_path)) {
    cat("Detected TSV format\n")
    variants <- read.table(file_path, header = TRUE, 
                           sep = "\t", stringsAsFactors = FALSE)
  } else {
    cat("Unknown format, attempting tab-delimited...\n")
    variants <- read.table(file_path, header = TRUE, 
                           stringsAsFactors = FALSE)
  }
  
  # Validate required columns
  required_cols <- c("CHROM", "POS", "REF", "ALT")
  missing_cols <- setdiff(required_cols, colnames(variants))
  
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", 
               paste(missing_cols, collapse = ", "),
               "\nRequired: CHROM, POS, REF, ALT"))
  }
  
  # Clean chromosome names (remove 'chr' prefix if present)
  variants$CHROM <- gsub("chr", "", variants$CHROM)
  
  # Ensure POS is numeric
  variants$POS <- as.numeric(variants$POS)
  
  cat(sprintf("Successfully loaded %d variants\n", nrow(variants)))
  cat(sprintf("Columns detected: %s\n", paste(colnames(variants), collapse = ", ")))
  
  return(variants)
}


# ========================================
# EXAMPLE: CREATE TEMPLATE FILES
# ========================================

create_template_csv <- function(output_file = "variant_template.csv") {
  template <- data.frame(
    Gene = c("EGFR", "TP53", "KRAS", "BRCA1"),
    CHROM = c("7", "17", "12", "17"),
    POS = c(55249071, 7577538, 25398284, 43094692),
    REF = c("G", "C", "C", "G"),
    ALT = c("A", "T", "T", "A"),
    Mutation_Type = c("Missense", "Nonsense", "Missense", "Silent"),
    Sample_ID = c("Sample1", "Sample2", "Sample1", "Sample3")
  )
  
  write.csv(template, output_file, row.names = FALSE)
  cat(sprintf("Template file created: %s\n", output_file))
  cat("Fill in your variant information following this format\n")
}


# ========================================
# CHECK YOUR FILE BEFORE RUNNING
# ========================================

validate_variant_file <- function(file_path) {
  
  cat("=== Validating Variant File ===\n\n")
  
  # Read file
  tryCatch({
    variants <- read_variant_file(file_path)
  }, error = function(e) {
    cat("ERROR reading file:", e$message, "\n")
    return(FALSE)
  })
  
  # Check required columns
  cat("✓ Required columns present\n")
  
  # Check data types
  if (!is.numeric(variants$POS)) {
    cat("WARNING: POS column is not numeric\n")
  } else {
    cat("✓ POS column is numeric\n")
  }
  
  # Check for missing values
  missing <- sapply(variants[, c("CHROM", "POS", "REF", "ALT")], 
                    function(x) sum(is.na(x)))
  if (sum(missing) > 0) {
    cat("WARNING: Missing values detected:\n")
    print(missing[missing > 0])
  } else {
    cat("✓ No missing values in required columns\n")
  }
  
  # Check chromosome format
  unique_chroms <- unique(variants$CHROM)
  cat(sprintf("✓ Chromosomes detected: %s\n", 
              paste(head(unique_chroms, 10), collapse = ", ")))
  
  # Check for duplicates
  dup_count <- sum(duplicated(paste(variants$CHROM, variants$POS, 
                                    variants$REF, variants$ALT)))
  if (dup_count > 0) {
    cat(sprintf("WARNING: %d duplicate variants detected\n", dup_count))
  } else {
    cat("✓ No duplicate variants\n")
  }
  
  # Sample preview
  cat("\n=== First 5 variants preview ===\n")
  print(head(variants[, c("CHROM", "POS", "REF", "ALT")], 5))
  
  cat("\n✓ File is ready for gnomAD annotation!\n")
  
  return(TRUE)
}


# ========================================
# CONVERT BETWEEN FORMATS
# ========================================

# Convert MAF (Mutation Annotation Format) to simple format
convert_maf_to_simple <- function(maf_file, output_file = "variants_simple.csv") {
  
  maf <- read.table(maf_file, header = TRUE, sep = "\t", 
                    comment.char = "#", stringsAsFactors = FALSE)
  
  # MAF standard columns
  simple <- data.frame(
    Gene = maf$Hugo_Symbol,
    CHROM = gsub("chr", "", maf$Chromosome),
    POS = maf$Start_Position,
    REF = maf$Reference_Allele,
    ALT = maf$Tumor_Seq_Allele2,  # or Tumor_Seq_Allele1
    Variant_Classification = maf$Variant_Classification,
    Sample = maf$Tumor_Sample_Barcode,
    stringsAsFactors = FALSE
  )
  
  write.csv(simple, output_file, row.names = FALSE)
  cat(sprintf("Converted %d variants from MAF to simple format\n", nrow(simple)))
  cat(sprintf("Output: %s\n", output_file))
  
  return(simple)
}


# Convert Annovar output to simple format
convert_annovar_to_simple <- function(annovar_file, output_file = "variants_simple.csv") {
  
  annovar <- read.table(annovar_file, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)
  
  simple <- data.frame(
    Gene = annovar$Gene.refGene,
    CHROM = gsub("chr", "", annovar$Chr),
    POS = annovar$Start,
    REF = annovar$Ref,
    ALT = annovar$Alt,
    Variant_Type = annovar$ExonicFunc.refGene,
    stringsAsFactors = FALSE
  )
  
  write.csv(simple, output_file, row.names = FALSE)
  cat(sprintf("Converted %d variants from Annovar to simple format\n", nrow(simple)))
  
  return(simple)
}


# ========================================
# QUICK START EXAMPLES
# ========================================

# Example 1: If you have a simple list
your_variants <- data.frame(
  CHROM = c("7", "17", "12"),
  POS = c(55249071, 7577538, 25398284),
  REF = c("G", "C", "C"),
  ALT = c("A", "T", "T")
)

# Save and annotate
# write.csv(your_variants, "my_variants.csv", row.names = FALSE)
# results <- run_complete_analysis("my_variants.csv")


# Example 2: If you have gene names and want to include them
your_variants_with_genes <- data.frame(
  Gene = c("EGFR", "TP53", "KRAS"),
  CHROM = c("7", "17", "12"),
  POS = c(55249071, 7577538, 25398284),
  REF = c("G", "C", "C"),
  ALT = c("A", "T", "T"),
  Mutation_Type = c("Missense", "Nonsense", "Missense")
)


# Example 3: Validate before running
# validate_variant_file("my_6000_variants.csv")
# results <- run_complete_analysis("my_6000_variants.csv")


# ========================================
# WHAT INFO COMES BACK FROM gnomAD
# ========================================

# After annotation, you'll get these columns ADDED:
# - gnomAD_genome_AF: Allele frequency in genome data
# - gnomAD_exome_AF: Allele frequency in exome data  
# - gnomAD_max_AF: Maximum of the two above
# - gnomAD_genome_AC: Allele count (how many people have it)
# - gnomAD_genome_AN: Allele number (total chromosomes tested)
# - in_gnomAD: TRUE/FALSE if variant found in gnomAD
# - variant_class: Classification (Common/Rare/Very_Rare/Not_in_gnomAD)

# Plus you keep all your original columns (Gene, Mutation_Type, etc.)