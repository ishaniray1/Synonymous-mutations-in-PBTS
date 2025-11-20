if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
packages <- c("VariantAnnotation", "GenomicRanges", "data.table", "dplyr", "vcfR", "biomaRt")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE))
    BiocManager::install(pkg)
}
library(VariantAnnotation)
library(GenomicRanges)
library(data.table)
library(dplyr)
library(vcfR)

# ========================================

annotate_with_gnomad <- function(your_variants_file, gnomad_vcf_file){
  # Read  variants
your_vcf <- read.vcfR(your_variants_file, verbose = FALSE)
your_df <- as.data.frame(your_vcf@fix)

# Read gnomAD data
gnomad_vcf <- read.vcfR(gnomad_vcf_file, verbose = FALSE)
gnomad_df <- as.data.frame(gnomad_vcf@fix)

# Extract gnomAD AF
gnomad_af <- extract.info(gnomad_vcf, element = "AF")
gnomad_df$gnomAD_AF <- as.numeric(gnomad_af)

# Create unique variant IDs
your_df$var_id <- paste(your_df$CHROM, your_df$POS, 
                        your_df$REF, your_df$ALT, sep = "_")
gnomad_df$var_id <- paste(gnomad_df$CHROM, gnomad_df$POS, 
                          gnomad_df$REF, gnomad_df$ALT, sep = "_")

# Merge with gnomAD
annotated <- your_df %>%
  left_join(gnomad_df %>% select(var_id, gnomAD_AF), 
            by = "var_id")

# Replace NA with 0 (variants not in gnomAD are likely v rare)
annotated$gnomAD_AF[is.na(annotated$gnomAD_AF)] <- 0

return(list(
  all_annotated = annotated
))
}

# ========================================

filter_with_bioconductor <- function(your_vcf_file, gnomad_vcf_file) {
  
  # Read VCF files
  your_vcf <- readVcf(your_vcf_file, "hg38")
  gnomad_vcf <- readVcf(gnomad_vcf_file, "hg38")
  
  # Extract gnomAD frequencies
  gnomad_af <- info(gnomad_vcf)$AF
  
  # Create GRanges object
  gnomad_gr <- rowRanges(gnomad_vcf)
  mcols(gnomad_gr)$AF <- gnomad_af
  
  your_gr <- rowRanges(your_vcf)
  
  # Find overlaps
  overlaps <- findOverlaps(your_gr, gnomad_gr)
  
  # Annotate  variants
  your_gr$gnomAD_AF <- NA
  your_gr$gnomAD_AF[queryHits(overlaps)] <- 
    gnomad_gr$AF[subjectHits(overlaps)]
  
  return(your_gr)
}

# ========================================

# ============================================================================
# gnomAD v4 Variant Annotation Script
# ============================================================================
# This script annotates variants with allele frequency data from gnomAD v4
# including population-specific frequencies, homozygote counts, and quality flags
# ============================================================================

# Install required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

packages <- c("VariantAnnotation", "GenomicRanges", "data.table", "dplyr", 
              "vcfR", "biomaRt", "httr", "jsonlite")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE))
    BiocManager::install(pkg)
}

library(VariantAnnotation)
library(GenomicRanges)
library(data.table)
library(dplyr)
library(vcfR)
library(httr)
library(jsonlite)

# Define null coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

# ============================================================================
# FUNCTION: Extract population-specific allele frequencies
# ============================================================================
extract_population_afs <- function(populations_list) {
  pops <- c("afr", "amr", "asj", "eas", "fin", "nfe", "sas", "mid", "ami", "remaining")
  result <- setNames(rep(0, length(pops)), paste0("gnomAD_", toupper(pops), "_AF"))
  
  if (!is.null(populations_list)) {
    for (pop in populations_list) {
      pop_id <- pop$id
      if (pop_id %in% pops) {
        result[paste0("gnomAD_", toupper(pop_id), "_AF")] <- pop$af %||% 0
      }
    }
  }
  
  return(as.list(result))
}

# ============================================================================
# FUNCTION: Annotate variants in batches with checkpoint support
# ============================================================================
annotate_variants_batch <- function(variants_df, 
                                    batch_size = 100,
                                    save_checkpoint = TRUE,
                                    checkpoint_file = "gnomad_checkpoint.rds") {
  
  # Check existing checkpoint
  start_idx <- 1
  if (save_checkpoint && file.exists(checkpoint_file)) {
    cat("Found checkpoint file. Resume? (y/n): ")
    response <- readline()
    if (tolower(response) == "y") {
      checkpoint <- readRDS(checkpoint_file)
      results <- checkpoint$results
      start_idx <- checkpoint$last_idx + 1
      cat(sprintf("Resuming from variant %d\n", start_idx))
    } else {
      results <- list()
    }
  } else {
    results <- list()
  }
  
  total <- nrow(variants_df)
  
  for (i in start_idx:total) {
    chrom <- gsub("chr", "", variants_df$CHROM[i])
    pos <- variants_df$POS[i]
    ref <- variants_df$REF[i]
    alt <- variants_df$ALT[i]
    
    # GraphQL query gnomAD v4
    query <- sprintf('{ 
      variant(variant_id: "%s-%s-%s-%s", dataset: gnomad_r4) {
        variant_id
        
        joint {
          af
          ac
          an
          homozygote_count
          hemizygote_count
          populations {
            id
            af
            ac
            an
            homozygote_count
          }
          faf95 {
            popmax
            popmax_population
          }
        }
        
        genome {
          af
          ac
          an
          homozygote_count
          populations {
            id
            af
            ac
            an
          }
          faf95 {
            popmax
            popmax_population
          }
        }
        
        exome {
          af
          ac
          an
          homozygote_count
          populations {
            id
            af
            ac
            an
          }
        }
        
        flags
        colocated_variants
      }
    }', chrom, pos, ref, alt)
    
    tryCatch({
      response <- POST(
        "https://gnomad.broadinstitute.org/api",
        body = list(query = query),
        encode = "json",
        timeout(10)
      )
      
      data <- content(response, "parsed")
      variant_data <- data$data$variant
      
      if (!is.null(variant_data)) {
        # Extract population frequencies
        joint_pop_afs <- extract_population_afs(variant_data$joint$populations)
        
        # Get frequencies
        joint_af <- variant_data$joint$af %||% 0
        genome_af <- variant_data$genome$af %||% 0
        exome_af <- variant_data$exome$af %||% 0
        
        # Calculate max AF across all datasets
        max_af <- max(joint_af, genome_af, exome_af, na.rm = TRUE)
        
        # Create result dataframe with all annotations
        results[[i]] <- data.frame(
          CHROM = variants_df$CHROM[i],
          POS = variants_df$POS[i],
          REF = ref,
          ALT = alt,
          
          # Joint frequencies (genome + exome combined - RECOMMENDED)
          gnomAD_AF = joint_af,
          gnomAD_AC = variant_data$joint$ac %||% 0,
          gnomAD_AN = variant_data$joint$an %||% 0,
          gnomAD_nhomalt = variant_data$joint$homozygote_count %||% 0,
          gnomAD_nhemialt = variant_data$joint$hemizygote_count %||% 0,
          
          # Separate genome data
          gnomAD_genome_AF = genome_af,
          gnomAD_genome_AC = variant_data$genome$ac %||% 0,
          gnomAD_genome_AN = variant_data$genome$an %||% 0,
          gnomAD_genome_nhomalt = variant_data$genome$homozygote_count %||% 0,
          
          # Separate exome data
          gnomAD_exome_AF = exome_af,
          gnomAD_exome_AC = variant_data$exome$ac %||% 0,
          gnomAD_exome_AN = variant_data$exome$an %||% 0,
          gnomAD_exome_nhomalt = variant_data$exome$homozygote_count %||% 0,
          
          # Maximum AF
          gnomAD_max_AF = max_af,
          
          # Population maximum (filtering AF at 95% confidence)
          gnomAD_AF_popmax = variant_data$joint$faf95$popmax %||% 0,
          gnomAD_popmax_population = variant_data$joint$faf95$popmax_population %||% NA,
          
          # Quality flags
          gnomAD_flags = paste(unlist(variant_data$flags), collapse = ";"),
          
          # Presence indicator
          in_gnomAD = TRUE,
          
          stringsAsFactors = FALSE
        ) %>%
          cbind(as.data.frame(joint_pop_afs))
        
      } else {
        # Variant not found in gnomAD
        zero_pops <- extract_population_afs(NULL)
        
        results[[i]] <- data.frame(
          CHROM = variants_df$CHROM[i],
          POS = variants_df$POS[i],
          REF = ref,
          ALT = alt,
          gnomAD_AF = 0,
          gnomAD_AC = 0,
          gnomAD_AN = 0,
          gnomAD_nhomalt = 0,
          gnomAD_nhemialt = 0,
          gnomAD_genome_AF = 0,
          gnomAD_genome_AC = 0,
          gnomAD_genome_AN = 0,
          gnomAD_genome_nhomalt = 0,
          gnomAD_exome_AF = 0,
          gnomAD_exome_AC = 0,
          gnomAD_exome_AN = 0,
          gnomAD_exome_nhomalt = 0,
          gnomAD_max_AF = 0,
          gnomAD_AF_popmax = 0,
          gnomAD_popmax_population = NA,
          gnomAD_flags = "",
          in_gnomAD = FALSE,
          stringsAsFactors = FALSE
        ) %>%
          cbind(as.data.frame(zero_pops))
      }
      
      # Rate limiting - be respectful to gnomAD servers
      Sys.sleep(0.05)
      
    }, error = function(e) {
      cat(sprintf("Error at variant %d: %s\n", i, e$message))
      
      zero_pops <- extract_population_afs(NULL)
      
      results[[i]] <- data.frame(
        CHROM = variants_df$CHROM[i],
        POS = variants_df$POS[i],
        REF = ref,
        ALT = alt,
        gnomAD_AF = NA,
        gnomAD_AC = NA,
        gnomAD_AN = NA,
        gnomAD_nhomalt = NA,
        gnomAD_nhemialt = NA,
        gnomAD_genome_AF = NA,
        gnomAD_genome_AC = NA,
        gnomAD_genome_AN = NA,
        gnomAD_genome_nhomalt = NA,
        gnomAD_exome_AF = NA,
        gnomAD_exome_AC = NA,
        gnomAD_exome_AN = NA,
        gnomAD_exome_nhomalt = NA,
        gnomAD_max_AF = NA,
        gnomAD_AF_popmax = NA,
        gnomAD_popmax_population = NA,
        gnomAD_flags = NA,
        in_gnomAD = NA,
        stringsAsFactors = FALSE
      ) %>%
        cbind(as.data.frame(lapply(zero_pops, function(x) NA)))
    })
    
    # Progress reporting and checkpointing
    if (i %% 100 == 0) {
      pct <- round(100 * i / total, 1)
      cat(sprintf("Progress: %d/%d (%s%%) - ETA: ~%.1f min\n", 
                  i, total, pct, (total - i) * 0.05 / 60))
      
      if (save_checkpoint) {
        saveRDS(list(results = results, last_idx = i), checkpoint_file)
      }
    }
  }
  
  cat("Combining results...\n")
  final_df <- bind_rows(results)
  
  # Clean up checkpoint file
  if (save_checkpoint && file.exists(checkpoint_file)) {
    file.remove(checkpoint_file)
    cat("Checkpoint file removed.\n")
  }
  
  return(final_df)
}

# ============================================================================
# FUNCTION: Complete annotation pipeline
# ============================================================================
run_complete_analysis <- function(input_file, 
                                  output_prefix = "gnomad_annotated") {
  
  cat("=== gnomAD v4 Annotation Pipeline ===\n\n")
  
  # Step 1: Load variants
  cat("Step 1: Loading variants...\n")
  
  # Auto-detect file format
  if (grepl("\\.vcf$|\\.vcf\\.gz$", input_file)) {
    # VCF format
    variants <- read.table(input_file, comment.char = "#", 
                           stringsAsFactors = FALSE)
    colnames(variants)[1:5] <- c("CHROM", "POS", "ID", "REF", "ALT")
  } else if (grepl("\\.csv$", input_file)) {
    variants <- read.csv(input_file, stringsAsFactors = FALSE)
  } else {
    variants <- read.table(input_file, header = TRUE, 
                           stringsAsFactors = FALSE)
  }
  
  cat(sprintf("Loaded %d variants\n\n", nrow(variants)))
  
  # Step 2: Annotate with gnomAD v4
  cat("Step 2: Annotating with gnomAD v4...\n")
  cat(sprintf("This will take approximately %.1f minutes\n\n", 
              nrow(variants) * 0.05 / 60))
  
  annotated <- annotate_variants_batch(variants)
  
  # Step 3: Save results
  output_file <- paste0(output_prefix, "_annotated.csv")
  write.csv(annotated, output_file, row.names = FALSE)
  cat(sprintf("\nSaved: %s\n", output_file))
  
  # Step 4: Summary statistics
  cat("\n=== Summary Statistics ===\n")
  cat(sprintf("Total variants: %d\n", nrow(annotated)))
  cat(sprintf("Found in gnomAD: %d (%.1f%%)\n", 
              sum(annotated$in_gnomAD, na.rm = TRUE),
              100 * sum(annotated$in_gnomAD, na.rm = TRUE) / nrow(annotated)))
  cat(sprintf("Not in gnomAD: %d (%.1f%%)\n", 
              sum(!annotated$in_gnomAD, na.rm = TRUE),
              100 * sum(!annotated$in_gnomAD, na.rm = TRUE) / nrow(annotated)))
  
  if (sum(annotated$in_gnomAD, na.rm = TRUE) > 0) {
    cat(sprintf("\nAllele Frequency Distribution (variants in gnomAD):\n"))
    cat(sprintf("  Mean AF: %.6f\n", 
                mean(annotated$gnomAD_AF[annotated$in_gnomAD], na.rm = TRUE)))
    cat(sprintf("  Median AF: %.6f\n", 
                median(annotated$gnomAD_AF[annotated$in_gnomAD], na.rm = TRUE)))
    cat(sprintf("  Rare (AF < 0.01): %d (%.1f%%)\n",
                sum(annotated$gnomAD_AF < 0.01 & annotated$in_gnomAD, na.rm = TRUE),
                100 * sum(annotated$gnomAD_AF < 0.01 & annotated$in_gnomAD, na.rm = TRUE) / 
                  sum(annotated$in_gnomAD, na.rm = TRUE)))
    cat(sprintf("  Common (AF >= 0.05): %d (%.1f%%)\n",
                sum(annotated$gnomAD_AF >= 0.05 & annotated$in_gnomAD, na.rm = TRUE),
                100 * sum(annotated$gnomAD_AF >= 0.05 & annotated$in_gnomAD, na.rm = TRUE) / 
                  sum(annotated$in_gnomAD, na.rm = TRUE)))
  }
  
  cat("\n=== Analysis Complete ===\n")
  
  return(annotated)
}

# ============================================================================
# USAGE EXAMPLES
# ============================================================================
# 
# Example 1: Annotate a VCF file
# result <- run_complete_analysis("my_variants.vcf", output_prefix = "my_output")
#
# Example 2: Annotate a CSV file
# result <- run_complete_analysis("my_variants.csv", output_prefix = "my_output")
#
# Example 3: Just annotate without the full pipeline
# variants_df <- read.csv("my_variants.csv")
# annotated <- annotate_variants_batch(variants_df)
#
# ============================================================================