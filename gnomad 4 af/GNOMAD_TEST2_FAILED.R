# ============================================================================
# FUNCTION: Extract population allele frequencies
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
    
    # FIXED GraphQL query for gnomAD v4 - Note the changes!
    # 1. variantId instead of variant_id
    # 2. dataset is separate parameter
    # 3. Different field structure
    query <- sprintf('{
      variant(variantId: "%s-%d-%s-%s", dataset: gnomad_r4) {
        variantId
        
        genome {
          ac
          an
          homozygoteCount
          hemizygoteCount
          populations {
            id
            ac
            an
            homozygoteCount
          }
          faf95 {
            popmax
            popmaxPopulation
          }
        }
        
        exome {
          ac
          an
          homozygoteCount
          populations {
            id
            ac
            an
            homozygoteCount
          }
        }
        
        flags
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
      
      # Check for errors
      if (!is.null(data$errors)) {
        stop(paste("API Error:", data$errors[[1]]$message))
      }
      
      variant_data <- data$data$variant
      
      if (!is.null(variant_data)) {
        # Calculate allele frequencies from AC/AN
        genome_af <- 0
        exome_af <- 0
        joint_af <- 0
        
        if (!is.null(variant_data$genome)) {
          genome_ac <- variant_data$genome$ac %||% 0
          genome_an <- variant_data$genome$an %||% 0
          genome_af <- if (genome_an > 0) genome_ac / genome_an else 0
        }
        
        if (!is.null(variant_data$exome)) {
          exome_ac <- variant_data$exome$ac %||% 0
          exome_an <- variant_data$exome$an %||% 0
          exome_af <- if (exome_an > 0) exome_ac / exome_an else 0
        }
        
        # Joint AF calculation (combined genome + exome)
        total_ac <- (variant_data$genome$ac %||% 0) + (variant_data$exome$ac %||% 0)
        total_an <- (variant_data$genome$an %||% 0) + (variant_data$exome$an %||% 0)
        joint_af <- if (total_an > 0) total_ac / total_an else 0
        
        # Calculate max AF
        max_af <- max(joint_af, genome_af, exome_af, na.rm = TRUE)
        
        # Extract population frequencies from genome data
        joint_pop_afs <- list()
        if (!is.null(variant_data$genome$populations)) {
          pops <- c("afr", "amr", "asj", "eas", "fin", "nfe", "sas", "mid", "ami", "remaining")
          for (pop_name in pops) {
            pop_data <- Filter(function(x) x$id == pop_name, variant_data$genome$populations)
            if (length(pop_data) > 0) {
              pop_ac <- pop_data[[1]]$ac %||% 0
              pop_an <- pop_data[[1]]$an %||% 0
              pop_af <- if (pop_an > 0) pop_ac / pop_an else 0
              joint_pop_afs[[paste0("gnomAD_", toupper(pop_name), "_AF")]] <- pop_af
            } else {
              joint_pop_afs[[paste0("gnomAD_", toupper(pop_name), "_AF")]] <- 0
            }
          }
        } else {
          joint_pop_afs <- extract_population_afs(NULL)
        }
        
        # Create result dataframe
        results[[i]] <- data.frame(
          CHROM = variants_df$CHROM[i],
          POS = variants_df$POS[i],
          REF = ref,
          ALT = alt,
          
          # Joint frequencies (calculated from genome + exome)
          gnomAD_AF = joint_af,
          gnomAD_AC = total_ac,
          gnomAD_AN = total_an,
          gnomAD_nhomalt = (variant_data$genome$homozygoteCount %||% 0) + 
            (variant_data$exome$homozygoteCount %||% 0),
          gnomAD_nhemialt = variant_data$genome$hemizygoteCount %||% 0,
          
          # Separate genome data
          gnomAD_genome_AF = genome_af,
          gnomAD_genome_AC = variant_data$genome$ac %||% 0,
          gnomAD_genome_AN = variant_data$genome$an %||% 0,
          gnomAD_genome_nhomalt = variant_data$genome$homozygoteCount %||% 0,
          
          # Separate exome data
          gnomAD_exome_AF = exome_af,
          gnomAD_exome_AC = variant_data$exome$ac %||% 0,
          gnomAD_exome_AN = variant_data$exome$an %||% 0,
          gnomAD_exome_nhomalt = variant_data$exome$homozygoteCount %||% 0,
          
          # Maximum AF
          gnomAD_max_AF = max_af,
          
          # Population maximum (filtering AF at 95% confidence)
          gnomAD_AF_popmax = variant_data$genome$faf95$popmax %||% 0,
          gnomAD_popmax_population = variant_data$genome$faf95$popmaxPopulation %||% NA,
          
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
  
  cat("=== gnomAD v4 Annotation Pipeline (FIXED) ===\n\n")
  
  # Step 1: Load variants
  cat("Step 1: Loading variants...\n")
  
  # Auto-detect file format
  if (grepl("\\.vcf$|\\.vcf\\.gz$", input_file)) {
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
# USAGE
# ============================================================================
# result <- run_complete_analysis("my_variants.vcf", output_prefix = "my_output")
# ============================================================================

# Test with first 5 variants
test_df <- head(file_renamed, 5)
test_result <- annotate_variants_batch(test_df)
print(test_result)
