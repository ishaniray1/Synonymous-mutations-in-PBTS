# FUNCTION: Robust retry with timeout handling and delays
retry_missing_variants_robust <- function(existing_results,
                                          max_retries = 5,
                                          delay_between_queries = 2,
                                          batch_size = 50,
                                          batch_delay = 30,
                                          output_file = NULL) {
  
  library(httr)
  library(dplyr)
  
  # Load existing results
  cat("Loading existing results...\n")
  if (is.data.frame(existing_results)) {
    final_df <- existing_results
  } else if (is.character(existing_results)) {
    if (grepl("\\.rds$", existing_results)) {
      final_df <- readRDS(existing_results)
    } else if (grepl("\\.csv$", existing_results)) {
      final_df <- read.csv(existing_results, stringsAsFactors = FALSE)
    } else {
      stop("File must be .rds or .csv")
    }
  } else {
    stop("Input must be a dataframe or file path")
  }
  
  cat(sprintf("Loaded %d total variants\n", nrow(final_df)))
  
  # Find missing variants (AF = 0 or NA)
  missing_mask <- final_df$gnomAD_4_AF == 0 | is.na(final_df$gnomAD_4_AF)
  missing_count <- sum(missing_mask, na.rm = TRUE)
  cat(sprintf("Found %d unannotated variants\n\n", missing_count))
  
  if (missing_count == 0) {
    cat("No missing variants to retry!\n")
    return(final_df)
  }
  
  # RETRY LOGIC with batching and delays
  for (retry_num in 1:max_retries) {
    cat(sprintf("=== RETRY PASS %d/%d ===\n", retry_num, max_retries))
    
    # Find variants with AF = 0 or NA
    retry_variants <- final_df %>%
      filter(gnomAD_4_AF == 0 | is.na(gnomAD_4_AF))
    
    if (nrow(retry_variants) == 0) {
      cat("All variants annotated! Done!\n")
      break
    }
    
    cat(sprintf("Retrying %d missing variants...\n", nrow(retry_variants)))
    cat(sprintf("Using %d second delay between queries\n", delay_between_queries))
    cat(sprintf("Taking %d second break every %d queries\n\n", batch_delay, batch_size))
    
    retry_results <- list()
    successful_updates <- 0
    failed_count <- 0
    
    for (j in 1:nrow(retry_variants)) {
      chrom <- gsub("chr", "", retry_variants$CHROM[j])
      pos <- retry_variants$POS[j]
      ref <- retry_variants$REF[j]
      alt <- retry_variants$ALT[j]
      
      query <- sprintf('{
        variant(variantId: "%s-%d-%s-%s", dataset: gnomad_r4) {
          variantId
          genome {
            ac
            an
            homozygote_count
            hemizygote_count
            populations {
              id
              ac
              an
              homozygote_count
            }
            faf95 {
              popmax
              popmax_population
            }
          }
          exome {
            ac
            an
            homozygote_count
            populations {
              id
              ac
              an
              homozygote_count
            }
          }
          flags
        }
      }', chrom, pos, ref, alt)
      
      success <- FALSE
      attempt <- 1
      max_attempts <- 3
      
      while (!success && attempt <= max_attempts) {
        tryCatch({
          response <- POST(
            "https://gnomad.broadinstitute.org/api",
            body = list(query = query),
            encode = "json",
            timeout(15)
          )
          
          if (status_code(response) == 200) {
            data <- content(response, "parsed")
            variant_data <- data$data$variant
            
            if (!is.null(variant_data)) {
              genome_af <- 0
              exome_af <- 0
              
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
              
              total_ac <- (variant_data$genome$ac %||% 0) + (variant_data$exome$ac %||% 0)
              total_an <- (variant_data$genome$an %||% 0) + (variant_data$exome$an %||% 0)
              joint_af <- if (total_an > 0) total_ac / total_an else 0
              
              # Only update if we found valid data
              if (total_an > 0) {
                max_af <- max(joint_af, genome_af, exome_af, na.rm = TRUE)
                
                # Population AFs
                joint_pop_afs <- list()
                if (!is.null(variant_data$genome$populations)) {
                  pops <- c("afr", "amr", "asj", "eas", "fin", "nfe", "sas", "mid", "ami", "remaining")
                  for (pop_name in pops) {
                    pop_data <- Filter(function(x) x$id == pop_name, variant_data$genome$populations)
                    if (length(pop_data) > 0) {
                      pop_ac <- pop_data[[1]]$ac %||% 0
                      pop_an <- pop_data[[1]]$an %||% 0
                      pop_af <- if (pop_an > 0) pop_ac / pop_an else 0
                      joint_pop_afs[[paste0("gnomAD_4_", toupper(pop_name), "_AF")]] <- pop_af
                    } else {
                      joint_pop_afs[[paste0("gnomAD_4_", toupper(pop_name), "_AF")]] <- 0
                    }
                  }
                } else {
                  pops <- c("afr", "amr", "asj", "eas", "fin", "nfe", "sas", "mid", "ami", "remaining")
                  for (pop_name in pops) {
                    joint_pop_afs[[paste0("gnomAD_4_", toupper(pop_name), "_AF")]] <- 0
                  }
                }
                
                retry_results[[length(retry_results) + 1]] <- data.frame(
                  CHROM = retry_variants$CHROM[j],
                  POS = retry_variants$POS[j],
                  REF = ref,
                  ALT = alt,
                  gnomAD_4_AF = joint_af,
                  gnomAD_4_AC = total_ac,
                  gnomAD_4_AN = total_an,
                  gnomAD_4_nhomalt = (variant_data$genome$homozygote_count %||% 0) + 
                    (variant_data$exome$homozygote_count %||% 0),
                  gnomAD_4_nhemialt = variant_data$genome$hemizygote_count %||% 0,
                  gnomAD_4_genome_AF = genome_af,
                  gnomAD_4_genome_AC = variant_data$genome$ac %||% 0,
                  gnomAD_4_genome_AN = variant_data$genome$an %||% 0,
                  gnomAD_4_genome_nhomalt = variant_data$genome$homozygote_count %||% 0,
                  gnomAD_4_exome_AF = exome_af,
                  gnomAD_4_exome_AC = variant_data$exome$ac %||% 0,
                  gnomAD_4_exome_AN = variant_data$exome$an %||% 0,
                  gnomAD_4_exome_nhomalt = variant_data$exome$homozygote_count %||% 0,
                  gnomAD_4_max_AF = max_af,
                  gnomAD_4_AF_popmax = variant_data$genome$faf95$popmax %||% 0,
                  gnomAD_4_popmax_population = variant_data$genome$faf95$popmax_population %||% NA,
                  gnomAD_4_flags = paste(unlist(variant_data$flags), collapse = ";"),
                  gnomAD_4_in_gnomAD = TRUE,
                  stringsAsFactors = FALSE
                ) %>%
                  cbind(as.data.frame(joint_pop_afs))
                
                successful_updates <- successful_updates + 1
              }
            }
            success <- TRUE
          }
          
        }, error = function(e) {
          if (attempt < max_attempts) {
            cat(sprintf("  Attempt %d failed for variant %s-%d-%s-%s, retrying...\n", 
                        attempt, chrom, pos, ref, alt))
            Sys.sleep(5)
          } else {
            failed_count <- failed_count + 1
          }
        })
        
        attempt <- attempt + 1
      }
      
      # Regular delay between queries
      Sys.sleep(delay_between_queries)
      
      # Batch delay - take a longer break
      if (j %% batch_size == 0) {
        cat(sprintf("  Progress: %d/%d (found %d, failed %d)\n", 
                    j, nrow(retry_variants), successful_updates, failed_count))
        cat(sprintf("  Taking %d second break...\n", batch_delay))
        Sys.sleep(batch_delay)
      }
    }
    
    # Update final_df with successful retries
    if (length(retry_results) > 0) {
      retry_df <- bind_rows(retry_results)
      
      # Remove old entries and add updated ones
      final_df <- final_df %>%
        anti_join(retry_df, by = c("CHROM", "POS", "REF", "ALT")) %>%
        bind_rows(retry_df) %>%
        arrange(CHROM, POS)
      
      cat(sprintf("\n  ✓ Updated %d variants with new data\n", nrow(retry_df)))
    } else {
      cat("\n  No new variants found in this pass\n")
    }
    
    # Calculate remaining
    remaining <- sum(final_df$gnomAD_4_AF == 0 | is.na(final_df$gnomAD_4_AF), na.rm = TRUE)
    cat(sprintf("  Remaining unannotated: %d\n", remaining))
    cat(sprintf("  Failed queries: %d\n\n", failed_count))
    
    # Save intermediate results
    if (!is.null(output_file)) {
      if (grepl("\\.rds$", output_file)) {
        saveRDS(final_df, output_file)
      } else {
        write.csv(final_df, output_file, row.names = FALSE)
      }
      cat(sprintf("  Saved intermediate results to %s\n\n", output_file))
    }
    
    # If no improvement, maybe stop
    if (length(retry_results) == 0 && retry_num > 2) {
      cat("No progress in this pass. Consider stopping or adjusting parameters.\n")
    }
  }
  
  # Final summary
  cat("\n=== RETRY COMPLETE ===\n")
  cat(sprintf("Total variants: %d\n", nrow(final_df)))
  cat(sprintf("Found in gnomAD (AF>0): %d\n", sum(final_df$gnomAD_4_AF > 0, na.rm = TRUE)))
  cat(sprintf("In gnomAD but rare (AF=0, AN>0): %d\n", 
              sum(final_df$gnomAD_4_AF == 0 & final_df$gnomAD_4_AN > 0, na.rm = TRUE)))
  cat(sprintf("Still missing: %d\n", 
              sum((final_df$gnomAD_4_AF == 0 | is.na(final_df$gnomAD_4_AF)) & 
                    (final_df$gnomAD_4_AN == 0 | is.na(final_df$gnomAD_4_AN)), na.rm = TRUE)))
  
  return(final_df)
}

# USAGE:
# Filter for missing variants 
annotated_again_v2_updated_2 <- read.csv("annotated_only_v2_updated_2.csv", header = TRUE, stringsAsFactors = FALSE)

missing_variants <- annotated_again_v2_updated_2 %>%
  filter(gnomAD_4_AF == 0 | is.na(gnomAD_4_AF))

cat(sprintf("Found %d missing variants to retry\n", nrow(missing_variants)))


annotated_v2_final_p1 <- read.csv("annotated_v2_final_p1.csv", header = TRUE, stringsAsFactors = FALSE)
missing_variants <- annotated_v2_final_p1 %>%
  filter(gnomAD_4_AF == 0 | is.na(gnomAD_4_AF))

cat(sprintf("Found %d missing variants to retry\n", nrow(missing_variants)))

# Run retry with conservative settings
final_results <- retry_missing_variants_robust(
  existing_results = missing_variants,
  max_retries = 1,              # Number of complete passes
  delay_between_queries = 3,     # 3 seconds between each query
  batch_size = 50,               # Take a break every 50 queries
  batch_delay = 60,              # 60 second break between batches
  output_file = "annotated_v2_final.csv"
)

silents_variants_df <- read.csv("silent_variants_df.csv", header = TRUE, stringsAsFactors = FALSE)
missing_only <- read.csv("annotated_v2_final.csv", header = TRUE, stringsAsFactors = FALSE)


# Save final results
write.csv(final_results, "annotated_v2_complete.csv", row.names = FALSE)