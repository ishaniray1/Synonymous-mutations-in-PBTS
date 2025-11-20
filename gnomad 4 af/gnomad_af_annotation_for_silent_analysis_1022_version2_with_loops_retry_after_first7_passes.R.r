# FUNCTION: Retry annotation for missing variants from existing results
retry_missing_variants <- function(existing_results,
                                   max_retries = 2,
                                   output_file = NULL) {
  
  # Load existing results (can be dataframe, file path, or variable name)
  cat("Loading existing results...\n")
  if (is.data.frame(existing_results)) {
    # Already a dataframe
    final_df <- existing_results
  } else if (is.character(existing_results)) {
    # File path
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
  missing_count <- sum(final_df$gnomAD_4_AF == 0 | is.na(final_df$gnomAD_4_AF), na.rm = TRUE)
  cat(sprintf("Found %d unannotated variants\n\n", missing_count))
  
  if (missing_count == 0) {
    cat("No missing variants to retry!\n")
    return(final_df)
  }
  
  # RETRY LOGIC
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
    
    retry_results <- list()
    successful_updates <- 0
    
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
          
          # Only update if we found non-zero AF OR valid AN (means variant exists but is rare)
          if (total_an > 0) {
            max_af <- max(joint_af, genome_af, exome_af, na.rm = TRUE)
            
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
              # Create empty population AFs
              pops <- c("afr", "amr", "asj", "eas", "fin", "nfe", "sas", "mid", "ami", "remaining")
              for (pop_name in pops) {
                joint_pop_afs[[paste0("gnomAD_", toupper(pop_name), "_AF")]] <- 0
              }
            }
            
            retry_results[[length(retry_results) + 1]] <- data.frame(
              CHROM = retry_variants$CHROM[j],
              POS = retry_variants$POS[j],
              REF = ref,
              ALT = alt,
              gnomAD_AF = joint_af,
              gnomAD_AC = total_ac,
              gnomAD_AN = total_an,
              gnomAD_nhomalt = (variant_data$genome$homozygote_count %||% 0) + 
                (variant_data$exome$homozygote_count %||% 0),
              gnomAD_nhemialt = variant_data$genome$hemizygote_count %||% 0,
              gnomAD_genome_AF = genome_af,
              gnomAD_genome_AC = variant_data$genome$ac %||% 0,
              gnomAD_genome_AN = variant_data$genome$an %||% 0,
              gnomAD_genome_nhomalt = variant_data$genome$homozygote_count %||% 0,
              gnomAD_exome_AF = exome_af,
              gnomAD_exome_AC = variant_data$exome$ac %||% 0,
              gnomAD_exome_AN = variant_data$exome$an %||% 0,
              gnomAD_exome_nhomalt = variant_data$exome$homozygote_count %||% 0,
              gnomAD_max_AF = max_af,
              gnomAD_AF_popmax = variant_data$genome$faf95$popmax %||% 0,
              gnomAD_popmax_population = variant_data$genome$faf95$popmax_population %||% NA,
              gnomAD_flags = paste(unlist(variant_data$flags), collapse = ";"),
              in_gnomAD = TRUE,
              stringsAsFactors = FALSE
            ) %>%
              cbind(as.data.frame(joint_pop_afs))
            
            successful_updates <- successful_updates + 1
          }
        }
        
        Sys.sleep(1)
        
      }, error = function(e) {
        # Silent error handling
      })
      
      if (j %% 100 == 0) {
        cat(sprintf("  Progress: %d/%d (found %d so far)\n", 
                    j, nrow(retry_variants), successful_updates))
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
      
      cat(sprintf("  ✓ Updated %d variants with new data\n", nrow(retry_df)))
    } else {
      cat("  No new variants found in this pass\n")
    }
    
    # Calculate remaining
    remaining <- sum(final_df$gnomAD_AF == 0 | is.na(final_df$gnomAD_AF), na.rm = TRUE)
    cat(sprintf("  Remaining unannotated: %d\n\n", remaining))
    
    # Save intermediate results
    if (!is.null(output_file)) {
      if (grepl("\\.rds$", output_file)) {
        saveRDS(final_df, output_file)
      } else {
        write.csv(final_df, output_file, row.names = FALSE)
      }
      cat(sprintf("  Saved intermediate results to %s\n\n", output_file))
    }
  }
  
  # Final summary
  cat("\n=== RETRY COMPLETE ===\n")
  cat(sprintf("Total variants: %d\n", nrow(final_df)))
  cat(sprintf("Found in gnomAD (AF>0): %d\n", sum(final_df$gnomAD_AF > 0, na.rm = TRUE)))
  cat(sprintf("In gnomAD but rare (AF=0, AN>0): %d\n", 
              sum(final_df$gnomAD_AF == 0 & final_df$gnomAD_AN > 0, na.rm = TRUE)))
  cat(sprintf("Still missing: %d\n", sum(final_df$gnomAD_AF == 0 & final_df$gnomAD_AN == 0, na.rm = TRUE)))
  
  return(final_df)
}



annotated_again_v2_updated <- retry_missing_variants(annotated_again_v2, max_retries = 10)
write.csv(annotated_again_v2_updated, "annotated_only_v2_updated.csv", row.names = FALSE)

annotated_again_v2_updated_2 <- retry_missing_variants(annotated_again_v2_updated, max_retries = 50)
write.csv(annotated_again_v2_updated_2, "annotated_only_v2_updated_2.csv", row.names = FALSE)
colnames(annotated_again_v2_updated_2)

annotated_again_v2_updated_2 <- annotated_again_v2_updated_2 %>%
  rename(
    gnomAD_4_AF = gnomAD_AF,
    gnomAD_4_AC = gnomAD_AC,
    gnomAD_4_AN = gnomAD_AN,
    gnomAD_4_nhomalt = gnomAD_nhomalt,
    gnomAD_4_nhemialt = gnomAD_nhemialt,
    gnomAD_4_genome_AF = gnomAD_genome_AF,
    gnomAD_4_genome_AC = gnomAD_genome_AC,
    gnomAD_4_genome_AN = gnomAD_genome_AN,
    gnomAD_4_genome_nhomalt = gnomAD_genome_nhomalt,
    gnomAD_4_exome_AF = gnomAD_exome_AF,
    gnomAD_4_exome_AC = gnomAD_exome_AC,
    gnomAD_4_exome_AN = gnomAD_exome_AN,
    gnomAD_4_exome_nhomalt = gnomAD_exome_nhomalt,
    gnomAD_4_max_AF = gnomAD_max_AF,
    gnomAD_4_AF_popmax = gnomAD_AF_popmax,
    gnomAD_4_popmax_population = gnomAD_popmax_population,
    gnomAD_4_flags = gnomAD_flags,
    gnomAD_4_in_gnomAD = in_gnomAD,
    gnomAD_4_AFR_AF = gnomAD_AFR_AF,
    gnomAD_4_AMR_AF = gnomAD_AMR_AF,
    gnomAD_4_ASJ_AF = gnomAD_ASJ_AF,
    gnomAD_4_EAS_AF = gnomAD_EAS_AF,
    gnomAD_4_FIN_AF = gnomAD_FIN_AF,
    gnomAD_4_NFE_AF = gnomAD_NFE_AF,
    gnomAD_4_SAS_AF = gnomAD_SAS_AF,
    gnomAD_4_MID_AF = gnomAD_MID_AF,
    gnomAD_4_AMI_AF = gnomAD_AMI_AF,
    gnomAD_4_REMAINING_AF = gnomAD_REMAINING_AF
  )
colnames(annotated_again_v2_updated_2)
write.csv(annotated_again_v2_updated_2, "annotated_only_v2_updated_2.csv", row.names = FALSE)

final_file_again_v2 <- file_renamed2 %>%
  left_join(
    annotated_again_v2_updated_2,
    by = c("Chromosome" = "CHROM", 
           "Start_Position" = "POS",
           "Reference_Allele" = "REF",
           "Tumor_Seq_Allele2" = "ALT"),
    relationship = "many-to-many"
  )
write.csv(final_file_again_v2, "final_file_again_v2_deduplicated.csv", row.names = FALSE)

final_file_with_all_rows <- file_renamed %>%
  left_join(
    annotated_again_v2_updated_2,
    by = c("Chromosome" = "CHROM", 
           "Start_Position" = "POS",
           "Reference_Allele" = "REF",
           "Tumor_Seq_Allele2" = "ALT"),
    relationship = "many-to-many"
  )

write.csv(final_file_with_all_rows, "final_file_with_all_rows.csv", row.names = FALSE)

colnames(final_file_with_all_rows)

comparison_stats <- data.frame(
  Version = c("gnomAD 3.1.1", "gnomAD 4"),
  Total_Rows = c(nrow(final_file_with_all_rows), nrow(final_file_with_all_rows)),
  NA_Count = c(
    sum(is.na(final_file_with_all_rows$gnomad_3_1_1_AF)),
    sum(is.na(final_file_with_all_rows$gnomAD_4_AF))
  ),
  Zero_Count = c(
    sum(final_file_with_all_rows$gnomad_3_1_1_AF == 0, na.rm = TRUE),
    sum(final_file_with_all_rows$gnomAD_4_AF == 0, na.rm = TRUE)
  )
)

comparison_stats <- comparison_stats %>%
  mutate(
    NA_Percentage = round((NA_Count / Total_Rows) * 100, 2),
    Zero_Percentage = round((Zero_Count / Total_Rows) * 100, 2),
    Non_NA_Non_Zero = Total_Rows - NA_Count - Zero_Count,
    Non_NA_Non_Zero_Percentage = round((Non_NA_Non_Zero / Total_Rows) * 100, 2)
  )

print(comparison_stats)

#------------

annotated_again_v2_updated_1030 <- retry_missing_variants(annotated_again_v2_updated_2, max_retries = 2)
write.csv(annotated_again_v2_updated_1030, "annotated_again_v2_updated_1030.csv", row.names = FALSE)











# USAGE EXAMPLES:
# 
# # Option 1: Pass a dataframe directly
# updated_results <- retry_missing_variants(my_df, max_retries = 10)
# 
# # Option 2: Load from RDS file
# updated_results <- retry_missing_variants("my_results.rds", max_retries = 10)
# 
# # Option 3: Load from CSV and save intermediate results
# updated_results <- retry_missing_variants(
#   "my_results.csv", 
#   max_retries = 10,
#   output_file = "updated_results.rds"
# )
# 
# # Then save final results
# saveRDS(updated_results, "final_annotated.rds")
# write.csv(updated_results, "final_annotated.csv", row.names = FALSE)