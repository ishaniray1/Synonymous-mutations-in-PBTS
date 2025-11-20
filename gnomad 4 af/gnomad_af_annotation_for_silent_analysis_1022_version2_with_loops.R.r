# ----------------------------------
# FUNCTION Extract population allele frequencies
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

# ----------------------------------------

# FUNCTION: Annotate variants with retry logic for missing variants
annotate_variants_batch_v2 <- function(variants_df, 
                                    batch_size = 100,
                                    max_retries = 7,
                                    save_checkpoint = TRUE,
                                    checkpoint_file = "gnomad_checkpoint.rds") {
  
  # Load existing checkpoint
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
  
  # Initial pass through all variants
  cat("=== INITIAL ANNOTATION PASS ===\n")
  for (i in start_idx:total) {
    chrom <- gsub("chr", "", variants_df$CHROM[i])
    pos <- variants_df$POS[i]
    ref <- variants_df$REF[i]
    alt <- variants_df$ALT[i]
    
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
      
      if (!is.null(data$errors)) {
        stop(paste("API Error:", data$errors[[1]]$message))
      }
      
      variant_data <- data$data$variant
      
      if (!is.null(variant_data)) {
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
        
        total_ac <- (variant_data$genome$ac %||% 0) + (variant_data$exome$ac %||% 0)
        total_an <- (variant_data$genome$an %||% 0) + (variant_data$exome$an %||% 0)
        joint_af <- if (total_an > 0) total_ac / total_an else 0
        
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
          joint_pop_afs <- extract_population_afs(NULL)
        }
        
        results[[i]] <- data.frame(
          CHROM = variants_df$CHROM[i],
          POS = variants_df$POS[i],
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
        
      } else {
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
    
    if (i %% batch_size == 0) {
      pct <- round(100 * i / total, 1)
      elapsed <- (i - start_idx + 1) * 0.05
      remaining <- (total - i) * 0.05
      cat(sprintf("Progress: %d/%d (%s%%) - ETA: ~%.1f min\n", 
                  i, total, pct, remaining / 60))
      
      if (save_checkpoint) {
        saveRDS(list(results = results, last_idx = i), checkpoint_file)
      }
    }
  }
  
  cat("\n=== INITIAL PASS COMPLETE ===\n")
  cat("Combining initial results...\n")
  final_df <- bind_rows(results)
  
  # RETRY LOGIC: Re-query variants with AF=0
  for (retry_num in 1:max_retries) {
    cat(sprintf("\n=== RETRY PASS %d/%d ===\n", retry_num, max_retries))
    
    # Find variants with AF = 0 (not found)
    retry_variants <- final_df %>%
      filter(gnomAD_AF == 0 & !is.na(gnomAD_AF))
    
    if (nrow(retry_variants) == 0) {
      cat("No variants with AF=0 to retry. Done!\n")
      break
    }
    
    cat(sprintf("Retrying %d variants with AF=0...\n", nrow(retry_variants)))
    
    retry_results <- list()
    
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
          
          # Only update if we found non-zero AF
          if (joint_af > 0) {
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
              joint_pop_afs <- extract_population_afs(NULL)
            }
            
            retry_results[[j]] <- data.frame(
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
          }
        }
        
        Sys.sleep(0.05)
        
      }, error = function(e) {
        # Silent error handling for retries
      })
      
      if (j %% batch_size == 0) {
        cat(sprintf("  Retry progress: %d/%d\n", j, nrow(retry_variants)))
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
      
      cat(sprintf("  Updated %d variants with new data\n", nrow(retry_df)))
    }
  }
  
  # Clean up checkpoint
  if (save_checkpoint && file.exists(checkpoint_file)) {
    file.remove(checkpoint_file)
    cat("\nCheckpoint file removed.\n")
  }
  
  cat("\n=== ANNOTATION COMPLETE ===\n")
  cat(sprintf("Total variants: %d\n", nrow(final_df)))
  cat(sprintf("Found in gnomAD: %d\n", sum(final_df$gnomAD_AF > 0, na.rm = TRUE)))
  cat(sprintf("Not found (AF=0): %d\n", sum(final_df$gnomAD_AF == 0, na.rm = TRUE)))
  
  return(final_df)
}

# Load file
file <- read_csv("silent_variants_df.csv")
colnames(file)

# copy with renamed columns for annotation
file_renamed <- file %>%
  mutate(
    CHROM = Chromosome,
    POS = Start_Position,
    REF = Reference_Allele,
    ALT = Tumor_Seq_Allele2
  )
# Verify req columns 
required_cols <- c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2")
if (!all(required_cols %in% colnames(file))) {
  stop("Missing required columns!")
}

file_renamed2 <-  file_renamed%>% 
  distinct(CHROM, POS, REF, ALT, .keep_all = TRUE) %>%  # keep only unique combos
  select(CHROM, POS, REF, ALT, Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2)

# Run annotation 
annotated_again_v2 <- annotate_variants_batch_v2(file_renamed2)
annotated_again_v2_backup <- annotated_again_v2
write.csv(annotated_again_v2, "annotated_only_v2.csv", row.names = FALSE)
