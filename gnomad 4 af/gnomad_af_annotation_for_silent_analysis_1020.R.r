# gnomAD v4 Variant Annotation Script
# ----------------------
# annotates variants with allele frequency data from gnomAD v4 
# including population-specific frequencies, homozygote counts, and quality flags
# -----------------------------

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
library(readr)

# null coalescing 
`%||%` <- function(x, y) if (is.null(x)) y else x

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
# FUNCTION Annotate variants in batches with checkpoint (ETA is wrong)
annotate_variants_batch <- function(variants_df, 
                                    batch_size = 50,
                                    save_checkpoint = TRUE,
                                    checkpoint_file = "gnomad_checkpoint.rds") {
  
  #  existing checkpoint
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
    
    # FIXED GraphQL query for gnomAD v4
    # variantId instead of variant_id
    # dataset is separate parameter
    # Different field structure
    # some others will note later
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
      
      # Check for errors
      if (!is.null(data$errors)) {
        stop(paste("API Error:", data$errors[[1]]$message))
      }
      
      variant_data <- data$data$variant
      
      if (!is.null(variant_data)) {
        #  allele frequencies from AC/AN
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
        
        # Joint AF  (combined genome + exome)
        total_ac <- (variant_data$genome$ac %||% 0) + (variant_data$exome$ac %||% 0)
        total_an <- (variant_data$genome$an %||% 0) + (variant_data$exome$an %||% 0)
        joint_af <- if (total_an > 0) total_ac / total_an else 0
        
        #  max AF
        max_af <- max(joint_af, genome_af, exome_af, na.rm = TRUE)
        
        # population frequencies from genome data
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
        
        # result dataframe
        results[[i]] <- data.frame(
          CHROM = variants_df$CHROM[i],
          POS = variants_df$POS[i],
          REF = ref,
          ALT = alt,
          
          # Joint frequencies ( from genome + exome)
          gnomAD_AF = joint_af,
          gnomAD_AC = total_ac,
          gnomAD_AN = total_an,
          gnomAD_nhomalt = (variant_data$genome$homozygote_count %||% 0) + 
            (variant_data$exome$homozygote_count %||% 0),
          gnomAD_nhemialt = variant_data$genome$hemizygote_count %||% 0,
          
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
          gnomAD_AF_popmax = variant_data$genome$faf95$popmax %||% 0,
          gnomAD_popmax_population = variant_data$genome$faf95$popmax_population %||% NA,
          
          # Quality 
          gnomAD_flags = paste(unlist(variant_data$flags), collapse = ";"),
          
          # Presence 
          in_gnomAD = TRUE,
          
          stringsAsFactors = FALSE
        ) %>%
          cbind(as.data.frame(joint_pop_afs))
        
      } else {
        # Variant not in gnomAD
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
      
      # Rate limiting 
      Sys.sleep(0.20)
      
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
    
    # Progress reporting, checkpoint
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
  
  # Clean checkpoint file
  if (save_checkpoint && file.exists(checkpoint_file)) {
    file.remove(checkpoint_file)
    cat("Checkpoint file removed.\n")
  }
  
  return(final_df)
}

# --------
# USAGE 
# variants_df <- read.csv("my_variants.csv")
# annotated <- annotate_variants_batch(variants_df)

# -----------

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
  dplyr::select(CHROM, POS, REF, ALT, Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2)

# Run annotation 
annotated_again <- annotate_variants_batch(file_renamed2)

annotated_subset <- annotated_again %>%
  filter(is.na(gnomAD_AF) | gnomAD_AF == 0)

annotated_subset <- annotated_again %>%
  filter(is.na(gnomAD_AF) | gnomAD_AF == 0) %>%
  select(CHROM, POS, REF, ALT)

annotated_again <- annotate_variants_batch(annotated_subset)


# Merge  to original 
final_file_again <- file %>%
  left_join(
    annotated_eg,
    by = c("Chromosome" = "CHROM", 
           "Start_Position" = "POS",
           "Reference_Allele" = "REF",
           "Tumor_Seq_Allele2" = "ALT"),
    suffix = c("", ".from_df1"),
    relationship = "many-to-many"
  )






write.csv(annotated, "annotated_only", 
          row.names = FALSE)
write.csv(annotated, "annotated_only.csv", row.names = FALSE)


write.csv(final_file, "silent_variants_df_with_gnomad_v4.csv", 
            row.names = FALSE)
# OR  as tab-delimited
write.table(final_file, "silent_variants_with_gnomad_v4.txt", 
             sep = "\t", quote = FALSE, row.names = FALSE)


#--Misc--
# How many were found?
table(annotated$in_gnomAD)

# some that were found
found_variants <- annotated[annotated$in_gnomAD == TRUE, ]
head(found_variants[, c("CHROM", "POS", "REF", "ALT", "gnomAD_AF", "gnomAD_AC")])
# Should see frequencies like 0.001234

colnames(final_file)

library(dplyr)
library(ggplot2)
final_file1 <- final_file %>%
  rename(
    gnomad_v4_AF = gnomAD_AF.y,
    gnomad_v3_AF = gnomad_3_1_1_AF
  )
str(final_file1[, c("gnomad_v3_AF", "gnomad_v4_AF")])
final_file1 <- final_file1 %>%
  mutate(
    gnomad_v3_AF = as.numeric(gnomad_v3_AF),
    gnomad_v4_AF = as.numeric(gnomad_v4_AF)
  )
final_file1 <- final_file1 %>%
  mutate(
    gnomad_v3_AF = ifelse(is.na(gnomad_v3_AF), 0, gnomad_v3_AF),
    gnomad_v4_AF = ifelse(is.na(gnomad_v4_AF), 0, gnomad_v4_AF)
  )
# Scatter plot
ggplot(df, aes(x = gnomad_v3_AF, y = gnomad_v4_AF)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  labs(
    title = "Correlation of gnomAD v3 vs v4 Allele Frequencies",
    x = "gnomAD v3 Allele Frequency",
    y = "gnomAD v4 Allele Frequency"
  ) +
  theme_minimal()
cor(df$gnomad_v3_AF, df$gnomad_v4_AF, use = "complete.obs", method = "pearson")

write.csv(final_file1, "silent_variants_df_with_gnomad_v4_2.csv", 
          row.names = FALSE)

#--------------using saved annotated files--------
final_file1 <- read_csv("silent_variants_df_with_gnomad_v4_2.csv")

annotated_eg <- read_csv("annotated_only.csv")

ggplot(final_file1, aes(x = log(gnomad_v3_AF), y = log(gnomad_v4_AF))) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  labs(
    title = "Correlation of gnomAD v3 vs v4 Allele Frequencies",
    x = "gnomAD v3 Allele Frequency",
    y = "gnomAD v4 Allele Frequency") 
  
cor(final_file1$gnomad_v3_AF, final_file1$gnomad_v4_AF, use = "complete.obs", method = "pearson")

final_file1 %>%
  distinct(Strand) %>%
  pull(Strand)
file %>%
  distinct(Strand) %>%
  pull(Strand)

#-----debugging NULL values----
colnames(file)
colnames(annotated_eg)
colnames(final_file1)

gnomad_v4_zero_v3_nonzero <- final_file1 %>%
  filter(gnomad_v4_AF == 0 & !is.na(gnomad_v3_AF) & gnomad_v3_AF != 0) %>%
  select(Chromosome, Start_Position, End_Position,
         Reference_Allele, Tumor_Seq_Allele2,
         Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode,
         Sequencer, gnomad_v4_AF, gnomad_v3_AF, gnomAD_AC, gnomAD_AF.x, gnomAD_AN)

gnomad_v4_zero_v3_nonzero_distinct <- final_file1 %>%
  filter(gnomad_v4_AF == 0 & !is.na(gnomad_v3_AF) & gnomad_v3_AF != 0) %>%
  select(Chromosome, Start_Position, End_Position,
         Reference_Allele, Tumor_Seq_Allele2,
         Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode,
         Sequencer, gnomad_v4_AF, gnomad_v3_AF, gnomAD_AC, gnomAD_AF.x, gnomAD_AN) %>%
  distinct(Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, .keep_all = TRUE)

gnomad_v4_zero_v3_nonzero_distinct_top10 <- final_file1 %>%
  filter(gnomad_v4_AF == 0 & !is.na(gnomad_v3_AF) & gnomad_v3_AF != 0) %>%
  select(Chromosome, Start_Position,
         Reference_Allele, Tumor_Seq_Allele2) %>%
  distinct(Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2, .keep_all = TRUE)


test1 <- gnomad_v4_zero_v3_nonzero_distinct_top10 %>%
  mutate(
    CHROM = Chromosome,
    POS = Start_Position,
    REF = Reference_Allele,
    ALT = Tumor_Seq_Allele2
  )

test2 <- gnomad_v4_zero_v3_nonzero_distinct_top10 %>%
  mutate(
    CHROM = Chromosome,
    POS = Start_Position,
    REF = Reference_Allele,
    ALT = Tumor_Seq_Allele2
  )

test3 <- gnomad_v4_zero_v3_nonzero_distinct_top10 %>%
  mutate(
    CHROM = Chromosome,
    POS = Start_Position,
    REF = Reference_Allele,
    ALT = Tumor_Seq_Allele2
  )

annotated <- annotate_variants_batch(test1)

annotated2 <- annotate_variants_batch(test2)

annotated3 <- annotate_variants_batch(test3)
annotated3 <- annotated3 %>%
  mutate(
    Chromosome = CHROM,
    Start_Position = POS,
    Reference_Allele = REF,
    Tumor_Seq_Allele2 = ALT
  )

final_file1 %>% select(ends_with(".y")) %>% colnames()
library(dplyr)
library(stringr)

final_file1 <- final_file1 %>%
  rename_with(~ str_replace(., "\\.y$", "_1"), ends_with(".y"))
colnames(final_file1)
final_file1 <- final_file1 %>%
  rename(
    gnomAD_4_AF = gnomad_v4_AF,
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
    gnomAD_4_AFR_AF = gnomAD_AFR_AF_1,
    gnomAD_4_AMR_AF = gnomAD_AMR_AF_1,
    gnomAD_4_ASJ_AF = gnomAD_ASJ_AF_1,
    gnomAD_4_EAS_AF = gnomAD_EAS_AF_1,
    gnomAD_4_FIN_AF = gnomAD_FIN_AF_1,
    gnomAD_4_NFE_AF = gnomAD_NFE_AF_1,
    gnomAD_4_SAS_AF = gnomAD_SAS_AF_1,
    gnomAD_4_MID_AF = gnomAD_MID_AF,
    gnomAD_4_AMI_AF = gnomAD_AMI_AF,
    gnomAD_4_REMAINING_AF = gnomAD_REMAINING_AF
  )

colnames(annotated3)
annotated3 <- annotated3 %>%
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
colnames(annotated3)

# Merge two dataframes with coalescing specific column

# Define key columns for matching
key_cols <- c("Chromosome", "Start_Position", "Reference_Allele", "Tumor_Seq_Allele2")

# Define columns to coalesce
coalesce_cols <- c(
  "gnomAD_4_AF", "gnomAD_4_AC", "gnomAD_4_AN", "gnomAD_4_nhomalt",
  "gnomAD_4_nhemialt", "gnomAD_4_genome_AF", "gnomAD_4_genome_AC", "gnomAD_4_genome_AN",
  "gnomAD_4_genome_nhomalt", "gnomAD_4_exome_AF", "gnomAD_4_exome_AC", "gnomAD_4_exome_AN",
  "gnomAD_4_exome_nhomalt", "gnomAD_4_max_AF", "gnomAD_4_AF_popmax", "gnomAD_4_popmax_population",
  "gnomAD_4_flags", "gnomAD_4_in_gnomAD", "gnomAD_4_AFR_AF", "gnomAD_4_AMR_AF",
  "gnomAD_4_ASJ_AF", "gnomAD_4_EAS_AF", "gnomAD_4_FIN_AF", "gnomAD_4_NFE_AF",
  "gnomAD_4_SAS_AF", "gnomAD_4_MID_AF", "gnomAD_4_AMI_AF", "gnomAD_4_REMAINING_AF"
)

# Perform left join on the key columns
merged_df <- final_file1 %>%
  left_join(annotated3, by = key_cols, suffix = c("", ".from_df1"))

# Coalesce the specified columns
# Replace 0 values in df2 columns with values from df1
for (col in coalesce_cols) {
  col_from_df1 <- paste0(col, ".from_df1")
  
  # Check if both columns exist in merged dataframe
  if (col %in% names(merged_df) && col_from_df1 %in% names(merged_df)) {
    # Replace 0 values with values from df1, keep non-zero values from df2
    merged_df[[col]] <- ifelse(
      merged_df[[col]] == 0 | is.na(merged_df[[col]]),
      merged_df[[col_from_df1]],
      merged_df[[col]]
    )

  
    # Filter for rows where gnomAD_4_AF is still 0
    zero_af_df <- merged_df %>%
      filter(gnomAD_4_AF == 0 | is.na(gnomAD_4_AF)) %>%
      select(Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2)
    
    # Calculate counts and percentages
    total_rows <- nrow(merged_df)
    zero_af_count <- nrow(zero_af_df)
    zero_af_percent <- (zero_af_count / total_rows) * 100
    
    # Print results
    cat("\n=== gnomAD_4_AF = 0 Analysis ===\n")
    cat(sprintf("Total rows in merged data: %d\n", total_rows))
    cat(sprintf("Rows with gnomAD_4_AF = 0: %d\n", zero_af_count))
    cat(sprintf("Percentage with gnomAD_4_AF = 0: %.2f%%\n", zero_af_percent))
    
    # Display the filtered dataframe
    print(zero_af_df)  
    # Remove the temporary column from df1
    merged_df[[col_from_df1]] <- NULL
  }
}

# Remove any remaining .from_df1 columns that weren't in coalesce_cols
cols_to_remove <- grep("\\.from_df1$", names(merged_df), value = TRUE)
if (length(cols_to_remove) > 0) {
  merged_df <- merged_df %>% select(-all_of(cols_to_remove))
}

# Result is in merged_df
print(paste("Merged dataframe has", nrow(merged_df), "rows"))
print(paste("Original df2 had", nrow(annotated3), "rows"))
