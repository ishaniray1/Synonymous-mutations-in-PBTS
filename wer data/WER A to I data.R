#================
#RM2Target http://rm2target.canceromics.org/#/download

#Validated---------------------------------------------
WER_validated_data <- read.csv("Validated_Homo_sapiens_m6A.csv", stringsAsFactors = FALSE)
head(WER_validated_data )
str(WER_validated_data )
colnames(WER_validated_data )
WER_validated_subset <- WER_validated_data[
  , c("rm2target_id",
      "wer",
      "wer_type",
      "target_gene",
      "cell_line",
      "pmid",
      "target_site",
      "target_region")
]
WER_validated_subset$wer_evidence <- "validated"
colnames(WER_validated_subset)
head(WER_validated_subset)


#Binding-----------------------------------------
path <- "C:/Users/isr4005/Downloads/Binding_Homo_sapiens_m6A (1)/"
binding_clip  <- read.csv(paste0(path, "Binding_CLIP_hg38_Homo_sapiens_m6A.csv"))
binding_ms    <- read.csv(paste0(path, "Binding_MS_Homo_sapiens_m6A.csv"))
binding_rip   <- read.csv(paste0(path, "Binding_RIP_Homo_sapiens_m6A.csv"))
binding_chip  <- read.csv(paste0(path, "Binding_CHIP_hg38_Homo_sapiens_m6A.csv"))
colnames(binding_clip)
colnames(binding_ms)
colnames(binding_rip)
colnames(binding_chip)

library(tidyr)
library(dplyr)
binding_chip_clean <- binding_chip %>%
  select(
    rm2target_id,
    wer,
    target_gene,
    cell_line,
    pmid,
    gse_id,
    wer_type,
    chr,
    strand,
    position
  ) %>%
  # evidence 
  mutate(evidence = "chip_binding") 
library(stringr)
binding_chip_clean <- binding_chip_clean %>%
  janitor::clean_names() %>%  
  mutate(position = as.character(position)) %>%
  # split row into multiple rows based on commas
  separate_rows(position, sep = ",") %>%
  # extract only numeric ranges
  mutate(position_clean = str_extract(position, "[0-9]+-[0-9]+")) %>%
  # split into start and end columns
  separate(position_clean, into = c("start", "end"), sep = "-", convert = TRUE) %>%
  # select 
  select(
    rm2target_id,
    wer,
    target_gene,
    cell_line,
    pmid,
    gse_id,
    wer_type,
    chr,
    strand,
    start,
    end,
    evidence
  )

binding_clip_clean <- binding_clip %>%
  select(
    rm2target_id,
    wer,
    target_gene,
    cell_line,
    pmid,
    gse_id,
    wer_type,
    chr,
    strand,
    position
  ) %>%
  # evidence 
  mutate(evidence = "clip_binding") 
binding_clip_clean <- binding_clip_clean %>%
  janitor::clean_names() %>%  
  mutate(position = as.character(position)) %>%
  # split row into multiple rows based on commas
  separate_rows(position, sep = ",") %>%
  # extract only numeric ranges
  mutate(position_clean = str_extract(position, "[0-9]+-[0-9]+")) %>%
  # split into start and end columns
  separate(position_clean, into = c("start", "end"), sep = "-", convert = TRUE) %>%
  # select 
  select(
    rm2target_id,
    wer,
    target_gene,
    cell_line,
    pmid,
    gse_id,
    wer_type,
    chr,
    strand,
    start,
    end,
    evidence
  )


binding_rip_clean <- binding_rip %>%
  janitor::clean_names() %>% 
  select(
    rm2target_id,
    wer,
    target_gene,
    cell_line,
    pmid,
    gse_id,
    log2fc,
    wer_type
  ) %>%
  mutate(evidence = "rip_binding")

binding_ms_clean <- binding_ms %>%
  janitor::clean_names() %>% 
  select(
    rm2target_id,
    wer,
    target_gene,
    cell_line,
    pmid,
    uniprot_id,
    uniprot_name,
    wer_type
  ) %>%
  mutate(evidence = "ms_binding")

#Perturbation--------------------------------------
path_p <- "C:/Users/isr4005/Downloads/Perturbation_Homo_sapiens_m6A (1)/"

perturb_modrip <- read.csv(paste0(path_p, "Perturbation_Mod-RIP_hg38_Homo_sapiens_m6A.csv"))
perturb_stability   <- read.csv(paste0(path_p, "Perturbation_Stability_Homo_sapiens_m6A.csv"))
perturb_te     <- read.csv(paste0(path_p, "Perturbation_TE_Homo_sapiens_m6A.csv"))
perturb_expr   <- read.csv(paste0(path_p, "Perturbation_Expr_Homo_sapiens_m6A.csv"))
perturb_as     <- read.csv(paste0(path_p, "Perturbation_AS_Homo_sapiens_m6A.csv"))
colnames(perturb_modrip)
colnames(perturb_stability)
colnames(perturb_te)
colnames(perturb_expr)
colnames(perturb_as)

perturb_modrip_clean <- perturb_modrip %>%
  clean_names() %>% 
  select(
    rm2target_id,
    wer,
    target_gene,
    cell_line,
    pmid,
    gse_id,
    peak,
    strand,
    log2fc,
    padj,
    direction,
    perturbation,
    wer_type
  ) %>%
  # split peak col into chr, start, end
  mutate(
    peak = str_replace_all(peak, "chr", "chr")  # just in case
  ) %>%
  separate(peak, into = c("chr", "range"), sep = ":", remove = TRUE) %>%
  separate(range, into = c("start", "end"), sep = "-", convert = TRUE) %>%
  # evidence column
  mutate(evidence = "perturbation_modrip")

perturb_stability_clean <- perturb_stability %>%
  clean_names() %>% 
  select(
    rm2target_id,
    wer,
    target_gene,
    cell_line,
    pmid,
    gse_id,
    log2fc,
    direction,
    perturbation,
    wer_type
  ) %>%
  mutate(evidence = "perturbation_stability")

library(dplyr)
library(janitor)

perturb_te_clean <- perturb_te %>%
  clean_names() %>% 
  select(
    rm2target_id,
    wer,
    target_gene,
    cell_line,
    pmid,
    gse_id,
    te,
    rpf_log2fc,
    direction,
    perturbation,
    wer_type
  ) %>%
  mutate(evidence = "perturbation_te")

library(dplyr)
library(janitor)

perturb_expr_clean <- perturb_expr %>%
  clean_names() %>%  
  select(
    rm2target_id,
    wer,
    target_gene,
    cell_line,
    pmid,
    gse_id,
    log2fc,
    padj,
    direction,
    perturbation,
    wer_type
  ) %>%
  mutate(evidence = "perturbation_expr")

perturb_as_clean <- perturb_as %>%
  clean_names() %>% 
  select(
    rm2target_id,
    wer,
    target_gene,
    cell_line,
    pmid,
    gse_id,
    p_value,
    perturbation,
    wer_type,
    inc_level_difference
  ) %>%
  mutate(evidence = "perturbation_as")

#All-------------------------------
library(dplyr)

#ensure all desired columns exist
ensure_cols <- function(df, cols) {
  missing_cols <- setdiff(cols, colnames(df))
  if(length(missing_cols) > 0){
    df[missing_cols] <- NA
  }
  df[, cols]
}

# unified columns
common_cols <- c(
  "wer",
  "wer_type",
  "target_gene",
  "chr",
  "strand",
  "start",
  "end",
  "uniprot_name",
  "log2fc",
  "inc_level_difference",
  "direction",
  "perturbation",
  "evidence"
)

# validation and binding datasets
binding_chip_clean2 <- ensure_cols(binding_chip_clean, common_cols)
binding_clip_clean2 <- ensure_cols(binding_clip_clean, common_cols)
binding_rip_clean2  <- ensure_cols(binding_rip_clean,  common_cols)
binding_ms_clean2   <- ensure_cols(binding_ms_clean,   common_cols)
WER_validated_subset2 <- ensure_cols(WER_validated_subset, common_cols)

# perturbation datasets
perturb_modrip_clean2 <- perturb_modrip_clean %>%
  ensure_cols(common_cols)

perturb_stability_clean2 <- perturb_stability_clean %>%
  ensure_cols(common_cols)

perturb_expr_clean2 <- perturb_expr_clean %>%
  ensure_cols(common_cols)

perturb_te_clean2 <- perturb_te_clean %>%
  dplyr::rename(log2fc = rpf_log2fc) %>%
  ensure_cols(common_cols)
colnames(perturb_te_clean)

perturb_as_clean2 <- perturb_as_clean %>%
  ensure_cols(common_cols)


# all binding datasets
binding_all <- bind_rows(
  binding_chip_clean2,
  binding_clip_clean2,
  binding_rip_clean2,
  binding_ms_clean2,
  WER_validated_subset2
)

#  all perturbation datasets
perturb_all <- bind_rows(
  perturb_modrip_clean2,
  perturb_stability_clean2,
  perturb_expr_clean2,
  perturb_te_clean2,
  perturb_as_clean2
)

#  everything
all_rm2target <- bind_rows(
  binding_all,
  perturb_all
)
head(all_rm2target)
colnames(all_rm2target)
library(dplyr)

# WERs that appear in multiple wer_types
wer_multiple_types <- all_rm2target %>%
  group_by(wer) %>%
  summarise(
    types = n_distinct(wer_type),
    wer_types_list = paste(unique(wer_type), collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(types > 1) %>%
  arrange(desc(types))

wer_multiple_types #Null (good)

# Distinct target genes per WER
distinct_genes_per_wer <- all_rm2target %>%
  group_by(wer) %>%
  summarise(
    distinct_genes = list(unique(target_gene)),
    num_genes = n_distinct(target_gene),
    .groups = "drop"
  )

distinct_genes_per_wer


#Summary--------------------------------------
library(dplyr)
# total rows and unique target_genes per wer
summary_wer <- all_rm2target %>%
  group_by(wer, wer_type) %>%
  summarise(
    total_rows = n(),
    unique_targets = n_distinct(target_gene),
    .groups = "drop"
  )
summary_wer

#distinct WERs per wer_type
distinct_counts <- all_rm2target %>%
  group_by(wer_type) %>%
  summarise(
    distinct_wer = n_distinct(wer),
    .groups = "drop"
  )
distinct_counts

# Writers
writers_table <- all_rm2target %>%
  filter(wer_type == "Writer") %>%
  group_by(wer) %>%
  summarise(
    total_rows = n(),
    unique_targets = n_distinct(target_gene),
    .groups = "drop"
  ) %>%
  arrange(desc(total_rows))

# Erasers
erasers_table <- all_rm2target %>%
  filter(wer_type == "Eraser") %>%
  group_by(wer) %>%
  summarise(
    total_rows = n(),
    unique_targets = n_distinct(target_gene),
    .groups = "drop"
  ) %>%
  arrange(desc(total_rows))

# Readers
readers_table <- all_rm2target %>%
  filter(wer_type == "Reader") %>%
  group_by(wer) %>%
  summarise(
    total_rows = n(),
    unique_targets = n_distinct(target_gene),
    .groups = "drop"
  ) %>%
  arrange(desc(total_rows))

writers_table
erasers_table
readers_table

# distinct WERs per target gene for each wer_type
wer_per_gene <- all_rm2target %>%
  group_by(target_gene, wer_type) %>%
  summarise(
    distinct_wer_count = n_distinct(wer),
    .groups = "drop"
  ) %>%
  arrange(desc(distinct_wer_count))
wer_per_gene

detailed_wer_per_gene <- all_rm2target %>%
  group_by(target_gene, wer_type) %>%
  summarise(
    distinct_wer_count = n_distinct(wer),
    distinct_wers = paste(sort(unique(wer)), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(distinct_wer_count))
detailed_wer_per_gene

# file paths
rds_path <- "C:/Users/isr4005/Downloads/all_rm2target.rds"
csv_path <- "C:/Users/isr4005/Downloads/all_rm2target.csv"
saveRDS(all_rm2target, file = rds_path)
write.csv(all_rm2target, file = csv_path, row.names = FALSE)

#long form --------------------------------

# Read RDS file
all_rm2target <- readRDS("C:/Users/isr4005/Downloads/all_rm2target.rds")
head(all_rm2target)
str(all_rm2target)

# Filter for writers (case-insensitive) 
witers_subset <- subset(all_rm2target, tolower(wer_type) == "writer")
write.csv(witers_subset, "C:/Users/isr4005/Downloads/witers_subset.csv", row.names = FALSE)
# Filter for readers (case-insensitive)
readers_subset <- subset(all_rm2target, tolower(wer_type) == "reader")
write.csv(readers_subset, "C:/Users/isr4005/Downloads/readers_subset.csv", row.names = FALSE)
# Filter for erasers (case-insensitive)
readers_subset <- subset(all_rm2target, tolower(wer_type) == "eraser")
write.csv(readers_subset, "C:/Users/isr4005/Downloads/erasers_subset.csv", row.names = FALSE)

# identify gene-only entries (no site info)
gene_only_entries <- all_rm2target %>%
  filter(is.na(chr) | is.na(start) | is.na(end)) %>%
  distinct(wer, wer_type, target_gene, .keep_all = TRUE) %>%
  select(wer, wer_type, target_gene, perturbation, direction, log2fc, 
         inc_level_difference, uniprot_name, evidence)

#  distinct gene-site combinations with their writers
gene_site_writers <- all_rm2target %>%
  filter(!is.na(chr) & !is.na(start) & !is.na(end)) %>%
  distinct(target_gene, chr, strand, start, end, wer, wer_type, .keep_all = TRUE)

# For each gene-site combination, add all the gene-only writers
expanded_gene_site <- gene_site_writers %>%
  # Get all distinct gene-site combinations
  distinct(target_gene, chr, strand, start, end) %>%
  # Cross join with gene-only writers for the same gene
  left_join(gene_only_entries, by = "target_gene", relationship = "many-to-many") %>%
  # Combine with the original site-specific writers
  bind_rows(
    gene_site_writers %>%
      select(target_gene, chr, strand, start, end, wer, wer_type, perturbation, 
             direction, log2fc, inc_level_difference, uniprot_name, evidence)
  ) %>%
  # DON'T remove duplicates - add evidence type columns
  mutate(
    validation = ifelse(grepl("validation", evidence, ignore.case = TRUE), TRUE, FALSE),
    binding = ifelse(grepl("binding", evidence, ignore.case = TRUE), TRUE, FALSE),
    perturbation_evidence = ifelse(grepl("perturbation", evidence, ignore.case = TRUE), TRUE, FALSE)
  ) %>%
  # Add source column to track where each entry came from
  group_by(target_gene, chr, strand, start, end, wer, wer_type) %>%
  mutate(
    source = case_when(
      n() > 1 & !is.na(first(na.omit(evidence))) ~ "both",
      !is.na(first(na.omit(evidence))) ~ "site_specific",
      TRUE ~ "gene_general"
    ),
    # Combine evidence from multiple entries
    evidence_combined = ifelse(n() > 1,
                               paste(na.omit(unique(evidence)), collapse = ";"),
                               first(evidence))
  ) %>%
  ungroup()

#  include pure gene-only rows (for genes that only have general)
pure_gene_only <- gene_only_entries %>%
  anti_join(gene_site_writers, by = "target_gene") %>%
  mutate(
    chr = NA_character_,
    strand = NA_character_,
    start = NA_integer_,
    end = NA_integer_,
    validation = ifelse(grepl("validation", evidence, ignore.case = TRUE), TRUE, FALSE),
    binding = ifelse(grepl("binding", evidence, ignore.case = TRUE), TRUE, FALSE),
    perturbation_evidence = ifelse(grepl("perturbation", evidence, ignore.case = TRUE), TRUE, FALSE),
    source = "gene_only",
    evidence_combined = evidence
  )

# Combine everything
final_long <- bind_rows(expanded_gene_site, pure_gene_only) %>%
  arrange(target_gene, chr, start, end, wer, wer_type) %>%
  # For the final output with numbered columns
  group_by(target_gene, chr, strand, start, end) %>%
  mutate(writer_number = row_number()) %>%
  ungroup() %>%
  
  # Create the numbered columns
  select(target_gene, chr, strand, start, end, writer_number, wer, wer_type, 
         perturbation, direction, log2fc, inc_level_difference, uniprot_name, 
         evidence_combined, validation, binding, perturbation_evidence, source) %>%
  
  # Pivot to wider format
  pivot_wider(
    id_cols = c(target_gene, chr, strand, start, end),
    names_from = writer_number,
    values_from = c(wer, wer_type, perturbation, direction, log2fc, 
                    inc_level_difference, uniprot_name, evidence_combined, 
                    validation, binding, perturbation_evidence, source),
    names_glue = "{.value}_{writer_number}"
  )

# Alternative: Keep in long format but with all evidence columns
final_long_alt <- bind_rows(expanded_gene_site, pure_gene_only) %>%
  arrange(target_gene, chr, start, end, wer, wer_type) %>%
  group_by(target_gene, chr, strand, start, end, wer, wer_type) %>%
  summarize(
    perturbation = first(na.omit(perturbation)),
    direction = first(na.omit(direction)),
    log2fc = first(na.omit(log2fc)),
    inc_level_difference = first(na.omit(inc_level_difference)),
    uniprot_name = first(na.omit(uniprot_name)),
    evidence = paste(na.omit(unique(evidence)), collapse = ";"),
    validation = any(validation, na.rm = TRUE),
    binding = any(binding, na.rm = TRUE),
    perturbation_evidence = any(perturbation_evidence, na.rm = TRUE),
    source = paste(unique(source), collapse = ";"),
    .groups = 'drop'
  ) %>%
  # Add writer number for ordering
  group_by(target_gene, chr, strand, start, end) %>%
  mutate(writer_number = row_number()) %>%
  ungroup()

# View the results
head(final_long)  # Wide format with numbered columns
head(final_long_alt)  # Long format with evidence type columns

#==============================
#REDIportal http://srv00.recas.ba.infn.it/atlas/ 



