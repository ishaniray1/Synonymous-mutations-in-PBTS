
# DATA LOADING & PREPARATION
histologies <- read.delim("histologies.tsv")
library(readr)
silent_variants_df <- read_csv("silent_variants_df.csv")
mut=silent_variants_df
mut=mut[mut$Tumor_Sample_Barcode %in% histologies$Kids_First_Biospecimen_ID, ]
histologies$Tumor_Sample_Barcode= histologies$Kids_First_Biospecimen_ID
histologies$Kids.First.Participant.ID= histologies$Kids_First_Participant_ID

library(stringr)
library(tidyverse)    # dplyr, ggplot2, tidyr, readr
library(pheatmap)     # heatmaps
library(ggplot2)
library(ggrepel)
library(viridis)
library(ComplexHeatmap)
library(circlize)
library(igraph)
library(ggraph)
library(tidyverse)

# FILTER FOR A→X MUTATIONS
str_detect(mut$HGVSc, "^c\\.[0-9]+[ACGT]>A") 
mut=mut[str_detect(mut$HGVSc, "^c\\.[0-9]+A>[ACGT]"),]

# CALCULATE VARIANT ALLELE FREQUENCIES
# ensure numeric columns are numeric
num_cols <- c("t_alt_count","t_depth","n_alt_count","n_depth","gnomad_3_1_1_AF")
mut <- mut %>% mutate(across(all_of(num_cols), ~ as.numeric(.)))
mut$gnomad_3_1_1_AF=as.numeric(mut$gnomad_3_1_1_AF)
mut$gnomad_3_1_1_AF[is.na(as.numeric(mut$gnomad_3_1_1_AF))]=0
#mut=mut[as.numeric(mut$gnomad_3_1_1_AF)<0.01, ]
#mut=mut[mut$gnomad_3_1_1_FILTER %in% c(".", "PASS", "AS_VQSR"),]
mut$fracN=as.numeric(mut$n_alt_count)/as.numeric(mut$n_depth)
mut$fracT=as.numeric(mut$t_alt_count)/as.numeric(mut$t_depth)
#mut=mut[mut$FILTER=="PASS",]
##mut=mut[mut$fracN==0, ]

# CREATE GENOMIC RANGES (for future m6A overlap
chroms <- sub("^chr", "", mut$Chromosome, ignore.case=TRUE)
library(GenomicRanges)
library(IRanges)
sum(is.na(mut$Start_Position))
sum(is.na(mut$End_Position))
mut[is.na(mut$Start_Position) | is.na(mut$End_Position), ]
mut <- mut[!is.na(mut$Start_Position) & !is.na(mut$End_Position), ]
chroms <- sub("^chr", "", mut$Chromosome, ignore.case=TRUE)
mut_gr <- GRanges(
  seqnames = chroms,
  ranges   = IRanges(
    start = as.integer(mut$Start_Position),
    end   = as.integer(mut$End_Position)
  )
)

# now these will overlap properly against m6A_intervals
#mut$in_m6A <- countOverlaps(mut_gr, atlas_gr) > 0 | countOverlaps(mut_gr, mirep_gr) > 0
#mut <- filter(mut, in_m6A)

#mut <- anti_join(mut, atoI,
                 #by = c("Chromosome"="seqname",
                        #"Start_Position"="pos"))


#mut=mut[mut$Hugo_Symbol %in% names(table(mut$Hugo_Symbol)[table(mut$Hugo_Symbol)>1]), ]
#mut=mut[mut$Hugo_Symbol %in% census$`Gene Symbol` ,]

#Creating unique mutation IDs and merging with histologies
mut2 <- mut %>%
  mutate(
    mut_id = paste(Tumor_Sample_Barcode,
      Hugo_Symbol,
      Chromosome,
      Start_Position,
      End_Position,
      Reference_Allele,
      Tumor_Seq_Allele2,
      sep = "_"
    )
  )

#Merge with histology to get participant IDs and diagnosis
mut2 <- mut2 %>%
  left_join(
    histologies%>% select(Tumor_Sample_Barcode, 
                    Kids.First.Participant.ID,
                    harmonized_diagnosis, cancer_group),
    by = "Tumor_Sample_Barcode"
  ) %>%
  filter(!is.na(Kids.First.Participant.ID))   # remove orphaned samples
#unique patient-gene combination identifier
mut2 <- mut2 %>%
  mutate(
    mut_id2= paste(Kids.First.Participant.ID,
                   Hugo_Symbol,
                   sep = "_"
    )
  )
# Keep first occurrence per patient-gene combination
mut3=mut2[!duplicated(mut2$mut_id2),]

# Before deduplication (mut2)
nrow(mut2)
# After deduplication (mut3)
nrow(mut3) 
# duplicates were removed
nrow(mut2) - nrow(mut3)
# Save as CSVs
write.csv(mut2, "mut2_all_mutations_histologies.csv", row.names = FALSE)
write.csv(mut3, "mut3_unique_patient_gene.csv", row.names = FALSE)

cat("Total A→X mutations:", nrow(mut), "\n")
cat("Unique patient-gene combinations:", nrow(mut3), "\n")
cat("  - Unique patients:", n_distinct(mut3$Kids.First.Participant.ID), "\n")
cat("  - Unique genes:", n_distinct(mut3$Hugo_Symbol), "\n")
cat("  - Unique tumor types:", n_distinct(mut3$harmonized_diagnosis), "\n")

#WER binding
#mut3 <- mut3 %>%
  #rowwise() %>%
  #mutate(
    #matched = list(
      #binding_expanded %>%
        #filter(chr == Chromosome, target_gene==Hugo_Symbol, Start_Position >= range_start, Start_Position <= range_end)
    #),
    #has_binder = ifelse(nrow(matched) > 0, "Yes", "No"),
    #binder_proteins = ifelse(has_binder == "Yes", paste(unique(matched$wer), collapse = ","), NA_character_),
    #3m6a_dependency = ifelse(has_binder == "Yes", paste(unique(matched$mod_dependent), collapse = ","), NA_character_)
  #) %>%
  #select(-matched) %>%
  #ungroup()

#mut3=mut3[!is.na(mut3$binder_proteins),]


# --- Step 1: Unique patient-gene-tumor combinations and add mutation metrics ---
mut_metrics <- mut3 %>%
  group_by(Hugo_Symbol, harmonized_diagnosis, cancer_group) %>%
  summarise(
    n_patients_with_mutation = n_distinct(Kids.First.Participant.ID),  # unique patients
    n_mutations = n(),
    n_distinct_mutations = n_distinct(paste(Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, sep = "_")),  # distinct mutations based on coordinates
    avg_gnomad_3_1_1_AF = mean(as.numeric(gnomad_3_1_1_AF), na.rm = TRUE),  # average popmax AF
    avg_in_tumor_expression = mean(fracT, na.rm = TRUE),  # average tumor expression
    avg_in_normal_expression = mean(fracN, na.rm = TRUE),  # average normal expression
    .groups = "drop"
  )

# --- Step 2: Get total unique patients per tumor type from histologies---
total_patients_per_tumor <- histologies%>%
  distinct(Kids.First.Participant.ID, harmonized_diagnosis) %>%
  group_by(harmonized_diagnosis) %>%
  summarise(
    total_patients = n(),
    .groups = "drop"
  )

# --- Step 3: Merge and calculate mutation frequency ---
mut_summary <- mut_metrics %>%
  left_join(total_patients_per_tumor, by = "harmonized_diagnosis") %>%
  mutate(
    freq_in_disease = n_patients_with_mutation / total_patients,
    percent_patients = freq_in_disease * 100
  ) %>%
  arrange(desc(freq_in_disease)) # optional: arrange by most frequently mutated genes first

# --- Step 4: View final table ---
glimpse(mut_summary)
mut_summary_recurrent=mut_summary[mut_summary$n_patients_with_mutation>1,] #filter for recurrent mutations?
cat("Gene×Tumor pairs with recurrent mutations:", nrow(mut_summary_recurrent), "\n")

write.csv(mut_summary, "mut_summary_all.csv", row.names = FALSE)
write.csv(mut_summary_recurrent, "mut_summary_recurrent.csv", row.names = FALSE)

# ----BY GENE ---
# Most frequently mutated genes (across all tumors)
gene_summary <- mut3 %>%
  group_by(Hugo_Symbol) %>%
  summarise(
    n_mutations = n(),
    n_distinct_mutations = n_distinct(paste(Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, sep = "_")),
    n_patients = n_distinct(Kids.First.Participant.ID),
    n_tumor_types = n_distinct(harmonized_diagnosis),
    median_fracT = median(fracT, na.rm=TRUE),
    median_fracN = median(fracN, na.rm=TRUE),
    mean_gnomad_AF = mean(gnomad_3_1_1_AF, na.rm=TRUE),
    tumor_types = paste(unique(harmonized_diagnosis), collapse = "; ")
  ) %>%
  arrange(desc(n_patients), desc(n_mutations))

print(gene_summary %>% select(Hugo_Symbol, n_patients, n_tumor_types) %>% head(5))
# Top genes 
top_genes <- gene_summary %>% head(30) %>% pull(Hugo_Symbol)

# Visualization Top 20 genes
a <- gene_summary %>%
  head(20) %>%
  ggplot(aes(x = reorder(Hugo_Symbol, n_patients), y = n_patients)) +
  geom_col(aes(fill = n_tumor_types)) +
  geom_text(aes(label = paste0(n_patients, " pts / ", n_distinct_mutations, " muts")),
            hjust = -0.05, size = 3) +
  coord_flip() +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  labs(title = "Top 20 Genes with A→X Mutations",
       subtitle = "Number of patients and distinct mutations per gene",
       x = "Gene", y = "Number of Patients",
       fill = "Tumor Types") +
  theme_minimal()
a

#--- Mutation type (nucleotide) breakdown---
mutation_types <- mut3 %>%
  mutate(mutation_type = paste0("A→", str_extract(Tumor_Seq_Allele2, "[ACGT]$"))) %>%
  count(mutation_type) %>%
  mutate(percent = n/sum(n)*100) %>%
  arrange(desc(n))
print(mutation_types)

mut3 %>%
  count(Reference_Allele, Tumor_Seq_Allele2) %>%
  arrange(desc(n))

# ---BY TUMOR TYPE---
tumor_burden <- mut3 %>%
  group_by(harmonized_diagnosis) %>%
  summarise(
    n_mutations = n(),
    n_patients_with_mut = n_distinct(Kids.First.Participant.ID),
    n_distinct_mutations = n_distinct(paste(Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2, sep = "_")),
    n_genes = n_distinct(Hugo_Symbol),
    median_fracT = median(fracT, na.rm=TRUE),
    median_fracN = median(fracN, na.rm=TRUE),
    mean_gnomad_AF = mean(gnomad_3_1_1_AF, na.rm=TRUE)
  ) %>%
  left_join(total_patients_per_tumor, by = "harmonized_diagnosis") %>%
  mutate(
    mutations_per_patient = round(n_mutations / n_patients_with_mut, 2),
    percent_patients_affected = round((n_patients_with_mut / total_patients) * 100, 1)
  ) %>%
  arrange(desc(mutations_per_patient))

print(tumor_burden %>% head(10), n = 10)

# Top tumors 
top_tumors <- tumor_burden %>%
  head(20) %>%
  pull(harmonized_diagnosis)
# Visualization 
tumor_burden %>%
  head(20) %>% 
  ggplot(aes(x = reorder(harmonized_diagnosis, mutations_per_patient), 
             y = mutations_per_patient)) +
  geom_col(aes(fill = percent_patients_affected)) +
  geom_text(aes(label = paste0(round(mutations_per_patient, 1), " / ", n_distinct_mutations)),
            hjust = -0.2, size = 3) +
  coord_flip() +
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(title = "A→X Mutation Burden by Tumor Type",
       subtitle = "Mutations per affected patient / Unique mutations",
       x = "Tumor Type", 
       y = "Mutations per Patient",
       fill = "% Patients\nAffected") +
  theme_minimal()


# ---VAF ----
summary(mut3$fracT)  # Tumor VAF
summary(mut3$fracN)  # Normal VAF
# tumor VAF
ggplot(mut3, aes(x = fracT)) +
  geom_histogram(binwidth = 0.01) +
  labs(x = "Tumor VAF (fracT)", y = "Count")
# VAF by tumor type
ggplot(mut3, aes(x = reorder(harmonized_diagnosis, fracT, FUN=median), y = fracT)) +
  geom_boxplot() + coord_flip() +
  labs(x = "", y = "Tumor VAF")
# VAF: Tumor vs Normal
p_vaf <- ggplot(mut3, aes(x = fracN, y = fracT)) +
  geom_point(alpha = 0.3, color = "darkred") +
  geom_abline(slope = 1, linetype = "dashed") +
  labs(
    title = "Tumor vs Normal VAF for A→X Mutations",
    x = "Normal VAF",
    y = "Tumor VAF"
  ) +
  theme_minimal()
p_vaf

# gnomAD AF distribution
summary(mut3$gnomad_3_1_1_AF)
cat("Mutations with AF=0 (not in gnomAD):", 
    sum(mut3$gnomad_3_1_1_AF == 0, na.rm=TRUE), 
    "(", round(mean(mut3$gnomad_3_1_1_AF == 0, na.rm=TRUE)*100, 1), "%)\n")


# ---RECURRENCE ---
# Hotspot mutations (same position in multiple patients)
hotspots <- mut3 %>%
  group_by(Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2) %>%
  summarise(
    n_patients = n_distinct(Kids.First.Participant.ID),
    genes = paste(unique(Hugo_Symbol), collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(n_patients > 2) %>%
  arrange(desc(n_patients))
print(hotspots)

#----
#Mutations exclusive to single tumor type
exclusive <- mut3 %>%
  group_by(Hugo_Symbol) %>%
  summarise(
    n_tumor_types = n_distinct(harmonized_diagnosis),
    n_patients = n_distinct(Kids.First.Participant.ID),
    tumor_types = paste(unique(harmonized_diagnosis), collapse = "; ")
  ) %>%
  mutate(
    exclusive = ifelse(n_tumor_types == 1, "Exclusive", 
                         ifelse(n_tumor_types <= 3, "Restricted", "Broad"))
  ) %>%
  arrange(desc(n_patients))
print(exclusive)
e <- exclusive %>% filter(exclusive == "Exclusive")
e

# all exclusive genes per tumor type
exclusive_genes <- mut3 %>%
  group_by(Hugo_Symbol) %>%
  filter(n_distinct(harmonized_diagnosis) == 1) %>%
  ungroup() %>%
  group_by(harmonized_diagnosis) %>%
  summarise(
    n_exclusive_genes = n_distinct(Hugo_Symbol),
    n_exclusive_mutations = n(),
    exclusive_genes = paste(unique(Hugo_Symbol), collapse = "; ")
  )
print(exclusive_genes)

tumor_specificity <- mut3 %>%
  group_by(Hugo_Symbol) %>%
  mutate(total_mutations_in_gene = dplyr::n()) %>%           
  group_by(Hugo_Symbol, harmonized_diagnosis) %>%
  summarise(
    n_in_tumor = n(),                                
    total_mutations = dplyr::first(total_mutations_in_gene), 
    specificity = n_in_tumor / total_mutations,      
    .groups = "drop"
  ) %>%
  arrange(desc(specificity))
highly_specific <- tumor_specificity %>%
  filter(specificity > 0.6, n_in_tumor > 2)

# --- genes  co-mutated in same patients?--
patient_gene_matrix <- mut3 %>%
  distinct(Kids.First.Participant.ID, Hugo_Symbol, harmonized_diagnosis, cancer_group) %>%
  group_by(Kids.First.Participant.ID) %>%
  mutate(n_mutations = n_distinct(Hugo_Symbol)) %>%
  # Keep only patients with 2+ mutations
  filter(n_mutations >= 2) %>%
  summarise(
    harmonized_diagnosis = paste(unique(harmonized_diagnosis), collapse = "; "),
    genes = list(Hugo_Symbol),
    cancer_group = dplyr::first(cancer_group),
    .groups = "drop"
  ) %>%
  unnest(genes) %>%
  # pivot
  mutate(present = 1) %>%
  pivot_wider(
    id_cols = c(Kids.First.Participant.ID, harmonized_diagnosis, cancer_group),
    names_from = genes,
    values_from = present,
    values_fill = 0
  )

patient_gene_matrix %>%
  group_by(Kids.First.Participant.ID) %>%
  filter(n() > 1) %>%
  arrange(Kids.First.Participant.ID)

patient_gene_matrix %>%
  count(Kids.First.Participant.ID, sort = TRUE)

# How many mutations does each patient have?
patient_burden <- patient_gene_matrix2 %>%
  mutate(total_mutations = rowSums(select(., -Kids.First.Participant.ID, -harmonized_diagnosis)))

ggplot(patient_burden, aes(x = total_mutations)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
  theme_minimal() +
  labs(title = "Distribution of Mutation Burden",
       x = "Number of Mutations per Patient",
       y = "Number of Patients")

# diagnoses with at least 10 patients
top_diagnoses <- patient_burden %>%
  count(harmonized_diagnosis) %>%
  filter(n >= 10) %>%
  pull(harmonized_diagnosis)

patient_burden_filtered <- patient_burden %>%
  filter(harmonized_diagnosis %in% top_diagnoses)

p_burden <- ggplot(patient_burden_filtered, 
                   aes(x = reorder(harmonized_diagnosis, total_mutations, median), 
                       y = total_mutations, 
                       fill = harmonized_diagnosis)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.2, size = 1) +
  coord_flip() +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = "Mutation Burden by Diagnosis (≥10 patients)",
    x = "Diagnosis",
    y = "Number of Mutations"
  )
p_burden

patient_burden_group <- patient_gene_matrix %>%
  mutate(total_mutations = rowSums(select(., -Kids.First.Participant.ID, -cancer_group, -harmonized_diagnosis)))
top_groups <- patient_burden %>%
  count(cancer_group) %>%
  filter(n >= 10) %>%
  pull(cancer_group)
patient_burden_filtered_group <- patient_burden_group %>%
  filter(cancer_group %in% top_groups)
p_burden_group <- ggplot(patient_burden_filtered_group, 
                         aes(x = reorder(cancer_group, total_mutations, median), 
                             y = total_mutations, 
                             fill = cancer_group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.2, size = 1) +
  coord_flip() +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = "Mutation Burden by Cancer Group (≥10 patients)",
    x = "Cancer Group",
    y = "Number of Mutations"
  )
p_burden_group

patient_gene_matrix2 <- mut3 %>%
  distinct(Kids.First.Participant.ID, Hugo_Symbol, harmonized_diagnosis, cancer_group) %>%
  group_by(Kids.First.Participant.ID) %>%
  mutate(n_mutations = n_distinct(Hugo_Symbol)) %>%
  filter(n_mutations >= 2) %>%
  summarise(
    harmonized_diagnosis = paste(unique(harmonized_diagnosis), collapse = "; "),
    genes = list(Hugo_Symbol),
    .groups = "drop"
  ) %>%
  unnest(genes) %>%
  mutate(present = 1) %>%
  pivot_wider(
    id_cols = c(Kids.First.Participant.ID, harmonized_diagnosis),
    names_from = genes,
    values_from = present,
    values_fill = 0
  )
#  pairwise co-occurrence
gene_cols <- setdiff(names(patient_gene_matrix2), c("Kids.First.Participant.ID", "harmonized_diagnosis"))
# co-occurrence matrix 
cooccur_df <- data.frame()
for(i in 1:(length(gene_cols)-1)) {
  for(j in (i+1):length(gene_cols)) {
    gene1 <- gene_cols[i]
    gene2 <- gene_cols[j]
    co_count <- sum(patient_gene_matrix2[[gene1]] == 1 & patient_gene_matrix[[gene2]] == 1)
    if(co_count > 0) {
      cooccur_df <- rbind(cooccur_df, data.frame(
        gene1 = gene1,
        gene2 = gene2,
        co_occurrence = co_count,
        stringsAsFactors = FALSE
      ))
    }
  }
}
#  20 most frequent co-mutations
top_pairs <- cooccur_df %>%
  arrange(desc(co_occurrence)) %>%
  head(20) %>%
  mutate(pair = paste(gene1, gene2, sep = " - "))

p_cooccurrence <- ggplot(top_pairs, aes(x = reorder(pair, co_occurrence), y = co_occurrence)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = co_occurrence), hjust = -0.2, size = 3) +
  coord_flip() +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = "Top 20 Co-mutated Gene Pairs",
    x = "Gene Pair",
    y = "Number of Patients with Both Mutations"
  )
p_cooccurrence

# full diagnosis for each pair
pair_diagnosis_full <- data.frame()

for(i in 1:nrow(top_pairs)) {
  gene1 <- top_pairs$gene1[i]
  gene2 <- top_pairs$gene2[i]
  
  patients_with_both <- patient_gene_matrix2 %>%
    filter(.data[[gene1]] == 1 & .data[[gene2]] == 1) %>%
    count(harmonized_diagnosis, name = "count")
  
  patients_with_both$pair <- paste(gene1, gene2, sep = " - ")
  pair_diagnosis_full <- rbind(pair_diagnosis_full, patients_with_both)
}

#  stacked
p_cooccurrence_diag <- ggplot(pair_diagnosis_full, 
                              aes(x = reorder(pair, count, sum), 
                                  y = count, 
                                  fill = harmonized_diagnosis)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = "Top 20 Co-mutated Gene Pairs by Diagnosis",
    x = "Gene Pair",
    y = "Number of Patients",
    fill = "Diagnosis"
  )
p_cooccurrence_diag

#  network
cooccur_filtered <- cooccur_df %>%
  filter(co_occurrence >= 3) 
#  network
g <- graph_from_data_frame(cooccur_filtered, directed = FALSE)
# Plot
p_network <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(width = co_occurrence, alpha = co_occurrence), color = "gray50") +
  geom_node_point(aes(size = degree(g)), color = "steelblue", alpha = 0.8) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3, max.overlaps = 20) +
  scale_edge_width(range = c(0.5, 3)) +
  scale_size(range = c(3, 10)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title = "Gene Co-mutation Network",
    subtitle = "Node size = number of connections, Edge width = co-occurrence count"
  )
p_network

#--Visualization

p_tile <- mut_summary %>%
  filter(
    Hugo_Symbol %in% top_genes,
    harmonized_diagnosis %in% top_tumors,
    n_patients_with_mutation > 1
  ) %>%
  ggplot(aes(x = harmonized_diagnosis, y = Hugo_Symbol)) +
  geom_tile(aes(fill = percent_patients), color = "white") +
  geom_text(aes(label = n_patients_with_mutation), size = 2.5, color = "black") +
  scale_fill_gradient(low = "white", high = "darkred", 
                      name = "% Patients") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = "A→X Mutation Frequency",
    subtitle = "Numbers indicate patient counts",
    x = "Tumor Type",
    y = "Gene"
  )
p_tile

#by harmonized diagnosis
freq_mat <- mut_summary %>%
  select(Hugo_Symbol, harmonized_diagnosis, n_patients_with_mutation) %>%
  pivot_wider(
    names_from = harmonized_diagnosis,
    values_from = n_patients_with_mutation,
    values_fill = 0
  )
mat <- as.matrix(freq_mat[, -1])
rownames(mat) <- freq_mat$Hugo_Symbol
top_row_idx <- order(apply(mat, 1, sum), decreasing = TRUE)[1:50]
pdf("heatmap_top50_genes_by_tumor.pdf", width = 16, height = 12)  # bigger
pheatmap(
  mat[top_row_idx, ],
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 8,   # slightly larger if enough space
  fontsize_col = 7,
  color = viridis(50),
  main = "Top 50 genes: A→X Mutations by Patient Count",
  display_numbers = FALSE,
  border_color = "grey90"
)
dev.off()

#by cancer group
freq_group_mat <- mut_summary %>%
  group_by(Hugo_Symbol, cancer_group) %>%
  summarise(
    n_patients_with_mutation = sum(n_patients_with_mutation, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = cancer_group,
    values_from = n_patients_with_mutation,
    values_fill = 0
  )
mat_group <- as.matrix(freq_group_mat[, -1])
rownames(mat_group) <- freq_group_mat$Hugo_Symbol
top_row_idx <- order(apply(mat_group, 1, sum), decreasing = TRUE)[1:50]
pdf("heatmap_top50_genes_by_cancer_group.pdf", width = 12, height = 10)  
pheatmap(
  mat_group[top_row_idx, ],
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 6,
  fontsize_col = 8,
  color = viridis(50),
  main = "A→X Mutations by Cancer Group (Top 50 Genes)",
  display_numbers = FALSE,
  border_color = "grey90"
)
dev.off()




