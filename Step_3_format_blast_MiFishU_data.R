# 02.07.25
# In this script I will filter my BLAST results

library(writexl)
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(phyloseq)

############################# Function #########################################

process_blast_summary <- function(file_path) {
  # Get the file stem (e.g., "RUN1_MiFishU" from "BLAST_summaries/RUN1_MiFishU_BLAST_summary.txt")
  file_stem <- tools::file_path_sans_ext(basename(file_path)) %>%
    str_remove("_BLAST_summary")
  
  # Read and process
  blast_data <- read.delim(file_path, header = FALSE, sep = "\t") %>%
    select(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V16, V17, V18, V20) %>%
    mutate(V1 = str_remove(V1, "_concatenated$")) %>%
    dplyr::rename(
      sample_rep_primer = V1,
      pident = V2,
      length = V3,
      mismatchs = V4,
      gapopens = V5,
      qstart = V6,
      qend = V7,
      sstart = V8,
      send = V9,
      evalue = V10,
      bitscore = V11,
      gen = V16,
      sp = V17,
      lib_ref = V18,
      seq = V20
    ) %>%
    mutate(
      species = str_c(gen, sp, sep = "_"),
      sample = str_split(sample_rep_primer, "_", simplify = TRUE)[, 5] %>%
        paste0("_", str_split(sample_rep_primer, "_", simplify = TRUE)[, 6])
    ) %>%
    select(-gen, -sp)
  
  # Filter by similarity and length
  blast_data <- blast_data %>%
    filter(pident > 99 & length > 150) %>%
    #    filter(pident > 98 & length > 150) %>%
    mutate(species = str_c(species, "_MiFishU"))
  
  # OTU table
  otu_tab <- blast_data %>%
    select(species, sample) %>%
    distinct() %>%
    mutate(present = 1) %>%
    pivot_wider(names_from = sample,
                values_from = present,
                values_fill = 0) %>%
    tibble::column_to_rownames("species")
  
  # Tax table
  tax_tab <- blast_data %>%
    select(species) %>%
    distinct()
  
  rownames(tax_tab) <- tax_tab$species
  
  tax_tab <- tax_tab %>%
    mutate(species = str_remove(species, "_MiFishU$"),
           species = str_replace_all(species, "_", " "))
  
  # Return a named list with dynamic naming
  return(setNames(
    list(
      blast_data = blast_data,
      otu_tab = otu_tab,
      tax_tab = tax_tab
    ),
    nm = paste0(file_stem, c("_data", "_otu_tab", "_tax_tab"))
  ))
}


RUN1_Mi_result <- process_blast_summary("BLAST_summaries/RUN1_MiFishU_BLAST_summary.txt")
RUN2_Mi_result <- process_blast_summary("BLAST_summaries/RUN2_MiFishU_BLAST_summary.txt")
RUN3_Mi_result <- process_blast_summary("BLAST_summaries/RUN3_MiFishU_BLAST_summary.txt")
RUN4_Mi_result <- process_blast_summary("BLAST_summaries/RUN4_MiFishU_BLAST_summary.txt")

############################### MiFishU otu table ##############################

RUN1_Mi_result$RUN1_MiFishU_otu_tab <- RUN1_Mi_result$RUN1_MiFishU_otu_tab %>% mutate(species = rownames(RUN1_Mi_result$RUN1_MiFishU_otu_tab))
RUN2_Mi_result$RUN2_MiFishU_otu_tab <- RUN2_Mi_result$RUN2_MiFishU_otu_tab %>% mutate(species = rownames(RUN2_Mi_result$RUN2_MiFishU_otu_tab))
RUN3_Mi_result$RUN3_MiFishU_otu_tab <- RUN3_Mi_result$RUN3_MiFishU_otu_tab %>% mutate(species = rownames(RUN3_Mi_result$RUN3_MiFishU_otu_tab))
RUN4_Mi_result$RUN4_MiFishU_otu_tab <- RUN4_Mi_result$RUN4_MiFishU_otu_tab %>% mutate(species = rownames(RUN4_Mi_result$RUN4_MiFishU_otu_tab))

# Merge them together by species
Mi_otu_tab <- RUN1_Mi_result$RUN1_MiFishU_otu_tab %>%
  full_join(RUN2_Mi_result$RUN2_MiFishU_otu_tab, by = "species") %>%
  full_join(RUN3_Mi_result$RUN3_MiFishU_otu_tab, by = "species") %>%
  full_join(RUN4_Mi_result$RUN4_MiFishU_otu_tab, by = "species") %>%
  replace(is.na(.), 0) %>%  # Replace NAs with 0 (species absent in that run)
  tibble::column_to_rownames("species")

Mi_otu_tab

Mi_otu_tab <- as.matrix(Mi_otu_tab)

############################### MiFishU tax table ##############################

RUN1_Mi_result$RUN1_MiFishU_tax_tab <- RUN1_Mi_result$RUN1_MiFishU_tax_tab %>% mutate(species = rownames(RUN1_Mi_result$RUN1_MiFishU_tax_tab))
RUN2_Mi_result$RUN2_MiFishU_tax_tab <- RUN2_Mi_result$RUN2_MiFishU_tax_tab %>% mutate(species = rownames(RUN2_Mi_result$RUN2_MiFishU_tax_tab))
RUN3_Mi_result$RUN3_MiFishU_tax_tab <- RUN3_Mi_result$RUN3_MiFishU_tax_tab %>% mutate(species = rownames(RUN3_Mi_result$RUN3_MiFishU_tax_tab))
RUN4_Mi_result$RUN4_MiFishU_tax_tab <- RUN4_Mi_result$RUN4_MiFishU_tax_tab %>% mutate(species = rownames(RUN4_Mi_result$RUN4_MiFishU_tax_tab))

RUN1_Mi_result$RUN1_MiFishU_tax_tab
RUN2_Mi_result$RUN2_MiFishU_tax_tab
RUN3_Mi_result$RUN3_MiFishU_tax_tab
RUN4_Mi_result$RUN4_MiFishU_tax_tab

# Merge them together by species
Mi_tax_tab <- RUN1_Mi_result$RUN1_MiFishU_tax_tab %>%
  full_join(RUN2_Mi_result$RUN2_MiFishU_tax_tab, by = "species") %>%
  full_join(RUN3_Mi_result$RUN3_MiFishU_tax_tab, by = "species") %>%
  full_join(RUN4_Mi_result$RUN4_MiFishU_tax_tab, by = "species") %>%
  replace(is.na(.), 0)# Replace NAs with 0 (species absent in that run)

Mi_tax_tab

Mi_tax_tab <- Mi_tax_tab %>%
  mutate(sp_primer = species) %>%
  mutate(species = str_replace(str_remove(species, "_MiFishU$"), "_", " "))

Mi_tax_tab

### Get tax info 

ref_tax_info <- read.csv("references.12s.miya.cleaned.v258.csv", header = TRUE)

ref_tax_info <- ref_tax_info %>%
  select(kingdom, phylum, class, order, family, genus, sciNameValid) %>%
  dplyr::rename(species = sciNameValid) %>%
  distinct()

Mi_tax_tab <- left_join(Mi_tax_tab, ref_tax_info, by = "species")

Mi_tax_tab <- Mi_tax_tab %>%
  column_to_rownames("sp_primer") %>%
  select(kingdom, phylum, class, order, family, genus, species)


Mi_tax_tab <- tax_table(as.matrix(Mi_tax_tab))

############################### sample data ####################################

Mi_sample_data <- read_xlsx("sample_data/met_MiFishU_sample_data.xlsx") %>%
  as.data.frame()                        # Convert to base R data frame

rownames(Mi_sample_data) <- Mi_sample_data$SampleID

Mi_sample_data <- sample_data(Mi_sample_data)

############################## Make phyloseq object ############################

Met_Mi <- phyloseq(tax_table(Mi_tax_tab), otu_table(Mi_otu_tab, taxa_are_rows = TRUE), Mi_sample_data)

Met_Mi

tax_table(Met_Mi)
otu_table(Met_Mi)
sample_data(Met_Mi)

save(Met_Mi, file = "working_phyloseq_data/Met_Mi.RData")









