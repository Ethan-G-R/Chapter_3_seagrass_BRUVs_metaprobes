# 02.07.25
# In this script I will filter my BLAST results
# At genus leve < 99%

library(writexl)
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(phyloseq)

############################# Function #########################################

process_blast_summary <- function(file_path) {
  # Get the file stem (e.g., "RUN1_Tele02" from "BLAST_summaries/RUN1_Tele02_BLAST_summary.txt")
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
    filter(pident > 97 & pident < 99 & length > 150) %>%
    #    filter(pident > 98 & length > 150) %>%
    mutate(species = str_c(species, "_Tele02_gen"))
  
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
    mutate(species = str_remove(species, "_Tele02$"),
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


RUN1_Tel_result <- process_blast_summary("BLAST_summaries/RUN1_Tele02_BLAST_summary.txt")
RUN2_Tel_result <- process_blast_summary("BLAST_summaries/RUN2_Tele02_BLAST_summary.txt")
RUN3_Tel_result <- process_blast_summary("BLAST_summaries/RUN3_Tele02_BLAST_summary.txt")
RUN4_Tel_result <- process_blast_summary("BLAST_summaries/RUN4_Tele02_BLAST_summary.txt")

############################### Tele02 otu table ##############################

RUN1_Tel_result$RUN1_Tele02_otu_tab <- RUN1_Tel_result$RUN1_Tele02_otu_tab %>% mutate(species = rownames(RUN1_Tel_result$RUN1_Tele02_otu_tab))
RUN2_Tel_result$RUN2_Tele02_otu_tab <- RUN2_Tel_result$RUN2_Tele02_otu_tab %>% mutate(species = rownames(RUN2_Tel_result$RUN2_Tele02_otu_tab))
RUN3_Tel_result$RUN3_Tele02_otu_tab <- RUN3_Tel_result$RUN3_Tele02_otu_tab %>% mutate(species = rownames(RUN3_Tel_result$RUN3_Tele02_otu_tab))
RUN4_Tel_result$RUN4_Tele02_otu_tab <- RUN4_Tel_result$RUN4_Tele02_otu_tab %>% mutate(species = rownames(RUN4_Tel_result$RUN4_Tele02_otu_tab))

# Merge them together by species
Tel_otu_tab <- RUN1_Tel_result$RUN1_Tele02_otu_tab %>%
  full_join(RUN2_Tel_result$RUN2_Tele02_otu_tab, by = "species") %>%
  full_join(RUN3_Tel_result$RUN3_Tele02_otu_tab, by = "species") %>%
  full_join(RUN4_Tel_result$RUN4_Tele02_otu_tab, by = "species") %>%
  replace(is.na(.), 0) %>%  # Replace NAs with 0 (species absent in that run)
  tibble::column_to_rownames("species")

Tel_otu_tab

Tel_otu_tab <- as.matrix(Tel_otu_tab)

############################### Tele02 tax table ##############################

RUN1_Tel_result$RUN1_Tele02_tax_tab <- RUN1_Tel_result$RUN1_Tele02_tax_tab %>% mutate(species = rownames(RUN1_Tel_result$RUN1_Tele02_tax_tab))
RUN2_Tel_result$RUN2_Tele02_tax_tab <- RUN2_Tel_result$RUN2_Tele02_tax_tab %>% mutate(species = rownames(RUN2_Tel_result$RUN2_Tele02_tax_tab))
RUN3_Tel_result$RUN3_Tele02_tax_tab <- RUN3_Tel_result$RUN3_Tele02_tax_tab %>% mutate(species = rownames(RUN3_Tel_result$RUN3_Tele02_tax_tab))
RUN4_Tel_result$RUN4_Tele02_tax_tab <- RUN4_Tel_result$RUN4_Tele02_tax_tab %>% mutate(species = rownames(RUN4_Tel_result$RUN4_Tele02_tax_tab))

RUN1_Tel_result$RUN1_Tele02_tax_tab
RUN2_Tel_result$RUN2_Tele02_tax_tab
RUN3_Tel_result$RUN3_Tele02_tax_tab
RUN4_Tel_result$RUN4_Tele02_tax_tab

# Merge them together by species
Tel_tax_tab <- RUN1_Tel_result$RUN1_Tele02_tax_tab %>%
  full_join(RUN2_Tel_result$RUN2_Tele02_tax_tab, by = "species") %>%
  full_join(RUN3_Tel_result$RUN3_Tele02_tax_tab, by = "species") %>%
  full_join(RUN4_Tel_result$RUN4_Tele02_tax_tab, by = "species") %>%
  replace(is.na(.), 0)# Replace NAs with 0 (species absent in that run)

Tel_tax_tab

Tel_tax_tab <- Tel_tax_tab %>%
  mutate(sp_primer = species) %>%
  mutate(species = str_replace(str_remove(species, "_Tele02_gen$"), "_", " "))

Tel_tax_tab

### Get tax info 

ref_tax_info <- read.csv("references.12s.miya.cleaned.v258.csv", header = TRUE)

ref_tax_info <- ref_tax_info %>%
  select(kingdom, phylum, class, order, family, genus, sciNameValid) %>%
  dplyr::rename(species = sciNameValid) %>%
  distinct()

Tel_tax_tab <- left_join(Tel_tax_tab, ref_tax_info, by = "species")

Tel_tax_tab <- Tel_tax_tab %>%
  column_to_rownames("sp_primer") %>%
  select(kingdom, phylum, class, order, family, genus) %>%
  mutate(species = paste0(genus, " spp."))


Tel_tax_tab <- tax_table(as.matrix(Tel_tax_tab))

Tel_tax_tab

############################### sample data ####################################

Tel_sample_data <- read_xlsx("sample_data/met_Tele02_sample_data.xlsx") %>%
  as.data.frame()                        # Convert to base R data frame

rownames(Tel_sample_data) <- Tel_sample_data$SampleID

Tel_sample_data <- sample_data(Tel_sample_data)

############################## Make phyloseq object ############################

Met_Tel <- phyloseq(tax_table(Tel_tax_tab), otu_table(Tel_otu_tab, taxa_are_rows = TRUE), Tel_sample_data)

Met_Tel

tax_table(Met_Tel)
otu_table(Met_Tel)
sample_data(Met_Tel)

save(Met_Tel, file = "working_phyloseq_data/Met_Tel_gen_99_97.RData")









