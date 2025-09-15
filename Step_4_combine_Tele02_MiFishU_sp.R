# 02.07.25
# In this script I will combine my Tel02 and MiFishU datasets

library(writexl)
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(phyloseq)

load("working_phyloseq_data/Met_Tel.RData")
load("working_phyloseq_data/Met_Mi.RData")

Met_Mi
Met_Tel

# Combine the two phyloseq data sets
Met_pseq <- merge_phyloseq(Met_Mi, Met_Tel)

save(Met_pseq, file = "working_phyloseq_data/Met_pseq.RData")

tax_table(Met_pseq)
otu_table(Met_pseq)
sample_data(Met_pseq)

Met_pseq

# Remove primer column
sample_data(Met_pseq)$Primer <- NULL

met_tax_table <- as.data.frame(tax_table(Met_pseq))

met_tax_table

met_tax_table <- met_tax_table %>%
  
  # Identical 12S
  mutate(species = if_else(species == "Melanogrammus aeglefinus", "Merlangius merlangus", species)) %>%
  mutate(genus = if_else(species == "Merlangius merlangus", "Merlangius", genus)) %>%

  # Almost identical 12S
  mutate(species = if_else(family == "Ammodytidae", "Ammodytidae spp.", species)) %>%
  mutate(genus = if_else(family == "Ammodytidae", "Ammodytidae gen.", genus)) %>%
  
  # Name change
  mutate(species = if_else(species == "Gobiusculus flavescens", "Pomatoschistus flavescens", species)) %>%
  mutate(genus = if_else(genus == "Gobiusculus", "Pomatoschistus", genus))


met_tax_table

met_tax_table <- tax_table(as.matrix(met_tax_table))

met_tax_table

Met_pseq_tax_updated <- phyloseq(tax_table(met_tax_table), otu_table(Met_pseq), sample_data(Met_pseq))

save(Met_pseq_tax_updated, file = "working_phyloseq_data/Met_pseq_tax_updated.RData")

Met_pseq_tax_updated_sp <- tax_glom(Met_pseq_tax_updated, taxrank = "species") 

tax_table(Met_pseq_tax_updated_sp)
otu_table(Met_pseq_tax_updated_sp)
sample_data(Met_pseq_tax_updated_sp)

# Convert to binary (i.e. make multiple detections using different primers
# from the same samples = 1)

otu_binary <- as(otu_table(Met_pseq_tax_updated_sp), "matrix")
otu_binary <- ifelse(otu_binary > 0, 1, 0)

# Recreate the OTU table
otu_table <- otu_table(otu_binary, taxa_are_rows = taxa_are_rows(Met_pseq_tax_updated_sp))

# Create a new phyloseq object with transformed OTU table
Met_pseq_tax_updated_sp <- phyloseq(otu_table, 
                                    tax_table(Met_pseq_tax_updated_sp), 
                                    sample_data(Met_pseq_tax_updated_sp))

save(Met_pseq_tax_updated_sp, file = "working_phyloseq_data/Met_pseq_tax_updated_glom_sp.RData")

load("working_phyloseq_data/Met_pseq_tax_updated_glom_sp.RData")

tax_table(Met_pseq_tax_updated_sp)
otu_table(Met_pseq_tax_updated_sp)
