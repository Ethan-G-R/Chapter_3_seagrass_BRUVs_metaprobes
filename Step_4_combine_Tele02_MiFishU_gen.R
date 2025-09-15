# 02.07.25
# In this script I will combine my Tel02 and MiFishU datasets at genus level
# <99% sim

library(writexl)
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(phyloseq)

load("working_phyloseq_data/Met_Tel_gen_99_97.RData")
load("working_phyloseq_data/Met_Mi_gen_99_97.RData")

tax_table(Met_Mi)
tax_table(Met_Tel)

Met_Mi_gen <- Met_Mi
Met_Tel_gen <- Met_Tel

# Combine the two phyloseq data sets
Met_pseq <- merge_phyloseq(Met_Mi_gen, Met_Tel_gen)

save(Met_pseq, file = "working_phyloseq_data/Met_gen_pseq.RData")

tax_table(Met_pseq)
otu_table(Met_pseq)
sample_data(Met_pseq)

# Remove primer column
sample_data(Met_pseq)$Primer <- NULL

# Update taxonomy

met_tax_table <- as.data.frame(tax_table(Met_pseq))

met_tax_table

### Add spp. for most genera but sp. for the genera which only have one species 
# in the UK

met_tax_table <- met_tax_table %>%
  
  # Identical 12S
  mutate(genus = if_else(genus == "Melanogrammus", "Merlangius", genus)) %>%
  mutate(species = if_else(genus == "Merlangius", "Merlangius spp.", species)) %>%

  # Name change
  mutate(species = if_else(species == "Gobiusculus spp.", "Pomatoschistus spp.", species)) %>%
  mutate(genus = if_else(genus == "Gobiusculus", "Pomatoschistus", genus)) %>%
  
  # Almost identical 12S
  mutate(genus = if_else(family == "Ammodytidae", "Ammodytidae gen.", genus)) %>%
  mutate(species = if_else(family == "Ammodytidae", "Ammodytidae spp.", species)) %>%
  
  # Almost identical 12S
  mutate(genus = if_else(family == "Clupeidae", "Clupeidae gen.", genus)) %>%
  mutate(species = if_else(family == "Clupeidae", "Clupeidae spp.", species))
 
met_tax_table

met_tax_table <- tax_table(as.matrix(met_tax_table))

met_tax_table

Met_pseq_tax_updated <- phyloseq(tax_table(met_tax_table), otu_table(Met_pseq), sample_data(Met_pseq))

save(Met_pseq_tax_updated, file = "working_phyloseq_data/Met_pseq_tax_updated_gen.RData")

# I aglom to species here since all the names are genus sp
Met_pseq_tax_updated_gen <- tax_glom(Met_pseq_tax_updated, taxrank = "species") 

tax_table(Met_pseq_tax_updated_gen)
otu_table(Met_pseq_tax_updated_gen)
sample_data(Met_pseq_tax_updated_gen)

# Convert to binary (i.e. make multiple detections using different primers
# from the same samples = 1)

otu_binary <- as(otu_table(Met_pseq_tax_updated_gen), "matrix")
otu_binary <- ifelse(otu_binary > 0, 1, 0)

# Recreate the OTU table
otu_table <- otu_table(otu_binary, taxa_are_rows = taxa_are_rows(Met_pseq_tax_updated_gen))

# Create a new phyloseq object with transformed OTU table
Met_pseq_tax_updated_gen <- phyloseq(otu_table, 
                                     tax_table(Met_pseq_tax_updated_gen), 
                                     sample_data(Met_pseq_tax_updated_gen))

save(Met_pseq_tax_updated_gen, file = "working_phyloseq_data/Met_pseq_tax_updated_glom_gen.RData")

Met_pseq_tax_updated_gen

