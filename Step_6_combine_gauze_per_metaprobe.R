# 02.07.2025
# In this script I will combine all gauze samples from the same metaprobe for
# plotting

# I will make a species only phyloseq object (only >99% sim)

# And a species nd genus phyloseq object (all >97% sim)

library(writexl)
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(phyloseq)
library(microbiome)
library(ggplot2)

########################## SPECIES LEVEL ONLY ##################################

load("working_phyloseq_data/Met_pseq_tax_updated_glom_sp.RData")

##################### Add genus hits to species hits ###########################

# Check genera have been added

tax_table(Met_pseq_tax_updated_sp)
otu_table(Met_pseq_tax_updated_sp)
sample_data(Met_pseq_tax_updated_sp)

# Make sample data sample into character
sample_data(Met_pseq_tax_updated_sp)$sample <- as.character(sample_data(Met_pseq_tax_updated_sp)$sample)

# The sample data will become NA after merging
merged_met_pseq <- merge_samples(Met_pseq_tax_updated_sp, group = "sample")

sample_data(merged_met_pseq)$SampleType <- "Metaprobe"

otu_table(merged_met_pseq)
sample_data(merged_met_pseq)

##### Change sample names

met_sample_data <- sample_data(merged_met_pseq)

rownames(met_sample_data) <- c("Met 1a", "Met 4a", "Met 4b", "Met 4c", "Met 1b",
                               "Met 1c", "Met 2a", "Met 2b", "Met 2c", "Met 3a",
                               "Met 3b", "Met 3c")

met_otu <- as.data.frame(otu_table(merged_met_pseq))


rownames(met_otu) <- c("Met 1a", "Met 4a", "Met 4b", "Met 4c", "Met 1b",
                       "Met 1c", "Met 2a", "Met 2b", "Met 2c", "Met 3a",
                       "Met 3b", "Met 3c")

merged_met_pseq <- phyloseq(sample_data(met_sample_data),
                            otu_table(met_otu, taxa_are_rows = FALSE),
                            tax_table(merged_met_pseq))

merged_met_pseq

otu_table(merged_met_pseq)
tax_table(merged_met_pseq)
sample_data(merged_met_pseq)


sample_data(merged_met_pseq)$deployment <- as.factor(sample_data(merged_met_pseq)$deployment)
sample_data(merged_met_pseq)$sample <- as.factor(sample_data(merged_met_pseq)$sample)

##################### Make phyloseq object for BRUV ############################

bruv_tax <- as.data.frame(read_xlsx("BRUV_data.xlsx", sheet = 1))
bruv_otu <- as.data.frame(read_xlsx("BRUV_data.xlsx", sheet = 2))
bruv_sam <-  as.data.frame(read_xlsx("BRUV_data.xlsx", sheet = 3))

rownames(bruv_tax) <- bruv_tax$unique
bruv_tax <- as.matrix(bruv_tax %>% select(-unique))

rownames(bruv_otu) <- bruv_otu$unique
bruv_otu <- as.matrix(bruv_otu %>% select(-unique))

rownames(bruv_sam) <- colnames(bruv_otu)

bruv_tax
bruv_otu
bruv_sam

bruv_pseq <- phyloseq(tax_table(bruv_tax), 
                      otu_table(bruv_otu, taxa_are_rows = TRUE),
                      sample_data(bruv_sam))

bruv_pseq

save(bruv_pseq, file = "working_phyloseq_data/bruv_pseq.RData")

##################### Combine BRUV and Metaprobe data ##########################

met_bruv_pseq <- merge_phyloseq(merged_met_pseq, bruv_pseq)

tax_table(met_bruv_pseq)
otu_table(met_bruv_pseq)
sample_data(met_bruv_pseq)

otu_transposed <- t(otu_table(met_bruv_pseq))  # transpose the matrix
otu_table(met_bruv_pseq) <- otu_table(otu_transposed, taxa_are_rows = TRUE)

##################### Agglomerate to species level #############################

met_bruv_pseq_sp <- tax_glom(met_bruv_pseq, taxrank = "species") 

save(met_bruv_pseq_sp, file = "working_phyloseq_data/met_bruv_pseq_glom_sp.RData")

tax_table(met_bruv_pseq_sp)
otu_table(met_bruv_pseq_sp)
sample_data(met_bruv_pseq_sp)

as.data.frame(tax_table(met_bruv_pseq_sp)) %>% 
  arrange(kingdom, phylum, class, order, family, genus, species)

################################################################################
########################## SPECIES AND GENUS LEVEL #############################
################################################################################

load("working_phyloseq_data/Met_pseq_tax_updated_glom_sp.RData")
load("working_phyloseq_data/Met_pseq_tax_updated_glom_gen.RData")

### Check differences

tax_table(Met_pseq_tax_updated_gen)
tax_table(Met_pseq_tax_updated_sp)


# # Select genera not in species phyloseq
# gen_to_keep <- c("Taurulus_bubalis_MiFishU", 
#                  "Conger_erebennus_Tele02", "Symphodus_melops_MiFishU")
# 
# Met_pseq_tax_updated_gen <- prune_taxa(gen_to_keep, Met_pseq_tax_updated_gen)

otu_table(Met_pseq_tax_updated_gen)
tax_table(Met_pseq_tax_updated_gen)
sample_data(Met_pseq_tax_updated_gen)

otu_table(Met_pseq_tax_updated_sp)
tax_table(Met_pseq_tax_updated_sp)
sample_data(Met_pseq_tax_updated_sp)

##################### Add genus hits to species hits ###########################

Met_pseq_tax_updated_sp <- merge_phyloseq(Met_pseq_tax_updated_sp, 
                                          Met_pseq_tax_updated_gen)

# Check genera have been added

tax_table(Met_pseq_tax_updated_sp)
otu_table(Met_pseq_tax_updated_sp)
sample_data(Met_pseq_tax_updated_sp)

# Make sample data sample into character
sample_data(Met_pseq_tax_updated_sp)$sample <- as.character(sample_data(Met_pseq_tax_updated_sp)$sample)

# The sample data will become NA after merging
merged_met_pseq <- merge_samples(Met_pseq_tax_updated_sp, group = "sample")

sample_data(merged_met_pseq)$SampleType <- "Metaprobe"

otu_table(merged_met_pseq)
sample_data(merged_met_pseq)
tax_table(merged_met_pseq)

##### Change sample names

met_sample_data <- sample_data(merged_met_pseq)

rownames(met_sample_data) <- c("Met 1a", "Met 4a", "Met 4b", "Met 4c", "Met 1b",
                               "Met 1c", "Met 2a", "Met 2b", "Met 2c", "Met 3a",
                               "Met 3b", "Met 3c")

met_otu <- as.data.frame(otu_table(merged_met_pseq))

rownames(met_otu) <- c("Met 1a", "Met 4a", "Met 4b", "Met 4c", "Met 1b",
                       "Met 1c", "Met 2a", "Met 2b", "Met 2c", "Met 3a",
                       "Met 3b", "Met 3c")

merged_met_pseq <- phyloseq(sample_data(met_sample_data),
                            otu_table(met_otu, taxa_are_rows = FALSE),
                            tax_table(merged_met_pseq))

merged_met_pseq

otu_table(merged_met_pseq)
tax_table(merged_met_pseq)
sample_data(merged_met_pseq)


sample_data(merged_met_pseq)$deployment <- as.factor(sample_data(merged_met_pseq)$deployment)
sample_data(merged_met_pseq)$sample <- as.factor(sample_data(merged_met_pseq)$sample)

##################### Get phyloseq object for BRUV ############################

load("working_phyloseq_data/bruv_pseq.RData")


##################### Combine BRUV and Metaprobe data ##########################

met_bruv_pseq <- merge_phyloseq(merged_met_pseq, bruv_pseq)

tax_table(met_bruv_pseq)
otu_table(met_bruv_pseq)
sample_data(met_bruv_pseq)

otu_transposed <- t(otu_table(met_bruv_pseq))  # transpose the matrix
otu_table(met_bruv_pseq) <- otu_table(otu_transposed, taxa_are_rows = TRUE)

##################### Agglomerate to species level #############################

met_bruv_pseq_sp <- tax_glom(met_bruv_pseq, taxrank = "species") 

save(met_bruv_pseq_sp, file = "working_phyloseq_data/met_bruv_pseq_glom_sp_with_gen.RData")

tax_table(met_bruv_pseq_sp)
otu_table(met_bruv_pseq_sp)
sample_data(met_bruv_pseq_sp)

as.data.frame(tax_table(met_bruv_pseq_sp)) %>% 
  arrange(kingdom, phylum, class, order, family, genus, species)







































