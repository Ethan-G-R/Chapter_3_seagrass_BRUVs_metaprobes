# 08.08.25
# Ethan Ross

# In this script I will carry out a restricted permutation PERMANOVA
# on my Metaprobe data

#install.packages("permute")

library(writexl)
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(permute)
library(vegan)


########################### RICHNESS LINEAR MODEL ##############################

# without BRUV data 

met_alpha_fishes <- read_xlsx("counts/alpha_diversity_per_metaprobe.xlsx")

met_alpha_fishes <- met_alpha_fishes %>%
  filter(sample_type == "Metaprobe")

m1 <- lm(Observed ~ deployment, data = met_alpha_fishes)

par(mfrow = c(2, 2))
plot(m1)

summary(m1)
anova(m1)


############################### PERMANOVA ######################################

load("working_phyloseq_data/met_bruv_pseq_glom_sp_with_gen.RData")

tax_table(met_bruv_pseq_sp)
otu_table(met_bruv_pseq_sp)
sample_data(met_bruv_pseq_sp)

samples_to_remove <- c("BRUV 1", "BRUV 2", "BRUV 3", "BRUV 4")

met_pseq_sp <- prune_samples(!(sample_names(met_bruv_pseq_sp) %in% samples_to_remove), met_bruv_pseq_sp)

tax_table(met_pseq_sp)
otu_table(met_pseq_sp)
sample_data(met_pseq_sp)

otu_data <- otu_table(met_pseq_sp)
otu_data <- otu_data[rowSums(otu_data) > 0, ]  # Remove rows with no data

# Make sample names rownames
otu_data <- t(otu_data)

otu_data

# Distance matrix
jaccard_dist <- vegdist(otu_data, method = "jaccard")

# Metadata
metadata <- as(sample_data(met_pseq_sp), "data.frame")

metadata

## Both work 

adonis2(jaccard_dist ~ deployment, data = metadata, permutations = 9999, strata = metadata$day)

adonis2(jaccard_dist ~ deployment, data = metadata, permutations = 9999)













