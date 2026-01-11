# 11.01.25
# In this script I will get the number of consensus sequences for each species
# As per viva feedback 

library(writexl)
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(phyloseq)

############################ Get phyloseq with old species names ###############

load("working_phyloseq_data/Met_pseq.RData")

Met_pseq

tax_table(Met_pseq)

############################ Get phyloseq with updated species names ###########

load("working_phyloseq_data/Met_pseq_tax_updated.RData")

Met_pseq_tax_updated

tax_table(Met_Tel_pseq_tax_updated)
