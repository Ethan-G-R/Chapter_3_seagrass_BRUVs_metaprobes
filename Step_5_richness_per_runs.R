# 02.07.2025
# Ethan Ross

# In this script I will get taxa richness for my MiFishU and Tele02 datasets

load("working_phyloseq_data/Met_Tel.RData")
load("working_phyloseq_data/Met_Mi.RData")

## Update tax

met_tax_table_Tel <- as.data.frame(tax_table(Met_Tel))

met_tax_table_Tel

met_tax_table_Tel <- met_tax_table_Tel %>%
  
  # Identical 12S
  mutate(species = if_else(species == "Melanogrammus_aeglefinus", "Merlangius merlangus", species)) %>%
  mutate(genus = if_else(species == "Merlangius merlangus", "Merlangius", genus)) %>%
  
  # Almost identical 12S
  mutate(species = if_else(family == "Ammodytidae", "Ammodytidae spp.", species)) %>%
  mutate(genus = if_else(family == "Ammodytidae", "Ammodytidae gen.", genus)) %>%
  
  # Name change
  mutate(species = if_else(species == "Gobiusculus flavescens", "Pomatoschistus flavescens", species)) %>%
  mutate(genus = if_else(genus == "Gobiusculus", "Pomatoschistus", genus))

met_tax_table_Tel  <- tax_table(as.matrix(met_tax_table_Tel))

met_tax_table_Tel

Met_Tel_pseq_tax_updated <- phyloseq(tax_table(met_tax_table_Tel), otu_table(Met_Tel), sample_data(Met_Tel))

tax_table(Met_Tel_pseq_tax_updated) 

save(Met_Tel_pseq_tax_updated, file = "working_phyloseq_data/Met_Tel_pseq_tax_updated.RData")

################

met_tax_table_Mi <- as.data.frame(tax_table(Met_Mi))

met_tax_table_Mi

met_tax_table_Mi <- met_tax_table_Mi %>%
  
  # Identical 12S
  mutate(species = if_else(species == "Melanogrammus aeglefinus", "Merlangius merlangus", species)) %>%
  mutate(genus = if_else(species == "Merlangius merlangus", "Merlangius", genus)) %>%
  
  # Almost identical 12S
  mutate(species = if_else(family == "Ammodytidae", "Ammodytidae spp.", species)) %>%
  mutate(genus = if_else(family == "Ammodytidae", "Ammodytidae gen.", genus)) %>%
  
  # Name change
  mutate(species = if_else(species == "Gobiusculus flavescens", "Pomatoschistus flavescens", species)) %>%
  mutate(genus = if_else(genus == "Gobiusculus", "Pomatoschistus", genus))

met_tax_table_Mi  <- tax_table(as.matrix(met_tax_table_Mi ))

met_tax_table_Mi

Met_Mi_pseq_tax_updated <- phyloseq(tax_table(met_tax_table_Mi), otu_table(Met_Mi), sample_data(Met_Mi))

tax_table(Met_Mi_pseq_tax_updated)

save(Met_Mi_pseq_tax_updated, file = "working_phyloseq_data/Met_Mi_pseq_tax_updated.RData")

################################################################################

Met_Mi_sp_glom <- tax_glom(Met_Mi_pseq_tax_updated, taxrank = "species")
Met_Tel_sp_glom <- tax_glom(Met_Tel_pseq_tax_updated, taxrank = "species")

sample_data(Met_Mi_sp_glom)

sample_data(Met_Mi_sp_glom)$deployment <- as.character(sample_data(Met_Mi_sp_glom)$deployment)
sample_data(Met_Tel_sp_glom)$deployment <- as.character(sample_data(Met_Tel_sp_glom)$deployment)

Met_Mi_tab <- as.data.frame(otu_table(Met_Mi_sp_glom))

Met_Mi_tab

# 1. Convert to long format
Met_Mi_tab_long <- Met_Mi_tab %>%
  as.data.frame() %>%
  rownames_to_column("species") %>%
  pivot_longer(-species, names_to = "sample", values_to = "abundance")

# 2. Extract deployment number and assign new group names
Met_Mi_tab_long <- Met_Mi_tab_long %>%
  mutate(deployment_num = as.numeric(str_extract(sample, "(?<=Met_)\\d+")),
         group = paste0("deployment_", ceiling(deployment_num / 3)))

# 3. Summarise by species and deployment group
Met_Mi_grouped_sp <- Met_Mi_tab_long %>%
  group_by(species, group) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = abundance)

########

Met_Tel_tab <- as.data.frame(otu_table(Met_Tel_sp_glom))

Met_Tel_tab

# 1. Convert to long format
Met_Tel_tab_long <- Met_Tel_tab %>%
  as.data.frame() %>%
  rownames_to_column("species") %>%
  pivot_longer(-species, names_to = "sample", values_to = "abundance")

# 2. Extract deployment number and assign new group names
Met_Tel_tab_long <- Met_Tel_tab_long %>%
  mutate(deployment_num = as.numeric(str_extract(sample, "(?<=Met_)\\d+")),
         group = paste0("deployment_", ceiling(deployment_num / 3)))

# 3. Summarise by species and deployment group
Met_Tel_grouped_sp <- Met_Tel_tab_long %>%
  group_by(species, group) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = abundance)

# View result
Met_Tel_grouped_sp

Met_Mi_grouped_sp


################################################################################


Mi_counts_ge_1 <- Met_Mi_grouped_sp %>%
  pivot_longer(
    cols = starts_with("deployment_"), # Select all columns starting with 'deployment_'
    names_to = "deployment",           # Name the new key column 'deployment'
    values_to = "value"                # Name the new value column 'value'
  ) %>%
  filter(value >= 1) %>%               # Keep only rows where the value is 1 or greater
  group_by(deployment) %>%             # Group by deployment type
  summarise(count_ge_1 = n()) %>%
  mutate(marker = "MiFishU")# Count the number of rows in each group

# Print the result
print(Mi_counts_ge_1)

Tel_counts_ge_1 <- Met_Tel_grouped_sp %>%
  pivot_longer(
    cols = starts_with("deployment_"), # Select all columns starting with 'deployment_'
    names_to = "deployment",           # Name the new key column 'deployment'
    values_to = "value"                # Name the new value column 'value'
  ) %>%
  filter(value >= 1) %>%               # Keep only rows where the value is 1 or greater
  group_by(deployment) %>%             # Group by deployment type
  summarise(count_ge_1 = n()) %>%          # Count the number of rows in each group
  mutate(marker = "Tele02")

# Print the result
print(Tel_counts_ge_1)


taxa_counts <- rbind(Mi_counts_ge_1, Tel_counts_ge_1)

taxa_counts 

write_xlsx(taxa_counts, "counts/total_taxa_per_run.xlsx")

################################# Unique to each ################################

Mi_fish <- Met_Mi_grouped_sp$species
Mi_fish <- sub("_MiFishU$", "", Mi_fish)

Tel_fish <- Met_Tel_grouped_sp$species
Tel_fish <- sub("_Tele02$", "", Tel_fish)

Mi_fish
Tel_fish

unique_to_Mi_fish <- setdiff(Mi_fish, Tel_fish)
unique_to_Mi_fish

unique_to_Tel_fish <- setdiff(Tel_fish, Mi_fish)
unique_to_Tel_fish

shared_values <- intersect(Mi_fish, Tel_fish)
shared_values
