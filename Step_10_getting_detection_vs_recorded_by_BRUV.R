# 06.08.25
# Ethan Ross

# I will get binary data of presence absence across deployments in prep for a 
# GLMM to test if detection by BRUV results in higher detection frequency

library(lme4)
library(tidyverse)
library(writexl)
library(readxl)
library(phyloseq)


load("working_phyloseq_data/Met_pseq_tax_updated_glom_sp.RData")

################################################################################

Met_pseq_tax_updated_sp

met_tax <- as.data.frame(otu_table(Met_pseq_tax_updated_sp))

met_tax

str(met_tax)

met_tax_combined <- met_tax %>%
  # Convert rownames to a column
  rownames_to_column(var = "rowname") %>%
  
  # Pivot longer for tidy reshaping
  pivot_longer(
    cols = -rowname,
    names_to = "met_group",
    values_to = "value"
  ) %>%
  
  # Extract the numeric part after "Met_"
  mutate(met_num = str_extract(met_group, "(?<=Met_)[0-9]+")) %>%
  
  # Group by row and met group number
  group_by(rowname, met_num) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = "drop") %>%
  
  # Pivot wider to restore to wide format
  pivot_wider(
    names_from = met_num,
    names_prefix = "Met_",
    values_from = value
  ) %>%
  
  # Restore rownames
  column_to_rownames(var = "rowname")

met_tax_combined

met_tax_combined <- met_tax_combined %>%
  mutate(across(everything(), ~ ifelse(. > 1, 1, .)))

met_tax_combined

################################################################################

sample_data <- as.matrix((sample_data(Met_pseq_tax_updated_sp)))

sample_data <- as.data.frame(sample_data) %>%
  select(sample, deployment, day) %>%
  mutate(sample_name = rownames(sample_data))

sample_data

sample_data_combined <- sample_data %>%
  mutate(sample_name = str_sub(sample_name, 1, -2)) %>%
  distinct()

sample_data_combined

################################################################################

recorded_sp <- c("Pollachius_virens_MiFishU", "Gobiusculus_flavescens_MiFishU", 
                 "Melanogrammus_aeglefinus_MiFishU", "Pollachius_pollachius_MiFishU",
                 "Gadus_morhua_Tele02", "Ctenolabrus_rupestris_Tele02",
                 "Hyperoplus_immaculatus_Tele02")

# Create a vector of species names (or IDs)
species_names <- rownames(met_tax)

recorded_sp_BRUV <- as.data.frame(species_names) %>%
  mutate(BRUV = ifelse(species_names %in% recorded_sp, 1, 0))

recorded_sp_BRUV

################################################################################

library(lme4)

str(met_tax)
str(sample_data)
str(recorded_sp_BRUV)

met_tax
met_tax_combined
sample_data
sample_data_combined
recorded_sp_BRUV

detection_df <- as.data.frame(met_tax_combined) %>%
  tibble::rownames_to_column(var = "species_names")

long_detection_data <- detection_df %>%
  pivot_longer(
    cols = starts_with("Met"),
    names_to = "sample_name",
    values_to = "detection"
  )

analysis_data <- left_join(long_detection_data, sample_data_combined, by = "sample_name") %>%
  left_join(recorded_sp_BRUV, by = "species_names")

analysis_data$day <- as.factor(analysis_data$day)
analysis_data$deployment <- as.factor(analysis_data$deployment)
analysis_data$sample <- as.factor(analysis_data$sample)
analysis_data$BRUV <- as.factor(analysis_data$BRUV)

str(analysis_data)

write_xlsx(analysis_data, "detection_frequency_data/metaprobe_detection_frequency_data.xlsx")

### The model

library(lme4)
library(readxl)

# The model formula is structured as follows:
# detection ~ known_present        -> `known_present` is a fixed effect
#           + (1 | day)            -> `day` is a random effect (intercept varies by day)
#           + (1 | deployment)     -> `deployment` is a random effect (intercept varies by deployment)
#
# The `family = binomial(link = "logit")` is used because the response variable
# (`detection`) is binary (0 or 1).

model <- glmer(
  detection ~ BRUV + (1 | deployment),
  data = analysis_data,
  family = binomial(link = "logit")
)

summary(model)

null_model <- glmer(
  detection ~ (1 | day) + (1 | deployment),
  data = analysis_data,
  family = binomial(link = "logit")
)

anova(model, null_model, test = "Chisq")














