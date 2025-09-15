# 03.07.25
# Ethan Ross

# In this script I will create species accumulation curves

library(phyloseq)
library(microbiome)
library(ggplot2)
library(dplyr)
library(stringr)
library(writexl)
library(readxl)
library(vegan)
library(tidyverse)

load("working_phyloseq_data/met_bruv_pseq_glom_sp.RData")

# Convert rownames to SampleID and tidy data
otu_table_df <- otu_table(met_bruv_pseq_sp) %>%
  t() %>%
  as.data.frame() %>%
  mutate(SampleID = rownames(.))

# Set the desired order and create the Group column
otu_table_df <- otu_table_df %>%
  mutate(Group = ifelse(grepl("Met", SampleID), "Metaprobe", "BRUV"))

desired_order <- c(
  "Met 1a", "Met 1b", "Met 1c",
  "Met 2a", "Met 2b", "Met 2c",
  "Met 3a", "Met 3b", "Met 3c",
  "Met 4a", "Met 4b", "Met 4c",
  "BRUV 1", "BRUV 2", "BRUV 3", "BRUV 4"
)

otu_table_df

# Set SampleID as a factor with the desired order
otu_table_df <- otu_table_df %>%
  mutate(SampleID = factor(SampleID, levels = desired_order)) %>%
  arrange(SampleID)

# Split data into Metaprobe and BRUV
otu_metaprobe <- otu_table_df %>% filter(Group == "Metaprobe") %>% select(-SampleID, -Group)
otu_bruv <- otu_table_df %>% filter(Group == "BRUV") %>% select(-SampleID, -Group)

# Run species accumulation for both groups
accum_metaprobe <- specaccum(otu_metaprobe > 0, method = "collector")
accum_bruv <- specaccum(otu_bruv > 0, method = "collector")

# Combine the results into a single data frame
accum_combined <- bind_rows(
  data.frame(Sample = 1:length(accum_metaprobe$sites), Richness = accum_metaprobe$richness, Group = "Metaprobe", SampleID = otu_table_df$SampleID[otu_table_df$Group == "Metaprobe"]),
  data.frame(Sample = 1:length(accum_bruv$sites), Richness = accum_bruv$richness, Group = "BRUV", SampleID = otu_table_df$SampleID[otu_table_df$Group == "BRUV"])
)

# Add deployment info
deployment_day <- c("1", "1", "1", "2", "2", "2", "3", "3", "3", "4", "4", "4",
                    "1", "2", "3", "4")

deployment_day <- data.frame(deployment = deployment_day)

accum_combined <- cbind(accum_combined, deployment_day)


accum_combined

accum_sp <- ggplot(accum_combined, aes(x = Sample, y = Richness, color = deployment, shape = Group)) +
  geom_line(size = 1) +
  geom_point(size = 5) +
  scale_shape_manual(values = c("Metaprobe" = 16, "BRUV" = 17)) +  # 16 is circle, 17 is triangle
  scale_color_manual(values = c(
    "1" = "#FF654B",
    "2" = "#EBCC2A",
    "3" = "#38C682",
    "4" = "#30788C"),
    labels = c(
      "1" = "1 (Day 1)",
      "2" = "2 (Day 1)",
      "3" = "3 (Day 2)",
      "4" = "4 (Day 2)"
    )) +
  #  geom_text(aes(label = SampleID), hjust = 0, vjust = 1.5, size = 4, check_overlap = TRUE) +
  geom_text(aes(label = SampleID), hjust = 1, vjust = -1, size = 4, check_overlap = TRUE) +
  scale_x_continuous(
    limits = c(0, 12),  # Set the x-axis range from 0 to 12
    breaks = 0:12) +       # Set the breaks to be every integer between 0 and 12
  labs(x = "Number of Samples", y = "Accumulated Taxonomic Richness", title = "", shape = "Sample\nMethod") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.key.size = unit(1, "cm"),
    legend.position = "right",
    text = element_text(size = 25)
  ) +
  guides(
    shape = guide_legend(override.aes = list(color = c("#38C682", "#30788C"), linewidth = 5)),  # Matching color to shape
    color = "none"  # Remove the color legend
  )

accum_sp

accum_combined

accum_sp <- ggplot(accum_combined, aes(x = Sample, y = Richness, colour = deployment, shape = Group)) +
  geom_line(aes(group = Group), color = "grey80", size = 1) +
  geom_line(size = 1) +
  geom_point(size = 5) +
  scale_shape_manual(values = c("Metaprobe" = 16, "BRUV" = 17)) +
  scale_colour_manual(values = c(
    "1" = "#FF654B",
    "2" = "#EBCC2A",
    "3" = "#38C682",
    "4" = "#30788C"),
    labels = c(
      "1" = "1 (Day 1)",
      "2" = "2 (Day 1)",
      "3" = "3 (Day 2)",
      "4" = "4 (Day 2)"
    )) +
  #  geom_text(aes(label = SampleID), hjust = 1, vjust = -1, size = 4, check_overlap = TRUE) +
  scale_x_continuous(
    limits = c(0, 12),
    breaks = 0:12) +
  ylim(0,24) +
  labs(x = "Number of Samples", y = "Accumulated Taxonomic Richness", title = "", shape = "Sample Type",
       colour = "Deployment") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.key.size = unit(1, "cm"),
    legend.position = "right",
    text = element_text(size = 25)
  )

accum_sp

ggsave(filename = "plots/accum/accum_sp_col.png", plot = accum_sp,
       device = "png", dpi = 600, units = "mm", height = 150, width = 200)











