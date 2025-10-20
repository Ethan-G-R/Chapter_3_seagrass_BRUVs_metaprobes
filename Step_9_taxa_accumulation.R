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

################################################################################
############################ with iNEXT ########################################
# 1. Install the iNEXT package
#install.packages("iNEXT")

library(iNEXT)

######## Get my data

load("working_phyloseq_data/met_bruv_pseq_glom_sp.RData")

tax_table(met_bruv_pseq_sp)
otu_table(met_bruv_pseq_sp)
sample_data(met_bruv_pseq_sp)

meta_data <- data.frame(t(otu_table(met_bruv_pseq_sp)))

meta_data <- meta_data %>%
  rownames_to_column(var = "Sample")

################################################################################

meta_sample_data <- data.frame(sample_data(met_bruv_pseq_sp))

meta_sample_data <- meta_sample_data %>%
  select(deployment, SampleType)

meta_sample_data <- meta_sample_data %>%
  rownames_to_column(var = "Sample")


#################### Combine ###################################################

meta_data_full <- full_join(meta_sample_data, meta_data, by = "Sample")

##### Species accumulation curves #####
# Function to convert a column to binary based on a threshold
binary_converter <- function(column, threshold) {
  binary_data <- ifelse(column >= threshold, 1, 0)
  return(binary_data)
}

str(meta_data_full)

meta_data_full_deployment <- meta_data_full %>%
  filter(SampleType == "Metaprobe") %>%
  select(-SampleType)

str(meta_data_full_deployment)

# Apply the binary_converter function to all columns
AccCurve_12S <- meta_data_full_deployment[,c(3:28)]
thresholds <- 1
AccCurve_12S <- data.frame(lapply(AccCurve_12S, binary_converter, threshold = thresholds))
AccCurve_12S <- data.frame(meta_data_full_deployment[,c(1,2)], AccCurve_12S)

## Accumulation curves 12S (iNEXT format)
# the second column (deployment) is removed since all values are the same after
# which selects for specific deployment
# Row names are then turned into a column

Accdep1 <- AccCurve_12S[which(AccCurve_12S$deployment == '1'),-2]
Accdep1 <- as.data.frame(t(Accdep1))
colnames(Accdep1) <- Accdep1[1,]
Accdep1 <- Accdep1[-1,]
Accdep1$species <- rownames(Accdep1)
rownames(Accdep1) <- NULL

Accdep2 <- AccCurve_12S[which(AccCurve_12S$deployment == '2'),-2]
Accdep2 <- as.data.frame(t(Accdep2))
colnames(Accdep2) <- Accdep2[1,]
Accdep2 <- Accdep2[-1,]
Accdep2$species <- rownames(Accdep2)
rownames(Accdep2) <- NULL

Accdep3 <- AccCurve_12S[which(AccCurve_12S$deployment == '3'),-2]
Accdep3 <- as.data.frame(t(Accdep3))
colnames(Accdep3) <- Accdep3[1,] 
Accdep3 <- Accdep3[-1,]
Accdep3$species <- rownames(Accdep3)
rownames(Accdep3) <- NULL

Accdep4 <- AccCurve_12S[which(AccCurve_12S$deployment == '4'),-2]
Accdep4 <- as.data.frame(t(Accdep4))
colnames(Accdep4) <- Accdep4[1,]
Accdep4 <- Accdep4[-1,]
Accdep4$species <- rownames(Accdep4)
rownames(Accdep4) <- NULL

# Turns fourth column (species names0 back into row names and makes sure
# the numbers are all integers

dep1 <- as.matrix(apply(Accdep1[,-4],2,as.integer))
row.names(dep1) <- Accdep1[,4]
dep2 <- as.matrix(apply(Accdep2[,-4],2,as.integer))
row.names(dep2) <- Accdep2[,4]
dep3 <- as.matrix(apply(Accdep3[,-4],2,as.integer))
row.names(dep3) <- Accdep3[,4]
dep4 <- as.matrix(apply(Accdep4[,-4],2,as.integer))
row.names(dep4) <- Accdep4[,4]

################################################################################

AccCurve12S = list(dep1 = dep1, dep2 = dep2, dep3 = dep3, dep4 = dep4)

str(AccCurve12S)

# End point is how many samples the accumulation curve is estimated up to

Acc12S_next_met <- iNEXT(AccCurve12S, datatype="incidence_raw", endpoint=12)

Acc12S_met <- ggiNEXT(Acc12S_next_met) +
  scale_fill_manual(values = c(
    "dep1" = "#FF654B",
    "dep2" = "#EBCC2A",
    "dep3" = "#38C682",
    "dep4" = "#30788C"),
    labels = c(
      "dep1" = "1",
      "dep2" = "2",
      "dep3" = "3",
      "dep4" = "4"
    )) +
  scale_color_manual(values = c(
    "dep1" = "#FF654B",
    "dep2" = "#EBCC2A",
    "dep3" = "#38C682",
    "dep4" = "#30788C"),
    labels = c(
      "dep1" = "1",
      "dep2" = "2",
      "dep3" = "3",
      "dep4" = "4"
    )) +
  scale_shape_manual(values = c(
    "dep1" = 15, # Solid circle
    "dep2" = 16, 
    "dep3" = 18, 
    "dep4" = 17),
    labels = c(
      "dep1" = "1",
      "dep2" = "2",
      "dep3" = "3",
      "dep4" = "4"
    )) +  
  xlab("Number of samples") + ylab("Species richness") +
  scale_x_continuous(breaks = seq(from = 3, to = 12, by = 3)) +
  scale_y_continuous(limits = c(2, 30)) +
  theme_classic() + 
  guides(
    color = guide_legend(title = "Metaprobe deployment", order = 1),
    fill = guide_legend(title = "Metaprobe deployment", order = 1),
    shape = guide_legend(title = "Metaprobe deployment", order = 1),
    # Target the guide responsible for the line/method (usually linetype)
    linetype = guide_legend(title = NULL, order = 2) 
  ) +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 17.5),
    legend.title.position = "top",
    legend.key.size = unit(1, "cm"),
    legend.position = "bottom"
  )


Acc12S_met

ggsave(filename = "plots/accum/Metaprobe_rare.png", plot = Acc12S_met,
       device = "png", dpi = 600, units = "mm", height = 125, width = 150)

################################################################################

### Metaprobe vs BRUV

meta_data_full_method <- meta_data_full %>%
  select(-deployment)

str(meta_data_full_method)

# Apply the binary_converter function to all columns
AccCurve_12S <- meta_data_full_method[,c(3:28)]
thresholds <- 1
AccCurve_12S <- data.frame(lapply(AccCurve_12S, binary_converter, threshold = thresholds))
AccCurve_12S <- data.frame(meta_data_full_method[,c(1,2)], AccCurve_12S)

## Accumulation curves 12S
Accmet <- AccCurve_12S[which(AccCurve_12S$SampleType == 'Metaprobe'),-2]
Accmet <- as.data.frame(t(Accmet))
colnames(Accmet) <- Accmet[1,]
Accmet <- Accmet[-1,]
Accmet$species <- rownames(Accmet)
rownames(Accmet) <- NULL

AccBRUV <- AccCurve_12S[which(AccCurve_12S$SampleType == 'BRUV'),-2]
AccBRUV <- as.data.frame(t(AccBRUV))
colnames(AccBRUV) <- AccBRUV[1,]
AccBRUV <- AccBRUV[-1,]
AccBRUV$species <- rownames(AccBRUV)
rownames(AccBRUV) <- NULL

met <- as.matrix(apply(Accmet[,-13],2,as.integer))
row.names(met) <- Accmet[,13]
BRUV <- as.matrix(apply(AccBRUV[,-5],2,as.integer))
row.names(BRUV) <- AccBRUV[,5]

################################################################################

AccCurve12S = list(met = met, BRUV = BRUV)

str(AccCurve12S)

Acc12S_next_met_vs_BRUV <- iNEXT(AccCurve12S, datatype="incidence_raw", endpoint=20)

Acc12S_met_vs_BRUV <- ggiNEXT(Acc12S_next_met_vs_BRUV) +
  scale_fill_manual(values = c(
    "met" = "#38C682",
    "BRUV" = "black"),
    labels = c(
      "met" = "Metaprobe",
      "BRUV" = "BRUV")) +
  scale_color_manual(values = c(
    "met" = "#38C682",
    "BRUV" = "black"),
    labels = c(
      "met" = "Metaprobe",
      "BRUV" = "BRUV")) +
  scale_shape_manual(values = c(
    "met" = 16,
    "BRUV" = 17),
    labels = c(
      "met" = "Metaprobe",
      "BRUV" = "BRUV")) +
  xlab("Number of samples") + ylab("Species richness") +
  scale_y_continuous(limits = c(2, 30)) +
  theme_classic() +
  guides(
    color = guide_legend(title = "Detection method", order = 1),
    fill = guide_legend(title = "Detection method", order = 1),
    shape = guide_legend(title = "Detection method", order = 1),
    # Target the guide responsible for the line/method (usually linetype)
    linetype = guide_legend(title = NULL, order = 2) 
  ) +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 17.5),
    legend.title.position = "top",
    legend.key.size = unit(1, "cm"),
    legend.position = "bottom"
  ) 


Acc12S_met_vs_BRUV

ggsave(filename = "plots/accum/Metaprobe_vs_BRUV_rare.png", plot = Acc12S_met_vs_BRUV,
       device = "png", dpi = 600, units = "mm", height = 125, width = 150)

################################################################################


AccCurve12S_all = list(met = met, BRUV = BRUV, dep1 = dep1, dep2 = dep2, dep3 = dep3, dep4 = dep4)

str(AccCurve12S_all)

Acc12S_next_all <- iNEXT(AccCurve12S_all, datatype="incidence_raw", endpoint=20)

Acc12S_all <- ggiNEXT(Acc12S_next_all) +
  scale_fill_manual(values = c(
    "met" = "#7F7F7F",
    "BRUV" = "black",
    "dep1" = "#FF654B",
    "dep2" = "#EBCC2A",
    "dep3" = "#38C682",
    "dep4" = "#30788C")) +
  scale_color_manual(values = c(
    "met" = "#7F7F7F",
    "BRUV" = "black",
    "dep1" = "#FF654B",
    "dep2" = "#EBCC2A",
    "dep3" = "#38C682",
    "dep4" = "#30788C")) +
  scale_shape_manual(values = c(
    "met" = 16,
    "BRUV" = 17, # Solid triangle
    "dep1" = 16, # Solid circle
    "dep2" = 16, 
    "dep3" = 16, 
    "dep4" = 16  
  )) +
  xlab("Number of Samples") + ylab("Species richness") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.key.size = unit(1, "cm"),
    legend.position = "bottom"
  ) 

Acc12S_all










