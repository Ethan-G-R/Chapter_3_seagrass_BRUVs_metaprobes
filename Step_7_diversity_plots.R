# 08.08.25
# Ethan Ross

library(writexl)
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(ggtext)

# In this script I will create diversity plots

load("working_phyloseq_data/met_bruv_pseq_glom_sp.RData")

tax_table(met_bruv_pseq_sp)
otu_table(met_bruv_pseq_sp)

tax_info <- as.data.frame(tax_table(met_bruv_pseq_sp))

tax_info <- tax_info %>% select(phylum, class, order, family, genus, species)

otu_binary <- as(otu_table(met_bruv_pseq_sp), "matrix")
otu_binary <- ifelse(otu_binary > 0, 1, 0)

heat_data <- as.data.frame(otu_binary)

heat_data

tax <- rownames(heat_data)

heat_data <- heat_data %>%
  mutate(`Met 1` = rowSums(select(., contains("Met 1")), na.rm = TRUE)) %>%
  mutate(`Met 2` = rowSums(select(., contains("Met 2")), na.rm = TRUE)) %>%
  mutate(`Met 3` = rowSums(select(., contains("Met 3")), na.rm = TRUE)) %>%
  mutate(`Met 4` = rowSums(select(., contains("Met 4")), na.rm = TRUE)) %>%
  select(`BRUV 1`, `Met 1`, `BRUV 2`, `Met 2`, `BRUV 3`, `Met 3`, `BRUV 4`, `Met 4`) %>%
  mutate(tax = tax)

heat_data$tax <- str_replace_all(heat_data$tax, "_MiFishU|_Tele02|_BRUV", "")
heat_data$tax <- str_replace_all(heat_data$tax, "_", " ")

heat_data <- cbind(heat_data, tax_info)

heat_data

## There is some hidden differences between Pomatoschistus flavescens
# this will resolve it
heat_data <- heat_data %>%
  mutate(species = if_else(str_detect(species, "flavescens"), 
                           "Pomatoschistus flavescens", 
                           species))

heat_data_long <- heat_data %>%
  pivot_longer(cols = matches("BRUV|Met"),  # Match columns with "BRUV" or "Met"
               names_to = c("method", "time"),
               names_pattern = "(.*) (.*)",
               values_to = "value")

heat_data_long <- heat_data_long %>%
  mutate(species = str_replace(species, " spp.", ""))

name_order <- heat_data_long %>%
  arrange(phylum, class, order, family, genus, species) %>%
  select(species) %>%
  unique()

name_order <- rev(name_order$species)

name_order

head(heat_data_long)

unique(heat_data_long$species)

# Filter data for method == "Met" for heatmap
heatmap_data <- heat_data_long %>%
  filter(method == "Met")

heatmap_data %>% filter(species == "Pomatoschistus flavescens")

heatmap_data

heat_sp <- ggplot(heatmap_data, aes(x = time, y = factor(species, levels = name_order), fill = value)) +
  geom_tile() +  # Create the heatmap
  scale_shape_manual(values = 4, name = "") +  # shape = 4 for 'x', blank title ("") or you can call it "Method"
  scale_fill_gradient(low = "white", high = "#38C682") +
  theme_classic() +
  guides(
    fill = guide_colorbar(order = 1),
    shape = guide_legend(order = 2)   # BRUV legend second
  ) +
  labs(title = "", fill = "Metaprobe\nDetection\nFrequency\n",
       x = "Deployment", y = "            Taxon\n\n") +
  geom_point(data = heat_data_long %>% filter(method == "BRUV" & value == 1), 
             aes(x = time, y = species, shape = "BRUV\nDetection"),  # map shape
             color = "black", size = 3, stroke = 1.5) +  # customize appearance
  theme(
    plot.title = element_text(size = 15),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    # Use element_markdown for conditional italicization
    axis.text.y = element_markdown(
      size = 12.5,
      # Use an `ifelse` statement or `case_when` to set the style
      # for each label. Labels not specified here will default to italic.
      # You might need to adjust 'name_order' to include these exact strings
      # if they are not already.
      # For "Labridae spp." and "Ammodytes spp." make them non-italic.
      # For all others, keep them italic.
      face = ifelse(levels(factor(heatmap_data$species, levels = name_order)) %in% c("Labridae", "Ammodytidae"), "plain", "italic")
    ),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.key.size = unit(1, "cm"),
    legend.position = "right",
    text = element_text(size = 25))

heat_sp

ggsave(filename = "plots/heat/heat_sp.png", plot = heat_sp,
       device = "png", dpi = 600, units = "mm", height = 200, width = 200)


ggsave(filename = "plots/heat/heat__long_sp.png", plot = heat_sp,
       device = "png", dpi = 600, units = "mm", height = 225, width = 200)

############################# Alpha diveristy ##################################

alpha_fishes <- estimate_richness(met_bruv_pseq_sp, measures = c("Observed"))

# 2. Convert to a data frame and add sample_type and deployment columns
alpha_fishes <- alpha_fishes %>%
  rownames_to_column(var = "sample_id") %>%
  mutate(
    sample_type = case_when(
      grepl("Met", sample_id) ~ "Metaprobe",
      grepl("BRUV", sample_id) ~ "BRUV",
      TRUE ~ NA_character_
    ),
    deployment = str_extract(sample_id, "[1-4]") # Extract the digit
  )

### Add day 
alpha_fishes <- alpha_fishes %>%
  mutate(
    Day = case_when(
      str_detect(sample_id, "1") | str_detect(sample_id, "2") ~ 1,
      str_detect(sample_id, "3") | str_detect(sample_id, "4") ~ 2,
      TRUE ~ NA_real_ # This is the "else" clause. Assign NA for any other cases.
      # Use NA_integer_ if you explicitly want integer NA.
    ),
    Day = as.factor(Day) # Convert the newly created Day column to a factor
  )



# 3. Convert grouping variables to factors
alpha_fishes$deployment <- as.factor(alpha_fishes$deployment)
alpha_fishes$sample_type <- as.factor(alpha_fishes$sample_type)
alpha_fishes

write_xlsx(alpha_fishes, "counts/alpha_diversity_per_metaprobe.xlsx")


#### Alpha plot

met_alpha_fishes <- read_xlsx("counts/alpha_diversity_per_metaprobe.xlsx")

met_alpha_fishes$deployment <- factor(met_alpha_fishes$deployment, levels = c("1", "2", "3", "4"))

met_alpha_fishes

# Calculate horizontal offsets for duplicates within each deployment
met_alpha_fishes_stacked <- met_alpha_fishes %>%
  group_by(deployment, Observed) %>%
  mutate(offset = (row_number() - mean(row_number())) * 0.15) %>%  # horizontal spacing
  ungroup() %>%
  mutate(deployment_offset = as.numeric(as.factor(deployment)) + offset)

alpha_point <- ggplot(met_alpha_fishes_stacked, 
                      aes(x = deployment_offset, y = Observed, 
                          colour = deployment, 
                          fill = deployment, 
                          shape = sample_type)) +
  geom_point(size = 5) +  
  # stat_summary(fun = mean, geom = "point", size = 3, shape = 4, color = "black", alpha = 0.5, stroke = 2) +
  # stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
  #               geom = "errorbar", width = 0.2, color = "black", alpha = 0.25, linewidth = 1) +
  theme_classic() +
  scale_x_continuous(breaks = 1:length(unique(met_alpha_fishes$deployment))) +
  scale_fill_manual(values = c(
    "1" = "#FF654B",
    "2" = "#EBCC2A",
    "3" = "#38C682",
    "4" = "#30788C"),
    labels = c(
      "1" = "1 (Day 1)",
      "2" = "2 (Day 1)",
      "3" = "3 (Day 2)",
      "4" = "4 (Day 2)")
  ) +
  scale_color_manual(values = c(
    "1" = "#FF654B",
    "2" = "#EBCC2A",
    "3" = "#38C682",
    "4" = "#30788C")
  ) +
  scale_shape_manual(
    values = c(
      "BRUV" = 17,  # triangle
      "Metaprobe" = 16  # circle
    )) +
  labs(x = "Deployment", y = "Species Richness", title = " ") +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.key.size = unit(1, "cm"),
    legend.position = "none"
  ) +
  ylim(0,17)

alpha_point

#################################### NMDS ######################################

load("working_phyloseq_data/met_bruv_pseq_glom_sp_with_gen.RData")

tax_table(met_bruv_pseq_sp)
otu_table(met_bruv_pseq_sp)
sample_data(met_bruv_pseq_sp)

# Calculate Jaccard distance
jaccard_dist <- phyloseq::distance(met_bruv_pseq_sp, method = "jaccard")
jaccard_dist

# Perform NMDS on the Jaccard distance matrix
jaccard_NMDS <- ordinate(met_bruv_pseq_sp, method = "NMDS", distance = jaccard_dist)

# Extract stress value
stress_value <- jaccard_NMDS$stress

# Create a label for the stress value
stress_label <- paste0("Stress: ", round(stress_value, 3))

sample_data(met_bruv_pseq_sp)

sample_data(met_bruv_pseq_sp)$deployment <- as.factor(sample_data(met_bruv_pseq_sp)$deployment)
sample_data(met_bruv_pseq_sp)$sample <- as.factor(sample_data(met_bruv_pseq_sp)$sample)
sample_data(met_bruv_pseq_sp)$day <- as.factor(sample_data(met_bruv_pseq_sp)$day)

met_bruv_pseq_sp

tax_table(met_bruv_pseq_sp)

# Plot NMDS
jaccard_NMDS_sp <- phyloseq::plot_ordination(met_bruv_pseq_sp, jaccard_NMDS,
                                             color = "deployment", shape = "SampleType") +
  geom_point(size = 5, alpha = 1, stroke = 0) +  # Set the size of the points
  #  geom_text(aes(label = sample_names(merged_pseq)), 
  #            size = 3, vjust = -1) +  # Add sample labels above points
  theme_classic() + # Optional: improve the theme
  #  stat_ellipse(aes(group = Sample_Type), level = 0.95) +
  ggtitle(paste0(jaccard_NMDS$method, stress_label)) +
  theme(
    plot.title = element_text(size = 15),,
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    legend.text = element_text(size = 15, lineheight = 2),
    legend.title = element_text(size = 20),
    legend.title.position = "top",
    legend.key.size = unit(1, "cm"),
    legend.position = "right",
    text = element_text(size = 20)) +
  labs(color = "Deployment", shape = "Sample Type") +
  # geom_text(aes(label = sample_names(met_bruv_pseq_sp)), 
  #            size = 4, vjust = -1.2, show.legend = FALSE)
  scale_shape_manual(
    values = c(
      "BRUV" = 17,  # solid circle
      "Metaprobe" = 16  # solid triangle
    )
  ) +
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
    ) 
  ) 
  # guides(
  #   color = guide_legend(ncol = 2, order = 1),
  #   shape = guide_legend(ncol = 2, order = 2)
  # )

#### THIS PLOT IS MADE FROM BOH GENUS AND SPECIES LEVEL ASSIGNMENTS

jaccard_NMDS_sp

################### Composit plot 

library(cowplot)

# Arrange the plots in a single column (ncol = 1)
composite_plot <- plot_grid(
  alpha_point,
  jaccard_NMDS_sp,
  ncol = 1,
  align = "v", # Align vertically
  axis = "l",   # Align on the left (y) axis
  rel_heights = c(1.5, 2.5),
  labels = c("(a)", "(b)"),
  label_size = 20,          # Make the labels bigger
  label_fontface = "plain"
)

composite_plot

ggsave(filename = "plots/composite_plots/alpha_with_NMDS_jacc_sp_and_gen.png", plot = composite_plot,
       device = "png", dpi = 600, units = "mm", height = 250, width = 200)

ggsave(filename = "plots/alpha_diversity/alpha_sp.png", plot = alpha_point,
       device = "png", dpi = 600, units = "mm", height = 125, width = 200)

ggsave(filename = "plots/beta_diversity/beta_NMDS_jacc_sp_and_gen.png", plot = jaccard_NMDS_sp,
       device = "png", dpi = 600, units = "mm", height = 200, width = 250)


























