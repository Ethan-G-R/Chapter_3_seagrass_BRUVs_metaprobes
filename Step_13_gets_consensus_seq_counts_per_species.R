# 11.01.25
# In this script I will get the number of consensus sequences for each species
# As per viva feedback 

library(writexl)
library(readxl)
library(phyloseq)
library(tidyverse)
library(ggtext)

############################ Get phyloseq with updated species names ###########

load("working_phyloseq_data/Met_pseq_tax_updated.RData")

Met_pseq_tax_updated

tax_Met_updated <- as.data.frame(as(tax_table(Met_pseq_tax_updated), "matrix"))

tax_Met_updated

tax_Met_updated <- tax_Met_updated %>%
  mutate(BLAST_species_primer = row.names(tax_Met_updated)) %>%
  mutate(BLAST_species = str_remove(BLAST_species_primer, "_[^_]+$")) %>%
  select(species, BLAST_species)

tax_Met_updated
  
################################################################################

### Get BLAST summaries for each consensus sequence

# 1. Get a list of all .txt files in the directory
file_list <- list.files(path = "BLAST_summaries", 
                        pattern = "\\.txt$", 
                        full.names = TRUE)

# 2. Read each file into a list of matrices
# Using an anonymous function to ensure consistent formatting
df_list <- lapply(file_list, function(x) {
  as.data.frame(read.table(x, header = FALSE, sep = "\t"))
})

# 3. Combine them all into one single df
# rbind "stacks" them vertically; use cbind if you want to join them side-by-side
BLAST_df <- do.call(rbind, df_list)

head(BLAST_df)

BLAST_df <- BLAST_df %>% 
  mutate(BLAST_species = paste(V16, V17, sep = "_"))

head(BLAST_df)

################################################################################

head(BLAST_df)

species_counts <- BLAST_df %>%
  count(BLAST_species, name = "number_of_consensus_sequences") %>%
  arrange(desc(number_of_consensus_sequences))

species_counts

tax_Met_updated <- left_join(tax_Met_updated, species_counts, by = "BLAST_species")

tax_Met_updated <- distinct(tax_Met_updated)

tax_Met_updated

# Great! Now I have a record of the number of consensus sequences attributed 
# to each unique BLAST hit

# There are some duplicate species (sequences which matched different species 
# but I have lumped together e.g. Hyperoplus_immaculatus and Ammodytes_tobianus
# which I just said was Ammodytidae spp.)

tax_Met_counts <- tax_Met_updated %>%
  select(species, number_of_consensus_sequences) %>%
  group_by(species) %>%
  summarise(number_of_consensus_sequences = sum(number_of_consensus_sequences, na.rm = TRUE)) %>%
  arrange(desc(number_of_consensus_sequences))

tax_Met_counts <- tax_Met_counts %>%
  mutate(primer = case_when(
    species == "Pollachius virens" ~ "Both",
    species == "Pomatoschistus flavescens" ~ "Both",
    species == "Scomber scombrus" ~ "Both",
    species == "Pholis gunnellus" ~ "Both",
    species == "Ammodytidae spp." ~ "Both",
    species == "Symphodus melops" ~ "Both",
    species == "Taurulus bubalis" ~ "Both",
    species == "Trisopterus minutus" ~ "Both",
    species == "Clupea harengus" ~ "MiFishU",
    species == "Pollachius pollachius" ~ "MiFishU",
    species == "Myoxocephalus scorpius" ~ "MiFishU",
    species == "Sprattus sprattus" ~ "Both",
    species == "Ctenolabrus rupestris" ~ "Both",
    species == "Labrus bergylta" ~ "MiFishU",
    species == "Merlangius merlangus" ~ "MiFishU",
    species == "Pomatoschistus minutus" ~ "Both",
    species == "Gadus morhua" ~ "Tele02",
    species == "Ciliata mustela" ~ "Both",
    species == "Conger conger" ~ "Tele02", 
    species == "Anguilla anguilla" ~ "Tele02",
    species == "Pomatoschistus pictus" ~ "Tele02",
    species == "Merlangius merlangus" ~ "MiFishU",
    species == "Pomatoschistus minutus" ~ "Both",
    species == "Spinachia spinachia" ~ "MiFishU",
    species == "Callionymus lyra" ~ "Tele02",
    species == "Parablennius gattorugine" ~ "Tele02", 
    TRUE ~ "Other" # This catches everything else
  ))

unique(tax_Met_counts$species)

tax_Met_counts

############################### Plot ###########################################

tax_Met_counts

library(ggtext)
library(dplyr)

tax_Met_counts <- tax_Met_counts %>%
  mutate(
    species_md = ifelse(
      species == "Ammodytidae spp.",
      species,
      paste0("*", species, "*")
    )
  )


consensus_counts <- ggplot(
  tax_Met_counts,
  aes(
    x = log10(number_of_consensus_sequences + 1),
    y = reorder(species, number_of_consensus_sequences),
    fill = primer
  )
) +
  labs(title = "", fill = "Primer(s)",
       x = "log10(Number of Consensus\nSequences + 1)", y = "Taxon",
       fill = "Primer(s)") +
  geom_col() +
  geom_text(
    aes(label = number_of_consensus_sequences),
    hjust = -0.1,
    size = 5
  ) +
  scale_fill_manual(
    values = c(
      Both = "#E0AC6F",
      MiFishU = "#F7EDAF",
      Tele02 = "#D17D7F"
    )
  ) +
  scale_x_continuous(
    limits = log10(c(0, 3500) + 1),
    breaks = log10(c(0, 1, 10, 100, 1000) + 1),
    labels = c("0", "1", "10", "100", "1000")
  ) +
  scale_y_discrete(
    labels = function(x) ifelse(
      x == "Ammodytidae spp.",
      x,
      paste0("*", x, "*")
    )
  ) +
  theme(
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_blank(),           # remove major grid
    panel.grid.minor = element_blank(),           # remove minor grid
    axis.line = element_line(color = "black"),    # draw axes
    axis.ticks = element_line(color = "black"),
    axis.text.y = element_markdown(size = 12.5),  # our markdown labels
    axis.text.x = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
    legend.key.size = unit(1, "cm"),
    legend.position = "top",
    plot.title = element_text(size = 15)
  )

consensus_counts

ggsave(filename = "plots/consensus_counts/consensus_counts.png", plot = consensus_counts,
       device = "png", dpi = 600, units = "mm", height = 250, width = 250)










