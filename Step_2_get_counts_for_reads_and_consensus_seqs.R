# 23.06.25
# Ethan Ross

# In this script I will get the read counts for my metaprobe samples
# I have created a .csv file from the fastq files on the cluster

library(writexl)
library(readxl)
library(dplyr)
library(stringr)
library(tidyr)

met_reads <- read.csv("fastq_info/cluster_fastq_read_counts.csv", header = T)
head(met_reads)

### Add columns 

met_reads <- met_reads %>%
  mutate(
    # Temporarily split the Filename into parts
    parts = str_split(Filename, "_"),
    
    # Extract 'metaprobe':
    # 1. Take the 5th element (e.g., "Met")
    # 2. Take the 6th element (e.g., "12a", "9b", "RUN1")
    # 3. Remove 'a', 'b', or 'c' from the END of the 6th element using regex '[abc]$'
    # 4. Join the (modified) 5th and 6th elements with '_'
    metaprobe = sapply(parts, function(x) {
      cleaned_sixth_part <- str_remove(x[6], "[abc]$")
      paste(x[5], cleaned_sixth_part, sep = "_")
    }),
    
    # Extract 'marker': 7th element, with '.fastq' removed
    marker = sapply(parts, function(x) str_remove(x[7], "\\.fastq$"))
  ) %>%
  # Remove the temporary 'parts' column
  select(-parts)

met_reads

# --- 2. Create the 'run' column using dplyr's mutate and case_when ---
met_reads <- met_reads %>%
  mutate(
    # Temporarily extract the numeric part from 'metaprobe'
    # This will be NA for 'metaprobe' values like "PCR_C", "Field_C", etc.
    numeric_metaprobe_val = str_extract(metaprobe, "\\d+$"),
    
    # Use case_when for conditional assignment to the 'run' column
    run = case_when(
      # Numeric ranges based on the extracted number
      !is.na(numeric_metaprobe_val) & as.numeric(numeric_metaprobe_val) >= 1 & as.numeric(numeric_metaprobe_val) <= 3 ~ "RUN1",
      !is.na(numeric_metaprobe_val) & as.numeric(numeric_metaprobe_val) >= 4 & as.numeric(numeric_metaprobe_val) <= 6 ~ "RUN2",
      !is.na(numeric_metaprobe_val) & as.numeric(numeric_metaprobe_val) >= 7 & as.numeric(numeric_metaprobe_val) <= 9 ~ "RUN3",
      !is.na(numeric_metaprobe_val) & as.numeric(numeric_metaprobe_val) >= 10 & as.numeric(numeric_metaprobe_val) <= 12 ~ "RUN4",
      
      # Default: If none of the above conditions are met (e.g., for "PCR_C", "Field_C", etc.)
      TRUE ~ NA_character_ # Assign NA
    )
  ) %>%
  # Remove the temporary helper column for cleanliness
  select(-numeric_metaprobe_val)

met_reads$metaprobe

met_reads <- met_reads %>%
  mutate(run = if_else(marker == "RUN1", "RUN1", run))

met_reads

total_reads_per_run <- met_reads %>%
  filter(!marker == "RUN1") %>%
  group_by(run, marker) %>%
  dplyr::summarise(
    Total_Reads = sum(Number_of_Reads, na.rm = TRUE)
  ) %>%
  ungroup()

write_xlsx(total_reads_per_run, "fastq_info/reads_per_run.xlsx")

################################################################################

library(Biostrings)

### Okay now I will get the total number of ZOTUs

### First I need to know which sequneces belong to fishes 

RUN1_MiFishU_sum <- read.delim("BLAST_summaries/RUN1_MiFishU_BLAST_summary.txt",
                               sep = "\t", header = F)

RUN1_MiFishU_fish_con_seqs <- RUN1_MiFishU_sum %>%
  filter(!is.na(V2)) %>%
  mutate(
    V20_cleaned = str_replace(V20, "^consensus_", ""),
    # Then combine V1 and the cleaned V20 with an underscore
    seq_id = str_c(V1, V20_cleaned, sep = "_")
  ) %>%
  select(seq_id)

head(RUN1_MiFishU_fish_con_seqs$seq_id)

RUN1_MiFishU_fasta <- readDNAStringSet("final_fasta/RUN1_MiFishU_All_barcodes_combined.fasta")

#Get the names from your DNAStringSet
RUN1_MiFishU_fasta_names <- names(RUN1_MiFishU_fasta)

#Get the list of IDs you want to keep
RUN1_MiFishU_fasta_to_keep <- RUN1_MiFishU_fish_con_seqs$seq_id

#Filter the DNAStringSet
# Use the %in% operator to check which fasta names are present in ids_to_keep
RUN1_MiFishU_fish_fasta <- RUN1_MiFishU_fasta[RUN1_MiFishU_fasta_names %in% RUN1_MiFishU_fasta_to_keep]

RUN1_MiFishU_fish_unique <- unique(RUN1_MiFishU_fish_fasta)

length(RUN1_MiFishU_fasta)
length(RUN1_MiFishU_fish_unique)

######################################## RUN1 ##################################

RUN1_Tele02_sum <- read.delim("BLAST_summaries/RUN1_Tele02_BLAST_summary.txt",
                              sep = "\t", header = F)

RUN1_Tele02_sum <- RUN1_Tele02_sum %>%
  filter(!is.na(V2))

RUN1_Tele02_fish_con_seqs <- RUN1_Tele02_sum %>%
  filter(!is.na(V2)) %>%
  mutate(
    V20_cleaned = str_replace(V20, "^consensus_", ""),
    # Then combine V1 and the cleaned V20 with an underscore
    seq_id = str_c(V1, V20_cleaned, sep = "_")
  ) %>%
  select(seq_id)

RUN1_Tele02_fasta <- readDNAStringSet("final_fasta/RUN1_Tele02_All_barcodes_combined.fasta")

#Get the names from your DNAStringSet
RUN1_Tele02_fasta_names <- names(RUN1_Tele02_fasta)

#Get the list of IDs you want to keep
RUN1_Tele02_fasta_to_keep <- RUN1_Tele02_fish_con_seqs$seq_id

#Filter the DNAStringSet
# Use the %in% operator to check which fasta names are present in ids_to_keep
RUN1_Tele02_fish_fasta <- RUN1_Tele02_fasta[RUN1_Tele02_fasta_names %in% RUN1_Tele02_fasta_to_keep]

RUN1_Tele02_fish_unique <- unique(RUN1_Tele02_fish_fasta)

length(RUN1_Tele02_fasta)
length(RUN1_Tele02_fish_unique)

#################################### RUN2 ######################################

RUN2_MiFishU_sum <- read.delim("BLAST_summaries/RUN2_MiFishU_BLAST_summary.txt",
                               sep = "\t", header = F)

RUN2_MiFishU_fish_con_seqs <- RUN2_MiFishU_sum %>%
  filter(!is.na(V2)) %>%
  mutate(
    V20_cleaned = str_replace(V20, "^consensus_", ""),
    # Then combine V1 and the cleaned V20 with an underscore
    seq_id = str_c(V1, V20_cleaned, sep = "_")
  ) %>%
  select(seq_id)

head(RUN2_MiFishU_fish_con_seqs$seq_id)

RUN2_MiFishU_fasta <- readDNAStringSet("final_fasta/RUN2_MiFishU_All_barcodes_combined.fasta")

#Get the names from your DNAStringSet
RUN2_MiFishU_fasta_names <- names(RUN2_MiFishU_fasta)

#Get the list of IDs you want to keep
RUN2_MiFishU_fasta_to_keep <- RUN2_MiFishU_fish_con_seqs$seq_id

#Filter the DNAStringSet
# Use the %in% operator to check which fasta names are present in ids_to_keep
RUN2_MiFishU_fish_fasta <- RUN2_MiFishU_fasta[RUN2_MiFishU_fasta_names %in% RUN2_MiFishU_fasta_to_keep]

RUN2_MiFishU_fish_unique <- unique(RUN2_MiFishU_fish_fasta)

length(RUN2_MiFishU_fasta)
length(RUN2_MiFishU_fish_unique)

################################################################################

RUN2_Tele02_sum <- read.delim("BLAST_summaries/RUN2_Tele02_BLAST_summary.txt",
                              sep = "\t", header = F)

RUN2_Tele02_sum <- RUN2_Tele02_sum %>%
  filter(!is.na(V2))

RUN2_Tele02_fish_con_seqs <- RUN2_Tele02_sum %>%
  filter(!is.na(V2)) %>%
  mutate(
    V20_cleaned = str_replace(V20, "^consensus_", ""),
    # Then combine V1 and the cleaned V20 with an underscore
    seq_id = str_c(V1, V20_cleaned, sep = "_")
  ) %>%
  select(seq_id)

RUN2_Tele02_fasta <- readDNAStringSet("final_fasta/RUN2_Tele02_All_barcodes_combined.fasta")

#Get the names from your DNAStringSet
RUN2_Tele02_fasta_names <- names(RUN2_Tele02_fasta)

#Get the list of IDs you want to keep
RUN2_Tele02_fasta_to_keep <- RUN2_Tele02_fish_con_seqs$seq_id

#Filter the DNAStringSet
# Use the %in% operator to check which fasta names are present in ids_to_keep
RUN2_Tele02_fish_fasta <- RUN2_Tele02_fasta[RUN2_Tele02_fasta_names %in% RUN2_Tele02_fasta_to_keep]

RUN2_Tele02_fish_unique <- unique(RUN2_Tele02_fish_fasta)

length(RUN2_Tele02_fasta)
length(RUN2_Tele02_fish_unique)

#################################### RUN3 ######################################

RUN3_MiFishU_sum <- read.delim("BLAST_summaries/RUN3_MiFishU_BLAST_summary.txt",
                               sep = "\t", header = F)

RUN3_MiFishU_fish_con_seqs <- RUN3_MiFishU_sum %>%
  filter(!is.na(V2)) %>%
  mutate(
    V20_cleaned = str_replace(V20, "^consensus_", ""),
    # Then combine V1 and the cleaned V20 with an underscore
    seq_id = str_c(V1, V20_cleaned, sep = "_")
  ) %>%
  select(seq_id)

head(RUN3_MiFishU_fish_con_seqs$seq_id)

RUN3_MiFishU_fasta <- readDNAStringSet("final_fasta/RUN3_MiFishU_All_barcodes_combined.fasta")

#Get the names from your DNAStringSet
RUN3_MiFishU_fasta_names <- names(RUN3_MiFishU_fasta)

#Get the list of IDs you want to keep
RUN3_MiFishU_fasta_to_keep <- RUN3_MiFishU_fish_con_seqs$seq_id

#Filter the DNAStringSet
# Use the %in% operator to check which fasta names are present in ids_to_keep
RUN3_MiFishU_fish_fasta <- RUN3_MiFishU_fasta[RUN3_MiFishU_fasta_names %in% RUN3_MiFishU_fasta_to_keep]

RUN3_MiFishU_fish_unique <- unique(RUN3_MiFishU_fish_fasta)

length(RUN3_MiFishU_fasta)
length(RUN3_MiFishU_fish_unique)

################################################################################

RUN3_Tele02_sum <- read.delim("BLAST_summaries/RUN3_Tele02_BLAST_summary.txt",
                              sep = "\t", header = F)

RUN3_Tele02_sum <- RUN3_Tele02_sum %>%
  filter(!is.na(V2))

RUN3_Tele02_fish_con_seqs <- RUN3_Tele02_sum %>%
  filter(!is.na(V2)) %>%
  mutate(
    V20_cleaned = str_replace(V20, "^consensus_", ""),
    # Then combine V1 and the cleaned V20 with an underscore
    seq_id = str_c(V1, V20_cleaned, sep = "_")
  ) %>%
  select(seq_id)

RUN3_Tele02_fasta <- readDNAStringSet("final_fasta/RUN3_Tele02_All_barcodes_combined.fasta")

#Get the names from your DNAStringSet
RUN3_Tele02_fasta_names <- names(RUN3_Tele02_fasta)

#Get the list of IDs you want to keep
RUN3_Tele02_fasta_to_keep <- RUN3_Tele02_fish_con_seqs$seq_id

#Filter the DNAStringSet
# Use the %in% operator to check which fasta names are present in ids_to_keep
RUN3_Tele02_fish_fasta <- RUN3_Tele02_fasta[RUN3_Tele02_fasta_names %in% RUN3_Tele02_fasta_to_keep]

RUN3_Tele02_fish_unique <- unique(RUN3_Tele02_fish_fasta)

length(RUN3_Tele02_fasta)
length(RUN3_Tele02_fish_unique)

#################################### RUN4 ######################################

RUN4_MiFishU_sum <- read.delim("BLAST_summaries/RUN4_MiFishU_BLAST_summary.txt",
                               sep = "\t", header = F)

RUN4_MiFishU_fish_con_seqs <- RUN4_MiFishU_sum %>%
  filter(!is.na(V2)) %>%
  mutate(
    V20_cleaned = str_replace(V20, "^consensus_", ""),
    # Then combine V1 and the cleaned V20 with an underscore
    seq_id = str_c(V1, V20_cleaned, sep = "_")
  ) %>%
  select(seq_id)

head(RUN4_MiFishU_fish_con_seqs$seq_id)

RUN4_MiFishU_fasta <- readDNAStringSet("final_fasta/RUN4_MiFishU_All_barcodes_combined.fasta")

#Get the names from your DNAStringSet
RUN4_MiFishU_fasta_names <- names(RUN4_MiFishU_fasta)

#Get the list of IDs you want to keep
RUN4_MiFishU_fasta_to_keep <- RUN4_MiFishU_fish_con_seqs$seq_id

#Filter the DNAStringSet
# Use the %in% operator to check which fasta names are present in ids_to_keep
RUN4_MiFishU_fish_fasta <- RUN4_MiFishU_fasta[RUN4_MiFishU_fasta_names %in% RUN4_MiFishU_fasta_to_keep]

RUN4_MiFishU_fish_unique <- unique(RUN4_MiFishU_fish_fasta)

length(RUN4_MiFishU_fasta)
length(RUN4_MiFishU_fish_unique)

################################################################################

RUN4_Tele02_sum <- read.delim("BLAST_summaries/RUN4_Tele02_BLAST_summary.txt",
                              sep = "\t", header = F)

RUN4_Tele02_sum <- RUN4_Tele02_sum %>%
  filter(!is.na(V2))

RUN4_Tele02_fish_con_seqs <- RUN4_Tele02_sum %>%
  filter(!is.na(V2)) %>%
  mutate(
    V20_cleaned = str_replace(V20, "^consensus_", ""),
    # Then combine V1 and the cleaned V20 with an underscore
    seq_id = str_c(V1, V20_cleaned, sep = "_")
  ) %>%
  select(seq_id)

RUN4_Tele02_fasta <- readDNAStringSet("final_fasta/RUN4_Tele02_All_barcodes_combined.fasta")

#Get the names from your DNAStringSet
RUN4_Tele02_fasta_names <- names(RUN4_Tele02_fasta)

#Get the list of IDs you want to keep
RUN4_Tele02_fasta_to_keep <- RUN4_Tele02_fish_con_seqs$seq_id

#Filter the DNAStringSet
# Use the %in% operator to check which fasta names are present in ids_to_keep
RUN4_Tele02_fish_fasta <- RUN4_Tele02_fasta[RUN4_Tele02_fasta_names %in% RUN4_Tele02_fasta_to_keep]

RUN4_Tele02_fish_unique <- unique(RUN4_Tele02_fish_fasta)

length(RUN4_Tele02_fasta)
length(RUN4_Tele02_fish_unique)

################################################################################

unique_seq_list <- list(
  RUN1_MiFishU = RUN1_MiFishU_fasta,
  RUN2_MiFishU = RUN2_MiFishU_fasta,
  RUN3_MiFishU = RUN3_MiFishU_fasta,
  RUN4_MiFishU = RUN4_MiFishU_fasta,
  RUN1_Tele02 = RUN1_Tele02_fasta,
  RUN2_Tele02 = RUN2_Tele02_fasta,
  RUN3_Tele02 = RUN3_Tele02_fasta,
  RUN4_Tele02 = RUN4_Tele02_fasta
)

unique_seq_list_fish <- list(
  RUN1_MiFishU = RUN1_MiFishU_fish_unique,
  RUN2_MiFishU = RUN2_MiFishU_fish_unique,
  RUN3_MiFishU = RUN3_MiFishU_fish_unique,
  RUN4_MiFishU = RUN4_MiFishU_fish_unique,
  RUN1_Tele02 = RUN1_Tele02_fish_unique,
  RUN2_Tele02 = RUN2_Tele02_fish_unique,
  RUN3_Tele02 = RUN3_Tele02_fish_unique,
  RUN4_Tele02 = RUN4_Tele02_fish_unique
)

# Process the list to create the summary table
unique_sequence_counts_table <- tibble(
  Run_Marker = names(unique_seq_list), # Gets the names from the list (e.g., "RUN1_MiFishU")
  Unique_Sequence_Count = sapply(unique_seq_list, length),
  Unique_Sequence_Count_Fish = sapply(unique_seq_list_fish, length) # Gets the length (number of sequences) of each DNAStringSet
) %>%
  # Optionally, you might want to separate "Run" and "Marker" into two columns
  separate(Run_Marker, into = c("Run", "Marker"), sep = "_", remove = FALSE)

# Print the resulting table
print(unique_sequence_counts_table)

write_xlsx(unique_sequence_counts_table, "fastq_info/fish_con_seqs_per_run.xlsx")

