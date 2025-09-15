# 25.06.25
# Ethan Ross

## In this script we will get the seaDNA database for the UK

## See https://github.com/genner-lab/meta-fish-li

# load packages (install if required)
library("tidyverse")
library("ape")

# load REMOTE references and cleaning scripts (requires internet connection)
source("https://raw.githubusercontent.com/genner-lab/meta-fish-lib/main/scripts/references-load-remote.R")
source("https://raw.githubusercontent.com/genner-lab/meta-fish-lib/main/scripts/references-clean.R")

# subset reference library table by metabarcode fragment (primer set) from the following options:
print(tibble(metabarcodes=c("coi.lerayxt","coi.ward","12s.miya","12s.riaz","12s.valentini","12s.taberlet","16s.berry","cytb.minamoto")))
# change 'metabarcode' argument as appropriate:

# Looks like these are identicle!
# reflib.sub.MiFishU <- subset_references(df=reflib.cleaned, metabarcode="12s.miya")
# reflib.sub.tele02 <- subset_references(df=reflib.cleaned, metabarcode="12s.taberlet")

reflib.sub.12S <- subset_references(df=reflib.cleaned, metabarcode="12s.miya")

# [OPTIONAL] taxonomically dereplicate and filter on sequence length
# 'proplen=0.5' removes sequences shorter than 50% of median sequence length
# 'proplen=0' retains all sequences
reflib.sub.12S <- derep_filter(df=reflib.sub.12S, derep=TRUE, proplen=0.5)

# write out reference library in various formats to current working directory
# currently supported formats are: [1] sintax, [2] dada2 (taxonomy), [3] dada2 (species), [4] Qiime2 (fasta and taxonomy), and [5] plain dbid (GenBank or BOLD database identifiers)
write_references_fasta(df=reflib.sub.12S)

