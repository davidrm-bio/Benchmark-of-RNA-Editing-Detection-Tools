# Library
library(JACUSA2helper)
library(tidyverse)
library(plyranges)
library(magrittr)

# Define variables

path = "path/to/file/file.out" # PATH to Jacusa Output

# Load data
jacusa <- read_result(path)

# Filter data
## Coverage >= 2
## Remove sites with more than 2 observed bases
## Retain only robust sites
## Remove Artefacts (D option of JACUSA)

jacusa$bc <- lapply_cond(jacusa$bases, function(b) { Reduce("+", b) } ) %>% Reduce("+", .) %>% base_count()
jacusa %>% 
  dplyr::filter(All(cov$cond1 >= 8) & All(cov$cond2 >= 8)) %>% 
  dplyr::filter (filter_artefact(filter, "D") == FALSE) %>% 
  dplyr::filter(bc <= 2) %>% 
  dplyr::filter(robust(bases)) -> jacusa.filt

# Create Table for python script
regions <- rep(jacusa.filt@seqnames@values, times =jacusa.filt@seqnames@lengths)

# Get WT changes
wt <- Reduce("+", jacusa.filt$bases$cond1)
ref2wt <- base_sub(wt, jacusa.filt$ref)
# Get ADAR changes
adarko <- Reduce("+", jacusa.filt$bases$cond2)
ref2adar <- base_sub(adarko, jacusa.filt$ref)
# Create df
jacusa.df <- data.frame(Region =regions, Position = jacusa.filt@ranges@start,  wt=ref2wt, adarko=ref2adar)

# Save to File
write.table(jacusa.df, "Jacusa_clean.tab", quote=F, sep="\t", row.names = F, col.names = T)