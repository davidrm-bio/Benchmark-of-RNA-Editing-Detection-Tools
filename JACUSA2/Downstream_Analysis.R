# Working process

# Load Library
library(JACUSA2helper)
library(plyranges)
library(magrittr)
library(tidyverse)

# Script for JACUSA2 downstream analysis

args = commandArgs(trailingOnly=TRUE) 

# Argument 1 - Jacusa2 output file
# Argument 2 - Name for condition 1
# Argument 3 - Name for condition 2

# Load Jacusa2 output file

jacusa <- read_result(arg[1], showProgress = T)

# Filter Jacusa2 as follow:
# - Removing sites with score < 2
# - Retaining sites with coverage >= 10 for all replicates
# - Remove sites with > 2 observed bases (excluding reference base)
# - Apply a filter to retain only robust sites (RNA editing must be present in all replicates)

filtered <- jacusa

filtered$bc <-
  lapply_cond(filtered$bases, function(b) {
    Reduce("+", b)
  }) %>% Reduce("+", .) %>% base_count()

filtered <- filtered %>%
  dplyr::filter(score >= 2) %>%
  dplyr::filter(All(cov$cond1 >= 10) & All(cov$cond2 > 10)) %>%
  dplyr::filter(bc <= 2) %>%
  dplyr::filter(robust(bases))

# Combine Counts from Replicates
cond1 <- Reduce("+", filtered$bases$cond1) %>%  mutate (Pos= filtered@ranges@start) 
cond2 <- Reduce("+", filtered$bases$cond2) %>%  mutate (Pos= filtered@ranges@start) 

# Compare Condition against the Reference
cond1VsRef <- base_sub(cond1[1:3], filtered$ref)
cond2VsRef <- base_sub(cond2[1:3], filtered$ref)

# Frequency Table - Just counts
table (cond1VsRef) -> cond1Table
table (cond2VsRef) -> cond2Table

ResultTable <-  cond1Table - cond2Table

# Specifically looking into A-to-G changes

cond1df <-
  data.frame("Changes" = cond1VsRef, "Pos" = cond1$Pos) %>%
  filter (Changes == "A->G")

cond2df <-
  data.frame("Changes" = cond2VsRef, "Pos" = cond2$Pos) %>%
  filter (Changes == "A->G")

jacusa_A_to_G <- cond1df %>%
  filter (!Pos %in% cond2df$Pos)

jacusa_A_to_G_v2 <- cond2df %>%
  filter (!Pos %in% cond1df$Pos)

print ("A->G changes in Condition 1 not present in Condition 2")
length (unique (jacusa_A_to_G$Pos))

print ("A->G changes in Condition 2 not present in Condition 1")
length (unique (jacusa_A_to_G_v2$Pos))


# Load REDIportal database
REDIportal <-
  read.table(
    "/binf-isilon/rennie/gsn480/data/compare_tools/dataset1/dbRNA-Editing/REDIportal_db_GRCh38.txt",
    sep = "\t",
    header = FALSE,
    skip = 1
  )

REDIportal %>% 
  filter (V3 == "A" & V4 == "G") -> REDIportal_A_to_G


jacusa_A_to_G %>%
  filter (Pos %in% REDIportal_A_to_G$V2)
  

# 26 not present
# 1415 present

# 3 present 

# 604 present
cond2df %>% 
  filter (!Pos %in% REDIportal_A_to_G$V2) -> tmp

cond1df %>%
  filter (!Pos %in% tmp$Pos)
  #filter (Pos %in% REDIportal_A_to_G$V2)


# Load DARNED hg19 & Adapt format for Lift Genome Annotations
darned_tmp <-
  read.table(
    "/binf-isilon/rennie/gsn480/data/compare_tools/dataset1/dbRNA-Editing/DARNED_hg19.txt",
    sep = "\t",
    header = T
  ) %>%
  mutate (chrom = gsub(" ", "", paste("chr", chrom))) %>% mutate (coordinate.1 = coordinate) %>%
  relocate(coordinate.1, .after =
             coordinate) %>%
  select (1:3)  %>% 
  mutate (name= apply(., 1, function(x) gsub (" ", "", paste0 (x[1],"-", x[2])))) 

# Save to run in Lift Genome Annotations
write.table(
  darned.tmp,
  "/binf-isilon/rennie/gsn480/data/compare_tools/dataset1/dbRNA-Editing/DARNED_hg19.tmp",
  quote = F,
  col.names = F,
  sep = "\t",
  row.names = F
)

# Update Database with new positions
darned_original <-
  read.table(
    "/binf-isilon/rennie/gsn480/data/compare_tools/dataset1/dbRNA-Editing/DARNED_hg19.txt",
    sep = "\t",
    header = T
  )

darned_converted <-
  read.table(
    "/binf-isilon/rennie/gsn480/data/compare_tools/dataset1/dbRNA-Editing/DARNED_hg19_to_GRCh38.txt",
    sep = "\t",
    header = F
  ) %>% select(1:2, 4)


darned_original %>% 
  mutate (Name = apply(.,1, function(x) gsub(" ", "", paste("chr", x[1], "-", x[2])))) -> darned_original


darned_original %>%  select (3:11) %>% relocate(Name, .before=strand) -> darned.tmp

names(darned_converted) <- c("chrom", "coordinate", "Name")

merge (darned_converted, darned.tmp, by="Name") -> darned_new 

# 371,176
# 333,215

darned_new[duplicated(darned_new), ]


