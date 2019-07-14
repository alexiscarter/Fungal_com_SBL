## Create phyloseq object

#### Packages ####
library(dplyr); packageVersion("dplyr") # ‘0.8.0.1’
library(vegan); packageVersion("vegan") # ‘2.5.4’
library(phyloseq);packageVersion("phyloseq") # '1.26.1'

#### Data import and processing ####

# SBL data
label <- read.csv("data/label_soil_sample.csv", sep = ";")
load("data/soil3.rda")

# Taxonomy table

load("dada2/saved_table/taxa.ITS.paired.13.05.dynamic_all_02.02.2019.minboot50.rda")
taxa_table <- taxa.paired

# ASV table
load("dada2/saved_table/seqtab.nochim.ITS.paired.10.05.rdata")
asv_table <- seqtab.nochim

# Discard soil samples
# S40 because M_AS_FS_03_Ah06 replicate of M_AS_FS_03_Ah
# all the Ah horizons: S70, S71, S57, S32, S38
# S19 as outlier
remove_sample <- c("S19", "S40", "S71", "S70", "S57", "S32", "S38")
asv_table <- asv_table[!rownames(asv_table) %in% remove_sample, ]

# Rename
soil3$sampleID <- gsub("1/2", "Mix", soil3$sampleID)

# Make a data.frame holding the sample data
soil <- soil3 %>%
  left_join(label, by = 'sampleID')

soil <- soil[!soil$sample %in% remove_sample, ]
soil$sample <- droplevels(soil$sample)

rownames(soil) <- soil$sample

# Standardize environmental data
chem <- c("totalC", "totalN", "CN", "totalP", "inorgP", "orgP", "BrayP", "CPo", "pHCaCl2", "Al", "Ca", "Fe", "K",  "Mg", "Mn", "Na", "TEB", "ECEC", "BS")
soil[,chem] <- decostand(soil[,chem], method = "standardize")

# Clean unuseful columns
remove_cols <- c("plot.x", "block.x", "suitable", "width", "succ", "watregime", "observer", "date_veg", "rawplot")
soil <- soil[,!colnames(soil) %in% remove_cols]

# Construct phyloseq object (straightforward from dada2 outputs)
ps <- phyloseq(otu_table(asv_table, taxa_are_rows=FALSE), 
               sample_data(soil), 
               tax_table(as.matrix(taxa_table)))

# Reorder factors
sample_data(ps)$horiz = factor(sample_data(ps)$horiz, levels = c('L', 'F', 'H', 'Ae', 'B'))
sample_data(ps)$myco = factor(sample_data(ps)$myco, levels = c('AM', 'mixed', 'ECM'))

## save ps object
# save(ps, file = "data/ps.paired.ITS.29.05.minboot50.rda")
