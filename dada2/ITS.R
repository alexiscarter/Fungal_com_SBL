## DADA2 workflowe for SBL Illumina sequences with ITS as marker gene
## Info
# Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., & Holmes, S. P. (2016). DADA2: high-resolution sample inference from Illumina amplicon data. Nature methods, 13(7), 581-583.
# http://benjjneb.github.io/dada2/tutorial.html

#### Load packages and data
library(dada2); packageVersion("dada2") # ‘1.11.1’ Feb-12-2019
path <- "/home/udem/Documents/dada2/ITS"
list.files(path)

#### Paths ####
# Read in the names of the fastq files, and perform some string manipulation to get lists of the forward and reverse fastq files in matched order:
# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(path, pattern="_R1.fastq"))
fnRs <- sort(list.files(path, pattern="_R2.fastq"))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(fnFs, "_"), `[`, 1)
sample.names 

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

#### Perform filtering and trimming ####
filt_path <- file.path(path, "filtered_pairedend") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Filter the forward and reverse reads:
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, truncQ=6, rm.phix=TRUE, trimLeft=c(18,20), compress=TRUE, multithread=TRUE,
                     truncLen = c(280,280), maxEE=c(3,3))

mean(out[,2])/mean(out[,1])

#### Learn the Error Rates ####
# It learns the error model from the data, by alternating estimation of the error rates and inference of sample
errF <- learnErrors(filtFs, nbases = 1e+10, multithread=TRUE)
errR <- learnErrors(filtRs, nbases = 1e+10, multithread=TRUE)

#### Dereplication ####
#combines all identical sequencing reads into into “unique sequences” with a corresponding “abundance”
#reduces computation time by eliminating redundant comparisons
derepFs <- derepFastq(filtFs,  n = 1e+06)
derepRs <- derepFastq(filtRs, n = 1e+06)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#### Sample Inference ####
# Infer the sequence variants in each sample
dadaFs <- dada(derepFs,  err=errF, BAND_SIZE = 16, KDIST_CUTOFF = 0.42, pool=TRUE, multithread=TRUE)

dadaRs <- dada(derepRs, err=errR,pool=TRUE, multithread=TRUE)

### Merging ####
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, maxMismatch = 0, returnRejects = FALSE, propagateCol = character(0), justConcatenate = FALSE, trimOverhang = FALSE,
                      minOverlap = 12)

#### Construct sequence table ####
# SV version of the OTU table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Inspect distribution of sequence lengths
seqtab.length.var.plot <- plot(table(nchar(getSequences(seqtab))))

#### Remove chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="pooled", multithread=TRUE, verbose=TRUE) 
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) #  percentage of the total sequence reads

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab.nochim)))

#### Assign taxonomy ####
taxa.paired <- assignTaxonomy(seqtab.nochim, "reference_database/sh_general_release_dynamic_all_02.02.2019.fasta", minBoot = 50, multithread=TRUE, verbose=TRUE)

#### Track reads through the pipeline ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, 
               sapply(dadaFs, getN), 
               sapply(dadaRs, getN), 
               sapply(mergers, getN), 
               rowSums(seqtab.nochim))

colnames(track) <- c("input", 
                     "filtered", 
                     "denoised F", 
                     "denoised R", 
                     "merged", 
                     "nonchim paired")
rownames(track) <- sample.names

