####
#### DADA2 analysis of Illumina fastq
#### Paired-end
#### 2024.04.04, R4.2.2
####

#--------------- DADA2 processing ---------------#
# Load library and functions
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(tidyverse); packageVersion("tidyverse")
source("funcs/F01_HelperFunctions.R")

# Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)
dir.create("00_SessionInfo")

# Generate output folder
output_folder <- "02_DADA2Out/Run2_Out"
dir.create("02_DADA2Out")
dir.create(output_folder)

# Load sequence reads
path <- paste0(getwd(), "/01_QualityFiltering_Run2Out")
system("rm 01_QualityFiltering_Run2Out/unknown*")
fnFs <- sort(list.files(path, pattern="R1.fastq.gz", full.names = T)) # Forward read files
fnRs <- sort(list.files(path, pattern="R2.fastq.gz", full.names = T)) # Reverse read files
# Get sample names
(sample_names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1))

# Visualize quality
#plotQualityProfile(fnFs[1:2])
#plotQualityProfile(fnRs[1:2])


# ------------------------------------ #
# Primer removal check
# ------------------------------------ #
# Mu31F - Dc320R
#GACACTGAAAATGTCTAGATGG
#TYAATCGTATGACCGCGGTG
# Dc671F - Dc1015R
#GCTACTYCAGTCTATATACC
#CACACYTTCCRGTAYGCTTACC
# Mu2084F - Dc2438R
#ATGAAYGGCCACACGAGGGTTTTA
#TGTCCTGATCCAACATCGAGG
# Mu9459F - Mu9822R
#CTGACTTCCAATCAGTTRGTTTCGG
#CATTCTARRCCYTYTTGRG

# Identify primers
FWD <- "ATGAAYGGCCACACGAGGGTTTTA"
REV <- "TGTCCTGATCCAACATCGAGG"
FWD_orients <- AllOrients(FWD)
REV_orients <- AllOrients(REV)

# Identify primers
seq_id <- 1
rbind(FWD.ForwardReads = sapply(FWD_orients, PrimerHits, fn = fnFs[[seq_id]]),
      FWD.ReverseReads = sapply(FWD_orients, PrimerHits, fn = fnRs[[seq_id]]), 
      REV.ForwardReads = sapply(REV_orients, PrimerHits, fn = fnFs[[seq_id]]), 
      REV.ReverseReads = sapply(REV_orients, PrimerHits, fn = fnRs[[seq_id]]))

# Performing filtering and trimming
filt_path <- file.path(path, "filtered") # Place filtered files in filtered/ subdirectory
filtFs <- file.path(filt_path, paste0(sample_names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample_names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                     truncLen=c(200, 200), # Remove low-quality bases at the end
                     maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = T,
                     compress = TRUE, multithread = TRUE) # On Windows set multithread=FALSE
#head(out)
out
#plotQualityProfile(filtFs[1:2])
#plotQualityProfile(filtRs[1:2])

# Exclude 0 seq samples, rename filtFs and filtRs
if(length(sample_names[out[,2]<1 | out[,1]<1]) > 0){
  filtFs <- file.path(filt_path, paste0(sample_names[out[,2]>0 & out[,1]>0], "_F_filt.fastq.gz"))
  filtRs <- file.path(filt_path, paste0(sample_names[out[,2]>0 & out[,1]>0], "_R_filt.fastq.gz"))
}

# Learn the error rates
min_nbases <- 200 * sum(out[,2]) # Use a small number of bases to speed up the analysis
errF <- learnErrors(filtFs, multithread=TRUE, randomize = TRUE, MAX_CONSIST = 20, nbases = min_nbases)
errR <- learnErrors(filtRs, multithread=TRUE, randomize = TRUE, MAX_CONSIST = 20, nbases = min_nbases)

# Visualize errors
#plotErrors(errF, nominalQ = T)
#plotErrors(errR, nominalQ = T)
ggsave(sprintf("%s/errF.pdf", output_folder), plotErrors(errF, nominalQ = T), width = 10, height = 10)
ggsave(sprintf("%s/errR.pdf", output_folder), plotErrors(errR, nominalQ = T), width = 10, height = 10)

# Dereplicatin
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample_names[out[,2]>0 & out[,1]>0]
names(derepRs) <- sample_names[out[,2]>0 & out[,1]>0]

# Sample inference
dadaFs <- dada(derepFs, err = errF, multithread = TRUE, pool = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE, pool = TRUE)
#dadaFs[[1]]

# Merging paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE, trimOverhang = TRUE)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab); sum(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
seqtab2 <- seqtab
table(nchar(getSequences(seqtab2)))

# Remove chimeras
seqtab_nochim <- removeBimeraDenovo(seqtab2, method = "consensus", multithread = TRUE, verbose = TRUE)
table(nchar(getSequences(seqtab_nochim)))
dim(seqtab_nochim)
sum(seqtab_nochim)/sum(seqtab2)

# Track reads thourhg the pipeline
out2 <- out[out[,2]>0 & out[,1]>0,]
getN <- function(x) sum(getUniques(x))
track <- cbind(out2, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab2), rowSums(seqtab_nochim),  rowSums(seqtab_nochim)/out2[,1])
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "tabled2", "nonchim", "prop(last/first)")
rownames(track) <- sample_names[out[,2]>0 & out[,1]>0]
head(track)

# Taxa output for claident tax assignment
seqs <- colnames(seqtab_nochim)
seqs_out <- as.matrix(c(rbind(sprintf(">ASV%05d", 1:length(seqs)), seqs)), ncol = 1)
seqtab_nochim_for_csv <- seqtab_nochim
colnames(seqtab_nochim_for_csv) <- sprintf("ASV%05d", 1:length(seqs))

# Save outputs
write.csv(seqs, paste0(output_folder, "/seq_only.csv"), row.names = colnames(seqtab_nochim))
write.table(seqs_out, paste0(output_folder, "/ASV_seqs.fa"), col.names = FALSE, row.names = FALSE, quote = FALSE)
write.csv(seqtab_nochim_for_csv, paste0(output_folder, "/seqtab_nochim.csv"))
write.csv(track, paste0(output_folder, "/track.csv"))

# Save workspace
rm(derepFs)
rm(derepRs)
rm(dadaFs)
rm(dadaRs)
save(list = ls(all.names = TRUE),
     file = paste0(output_folder, "/DADA2Out.RData"))

# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("00_SessionInfo/02_DADA2Out_Run1_", substr(Sys.time(), 1, 10), ".txt"))


