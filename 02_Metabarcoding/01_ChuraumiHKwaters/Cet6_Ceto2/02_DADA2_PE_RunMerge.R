####
#### DADA2 analysis of Illumina fastq: Merge multiple runs
#### 2025.02.03, R4.3.0
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
output_folder <- "02_DADA2Out"
dir.create(output_folder)


# -------------------------------------------------------------- #
# Load DADA2 results from individual sequence runs
# -------------------------------------------------------------- #
# Run 1
seqtab1 <- read.csv("02_DADA2Out/Run1_Out/seqtab_nochim.csv", row.names = 1) %>% as.matrix
asv1 <- read.csv("02_DADA2Out/Run1_Out/seq_only.csv", row.names = 1)
track1 <- read.csv("02_DADA2Out/Run1_Out/track.csv", row.names = 1)
## Check data dimension
dim(seqtab1); dim(asv1); dim(track1)
## Replace ASV table column names
colnames(seqtab1) <- asv1$x


# -------------------------------------------------------------- #
# Merge DADA2 tables
# -------------------------------------------------------------- #
# Rename objects
seqtab_nochim <- seqtab1
track <- track1
rm(asv1); rm(track1); rm(seqtab1)


# -------------------------------------------------------------- #
# Save merged DADA2 tables
# -------------------------------------------------------------- #
# Taxa output for claident tax assignment
seqs <- colnames(seqtab_nochim)
seqs_out <- as.matrix(c(rbind(sprintf(">ASV%05d", 1:length(seqs)), seqs)), ncol = 1)
seqtab_nochim_for_csv <- seqtab_nochim
colnames(seqtab_nochim_for_csv) <- sprintf("ASV%05d", 1:length(seqs))

# Save outputs
write.csv(seqs, paste0(output_folder, "/seq_only.csv"), row.names = colnames(seqtab_nochim_for_csv))
write.table(seqs_out, paste0(output_folder, "/ASV_seqs.fa"), col.names = FALSE, row.names = FALSE, quote = FALSE)
write.csv(seqtab_nochim_for_csv, paste0(output_folder, "/seqtab_nochim.csv"))
write.csv(track, paste0(output_folder, "/track.csv"))

# Save workspace
save(list = ls(all.names = TRUE),
     file = paste0(output_folder, "/", output_folder, "_Merge.RData"))

# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("00_SessionInfo/", output_folder, "_Merge_", substr(Sys.time(), 1, 10), ".txt"))

