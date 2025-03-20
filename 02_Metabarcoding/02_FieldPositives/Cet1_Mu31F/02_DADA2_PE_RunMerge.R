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

# FALSE if there is only one run
multiple_runs <- TRUE

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

if (multiple_runs) {
  # Run 2
  seqtab2 <- read.csv("02_DADA2Out/Run2_Out/seqtab_nochim.csv", row.names = 1) %>% as.matrix
  asv2 <- read.csv("02_DADA2Out/Run2_Out/seq_only.csv", row.names = 1)
  track2 <- read.csv("02_DADA2Out/Run2_Out/track.csv", row.names = 1)
  ## Check data dimension
  dim(seqtab2); dim(asv2); dim(track2)
  ## Replace ASV table column names
  colnames(seqtab2) <- asv2$x
  
  # Run 3
  seqtab3 <- read.csv("02_DADA2Out/Run3_Out/seqtab_nochim.csv", row.names = 1) %>% as.matrix
  asv3 <- read.csv("02_DADA2Out/Run3_Out/seq_only.csv", row.names = 1)
  track3 <- read.csv("02_DADA2Out/Run3_Out/track.csv", row.names = 1)
  ## Check data dimension
  dim(seqtab3); dim(asv3); dim(track3)
  ## Replace ASV table column names
  colnames(seqtab3) <- asv3$x
}

# -------------------------------------------------------------- #
# Merge DADA2 tables
# -------------------------------------------------------------- #
if (multiple_runs) {
  # Merge DADA2 tables by dada2::mergeSequenceTables
  seqtab_nochim <- mergeSequenceTables(seqtab1, seqtab2, seqtab3,
                                       repeats = "sum", orderBy = "abundance")
  # Merge track data
  ## Check sample names
  track_vars <- unique(c(colnames(track1), colnames(track2), colnames(track3)))
  track_names <- sort(unique(c(rownames(track1), rownames(track2), rownames(track3))))
  if (all(rownames(seqtab_nochim) %in% track_names)) {
    track_names <- rownames(seqtab_nochim)  
  } else {
    stop("Check rownames!")
  }
  track <- matrix(0, nrow = length(track_names), ncol = length(track_vars)) %>% as.data.frame
  rownames(track) <- track_names
  colnames(track) <- track_vars
  ## Assign new values
  vars_tmp <- track_vars[track_vars != "prop.last.first."]
  track[track_names %in% rownames(track1), vars_tmp] <- track[track_names %in% rownames(track1), vars_tmp] + track1[,vars_tmp]
  track[track_names %in% rownames(track2), vars_tmp] <- track[track_names %in% rownames(track2), vars_tmp] + track2[,vars_tmp]
  track[track_names %in% rownames(track3), vars_tmp] <- track[track_names %in% rownames(track3), vars_tmp] + track3[,vars_tmp]
  track$prop.last.first. <- track$nonchim / track$input
  
  # Delete specific objects to "02_MergeDADA2.R"
  rm(track_names)
  rm(track_vars)
  rm(vars_tmp)
  rm(asv1); rm(track1); rm(seqtab1)
  rm(asv2); rm(track2); rm(seqtab2)
  rm(asv3); rm(track3); rm(seqtab3)
} else {
  # Rename objects
  seqtab_nochim <- seqtab1
  track <- track1
  rm(asv1); rm(track1); rm(seqtab1)
}

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

