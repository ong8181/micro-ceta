####
#### Summarize results
#### 2025.02.03, R4.3.0
####

# Load library and functions
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(cowplot); packageVersion("cowplot")
library(ggsci); packageVersion("ggsci")
theme_set(theme_cowplot())
source("funcs/F02_HelperFunctions.R") # Helper function for visualization

# Set random seeds (for reproduction)
ran.seed <- 1234
set.seed(ran.seed)
wdir <- basename(rstudioapi::getSourceEditorContext()$path)
(output_folder <- paste0(str_sub(wdir, end = -3), "Out")); rm(wdir)
dir.create(output_folder)


# <-----------------------------------------------------> #
#  Load data
# <-----------------------------------------------------> #
# Load sample data
sample_sheet <- read.csv("sampledata/SampleSheet.csv")
## OTU-based analysis
seqtab_data <- read.csv("03_OTUClusteringOut/otu_table.csv", row.names = 1)
tax_sheet <- read.delim("04_TaxaAssignmentOut/final_classigntax_otu")
tax_seq <- read.csv("03_OTUClusteringOut/otu_only.csv", row.names = 1)
seq_i <- 2
## ASV-based analysis
# seqtab_data <- read.csv("02_DADA2Out/seqtab_nochim.csv", row.names = 1)
# tax_sheet <- read.delim("04_TaxaAssignmentOut/TaxaAssignmentOut_ASV/final_classigntax_asv")
# tax_seq <- read.csv("02_DADA2Out/seq_only.csv", row.names = 1)
#seq_i <- 1

# Check structure
dim(seqtab_data)
dim(sample_sheet)
dim(tax_sheet); dim(tax_seq)


# <-----------------------------------------------------> #
#  Compile data
# <-----------------------------------------------------> #
# Check whether there is 0 sequences samples
(zero_sample <- sample_sheet$Sample_Name2[is.na(match(sample_sheet$Sample_Name2, rownames(seqtab_data)))])
if(length(zero_sample) > 0){
  # Generate dummy data frame
  seqtab_data_tmp <- matrix(0, ncol = ncol(seqtab_data), nrow = nrow(sample_sheet))
  rownames(seqtab_data_tmp) <- sample_sheet$Sample_Name2
  colnames(seqtab_data_tmp) <- colnames(seqtab_data)
  # Add object
  seqtab_data_tmp[match(rownames(seqtab_data), rownames(seqtab_data_tmp)),] <-
    as.matrix(seqtab_data)
  # Replace seq tables
  seqtab_data <- as.data.frame(seqtab_data_tmp)
}


# <-----------------------------------------------------> #
#  Import to phyloseq
# <-----------------------------------------------------> #
# Check names and orders
rownames(seqtab_data) == sample_sheet$Sample_Name2
sort(rownames(seqtab_data)) == sort(sample_sheet$Sample_Name2)
colnames(seqtab_data)
rownames(tax_seq); rownames(tax_sheet); tax_sheet$query
# Adjust row- and col-names
rownames(sample_sheet) <- rownames(seqtab_data) <- sample_sheet$Sample_Name2
tax_sheet$seq <- tax_seq[,seq_i]
tax_sheet$seqlen <- nchar(tax_seq[,seq_i])
rownames(tax_sheet) <- colnames(seqtab_data)

# Species assignment using assignSpecies()
## Load library
library(dada2); packageVersion("dada2") # v1.30.0, 2024.7.1
## Assign species
spp_assign <- assignSpecies(tax_sheet$seq, "../custom_db/marine_other_mammals.fasta.gz", allowMultiple = TRUE)
rownames(spp_assign) <- rownames(tax_sheet)
tax_sheet$assignSp_Genus <- spp_assign[,"Genus"]
tax_sheet$assignSp_Species <- spp_assign[,"Species"]
## Single species assignment or not
tax_sheet$assignSp_Single <- !is.na(spp_assign[,"Species"]) & !str_detect(spp_assign[,"Species"], "/")
sum(tax_sheet$assignSp_Single)

# Import to phyloseq
ps_all <- phyloseq(otu_table(seqtab_data, taxa_are_rows = FALSE),
                   sample_data(sample_sheet),
                   tax_table(as.matrix(tax_sheet)))

# Visualize
target_taxa <- "family"
ps_sub <- taxa_name_summarize(ps_all, target_taxa, top_taxa_n = 10)
ps_m1 <- speedyseq::psmelt(ps_sub)
ps_m2 <- ps_m1 %>% group_by_at(c("Sample", target_taxa)) %>%
  summarize(sequence_reads = sum(Abundance))
ps_m3 <- ps_m1 %>% group_by(Sample, rep_tax) %>%
  summarize(sequence_reads = sum(Abundance))

# Figures
f1 <- ggplot(ps_m2, aes_(x = as.name("Sample"), y = as.name("sequence_reads"),
                         fill = as.name(target_taxa))) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  scale_fill_igv() +
  xlab(NULL) + ylab("Sequence reads")

f2 <- ggplot(ps_m3, aes(x = Sample, y = sequence_reads, fill = rep_tax)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  xlab(NULL) + ylab("Sequence reads") +
  scale_fill_igv(name = target_taxa) +
  NULL

f3 <- plot_richness(ps_all, measures = "Observed") + xlab(NULL) + theme(axis.text.x = element_text(size = 6))

# Assign single species?
ps_m4 <- ps_m1 %>% group_by(Sample, assignSp_Single) %>%
  summarize(sequence_reads = sum(Abundance)) #%>% filter(sequence_reads > 0)
f4 <- ggplot(ps_m4, aes(x = Sample, y = sequence_reads, fill = assignSp_Single)) +
  geom_bar(stat = "identity", colour = NA) +
  theme(axis.text.x = element_text(angle = -90, hjust = 1, vjust = 0.5, size = 6)) + 
  xlab(NULL) + ylab("Sequence reads") +
  scale_fill_manual(values = c("gray80", "red3"), name = "Single species assigned?") +
  NULL


# --------------------------------------------- #
# Save results
# --------------------------------------------- #
# Output figures
pdf(file = sprintf("%s/Summary.pdf", output_folder), width = 18, height = 6)
plot_grid(f2, f3, ncol = 2, rel_widths = c(1,0.7))
dev.off()
# Single species assigned?
ggsave(sprintf("%s/SingleSp.pdf", output_folder), plot = f4, width = 10, height = 6)

# Re-output data
write.csv(otu_table(ps_all), sprintf("%s/otu_table.csv", output_folder))
write.csv(data.frame(sample_data(ps_all)), sprintf("%s/sample_data.csv", output_folder))
write.csv(data.frame(tax_table(ps_all)), sprintf("%s/tax_table.csv", output_folder))
saveRDS(ps_all, sprintf("%s/ps_all.obj", output_folder))

# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("00_SessionInfo/", output_folder, "_", substr(Sys.time(), 1, 10), ".txt"))
