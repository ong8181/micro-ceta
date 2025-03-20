####
#### Coverage-based rarefaction
#### R 4.3.2
#### 2024.09.23 Ushio
####


# Load library and functions
library(tidyverse); packageVersion("tidyverse") # 2.0.0, 2024.09.23
library(phyloseq); packageVersion("phyloseq") # 1.46.0, 2024.09.23
library(cowplot); packageVersion("cowplot") # 1.1.3, 2024.09.23
theme_set(theme_cowplot())

# Load custom R package
# https://github.com/ong8181/macam
library(macam); packageVersion("macam") # 0.1.4

# Generate output folder
(output_folder <- macam::outdir_create())


# --------------------------------------------- #
# Load data
# --------------------------------------------- #
ps_all <- readRDS("05_SummarizeOut/ps_all.obj")
sample_sums(ps_all)
sample_data(ps_all)


# --------------------------------------------- #
# Extract non-NC samples & samples with > 0 seqs
# !! You may change depending on the purpose and sample characteristics
# --------------------------------------------- #
ps_sample <- subset_samples(ps_all, sample_nc == "sample") %>% 
  prune_samples(sample_sums(.) > 0, .) # Remove samples with no reads


# --------------------------------------------- #
# Coverage-based rarefaction for samples with sufficient reads
# --------------------------------------------- #
macam::coverage_info(ps_sample) # Check coverage info

## Coverage-based rarefaction
system.time(
  ps_rare_raw <- macamseq::rarefy_even_coverage(ps_sample,
                                                coverage = 0.999,
                                                knots = 100,
                                                n_rarefy_iter = 10,
                                                include_iNEXT_results = TRUE)
)
ps_rare <- ps_rare_raw[[1]] # Extract phyloseq object


# ----------------------------------------------- #
# Visualize rarefaction curve
# ----------------------------------------------- #
r1 <- plot_rarefy(ps_rare_raw, plot_rarefied_point = FALSE) +
  theme(legend.position = "none")
r2 <- plot_rarefy(ps_rare_raw, plot_rarefied_point = TRUE, plot_slope = TRUE) +
  theme(legend.position = "none")

# Remove ASVs with 0 reads and transform to the relative abundance
ps_rare <- ps_rare %>%
  prune_taxa(taxa_sums(.) > 0, .)
# Replace ASV less than 0.5% with 0
ps_rel <- ps_rare %>%
  transform_sample_counts(function(x) x/sum(x))
otu_table(ps_rel)[otu_table(ps_rel) < 0.005] <- 0
ps_rel <- ps_rel %>% prune_taxa(taxa_sums(.) > 0, .)

# --------------------------------------------- #
# Save results
# --------------------------------------------- #
# Save rarefied ps_object
saveRDS(ps_rare_raw, sprintf("%s/ps_rare_raw.obj", output_folder))
saveRDS(ps_rare, sprintf("%s/ps_rare.obj", output_folder))
saveRDS(ps_rel, sprintf("%s/ps_rel.obj", output_folder))
# Re-output data
write.csv(otu_table(ps_rel), sprintf("%s/otu_table_rel.csv", output_folder))
write.csv(data.frame(sample_data(ps_rel)), sprintf("%s/sample_data_rel.csv", output_folder))
write.csv(data.frame(tax_table(ps_rel)), sprintf("%s/tax_table_rel.csv", output_folder))

# Save plot
ggsave(sprintf("%s/Rarefaction_curve.pdf", output_folder), plot = r1, width = 10, height = 6)
ggsave(sprintf("%s/Rarefaction_curve_w_slope.pdf", output_folder), plot = r2, width = 10, height = 6)

# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("00_SessionInfo/", output_folder, "_", substr(Sys.time(), 1, 10), ".txt"))

save(list = ls(all.names = TRUE),
     file = paste0(output_folder, "/", output_folder, ".RData"))

