####
#### Evaluate interspecific variations for marine mammals
#### 2023.07.20 Ushio
####

# Load libraries
library(tidyverse); packageVersion("tidyverse") # 1.3.2, 2023.7.18
library(cowplot); packageVersion("cowplot") # 1.1.1, 2023.7.18
library(Biostrings); packageVersion("Biostrings") # 2.66.0, 2023.7.18
library(Rtsne); packageVersion("Rtsne") # 0.16, 2023.7.18
#library(macam); packageVersion("macam") # 0.1.4, 2023.7.18
theme_set(theme_cowplot())

# Create output directory
set.seed(1234)
(output_folder <- macam::outdir_create())

# Load sequence metadata
acc_taxdb <- read.csv("02_CompileTaxaOut/acc_tax_db_ed.csv")
# Modify the primer order
primer_order <- c("MiMammal", "MarVer1", "MarVer2", "MarVer3", "Ceto2", "Riaz12S",
                  "Mu31F_Dc320R", "Mu31F_MiMammalR", "Dc321F_Dc495R", "Dc494F_Mu643R",
                  "Mu642F_Dc1015R", "Mu642F_Mu1021R", "Dc671F_Dc1015R", "Dc671F_Mu1021R",
                  "Dc1458_1638", "Dc1465F_Dc1646R", "Mu2084F_Dc2438R", "Dc2173_2580",
                  "Mu2187F_Dc2438R", "Mu2187F_Mu2563R", "Dc2385_2580", "Dc2430_2580",
                  "Dc2438F_Mu2563R", "Dc2505_2580", "Mu9459F_Mu9822R", "Mu9459F_Mu9834R")


# ------------------------------------------------ #
# Load primers
# ------------------------------------------------ #
# Load tested primer names
primer_tested <- list.files("01_usearchOut/1_EachPrimerSet/")
primer_tested <- primer_tested[!primer_tested %in% c("CheckEditDistance_Cetacea", "CheckEditDistance_MarineMammal")]
# Drop some primer candidates
## There are dropped due to unexpected amplifications
primer_tested <- primer_tested[primer_tested != "Mu2743F_Mu2950R" & primer_tested != "Mu15826F_Mu15796R"]
## Which primers should we test?
primer_tested <- "MiMammal"

# Load compiled results
match_seq_detail <- readRDS("03_SummarizeInSilicoPCROut/match_seq_detail.obj")

# Set the number of samples to be analyzed
N_SAMPLE <- 5000


# ------------------------------------------------ #
# Dimension reduction
# ------------------------------------------------ #
#for (primer_i in primer_tested) {
for (primer_i in c("Mu31F_Dc320R", "Mu2084F_Dc2438R", "Mu2187F_Dc2438R", "Mu9459F_Mu9822R")) {
  # Record start time
  time_start <- proc.time()[3]
  # Load sequences
  primer <- readDNAStringSet(sprintf("data_primerset/each_primerset/%s.fa", primer_i))
  amplicon <- readDNAStringSet(sprintf("01_usearchOut/1_EachPrimerSet/%s/hits_all.fa", primer_i))
  match_df <- match_seq_detail[[primer_i]] %>%
    filter(forward_primer == names(primer)[1] & reverse_primer == names(primer[2])) %>% 
    filter(strand_fprimer == "+" & strand_rprimer == "-")
  
  # Downsample for a trial
  set.seed(1234)
  rand_id <- runif(N_SAMPLE, min = 1, max = nrow(match_df))
  #rand_acc <- match_df$accession_id
  rand_acc <- match_df$accession_id[rand_id]
  
  # Extract meta data
  amplicon <- amplicon[match(rand_acc, names(amplicon) %>% str_split(" ") %>% sapply('[', 1))]
  md_sub <- acc_taxdb[match(rand_acc, acc_taxdb$accession_id) %>% sort,]
  # Sort sequence so that it matches metadata
  amplicon <- amplicon[match(md_sub$accession_id, names(amplicon) %>% str_split(" ") %>% sapply('[', 1))]
  
  # Remove primers
  n_fp <- width(primer)[1]
  n_rp <- width(primer)[2]
  amplicon <- subseq(amplicon, start = n_fp+1, end = -n_rp-1)
  writeXStringSet(amplicon, sprintf("%s/amp_tmp.fa", output_folder))
  
  # Alignment by MAFFT
  system(sprintf("mafft %s/amp_tmp.fa > %s/amp_aligned.fa", output_folder, output_folder))
  amplicon_aln <- readDNAStringSet(sprintf("%s/amp_aligned.fa", output_folder))
  ## Change sequence names
  names(amplicon_aln) <- md_sub$species
  system(sprintf("rm %s/amp_tmp.fa %s/amp_aligned.fa", output_folder, output_folder))

  # ------------------------------------------------ #
  # Calculate distance (edit distance)
  # ------------------------------------------------ #
  # Biostrings package
  time_start2 <- proc.time()[3]
  amplicon_mt <- stringDist(amplicon_aln, method = "levenshtein")
  message(sprintf("Distance calculated: %.1f sec elapsed\n", proc.time()[3] - time_start2))
  
  # ------------------------------------------------ #
  # t-SNE
  # ------------------------------------------------ #
  ## Perform t-SNE
  tsne <- Rtsne(amplicon_mt, is_distance = TRUE, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500)
  ## Compile results
  tsne_df <- cbind(md_sub, tsne$Y)
  colnames(tsne_df)[(ncol(tsne_df)-1):ncol(tsne_df)] <- c("x", "y")
  saveRDS(list(tsne_df, tsne), sprintf("%s/tsne_res_%s_n%s.obj", output_folder, primer_i, N_SAMPLE))
  saveRDS(amplicon_mt, sprintf("%s/amplicon_mts_%s_n%s.obj", output_folder, primer_i, N_SAMPLE))

  # ------------------------------------------------ #
  # Visualize
  # ------------------------------------------------ #
  g1 <- tsne_df %>% ggplot(aes(x = x, y = y, col = cat)) +
    geom_point(size = 3, alpha = 0.5) +
    scale_color_manual(values = c("deepskyblue1", "navyblue", "orangered2"), name = "Taxa") +
    ggtitle(sprintf("t-SNE visualization of %s sequences amplified by %s", N_SAMPLE, primer_i)) +
    xlab("Axis 1") + ylab("Axis 2") +
    theme_linedraw() + theme(panel.grid = element_blank())
  
  
  # ------------------------------------------------ #
  # Save results
  # ------------------------------------------------ #
  ggsave(sprintf("%s/%s_n%s.pdf", output_folder, primer_i, N_SAMPLE), g1,
         width = 10, height = 9)
  
  message(sprintf("\`%s\` analyzed: %.1f sec elapsed", primer_i, proc.time()[3] - time_start))
}


# ------------------------------------------------ #
# Save results
# ------------------------------------------------ #
#save(list = ls(all.names = TRUE), file = sprintf("%s/%s.RData", output_folder, output_folder))

# Save sessioninfo
macam::save_session_info()
