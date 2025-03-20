####
#### Evaluate interspecific variations for marine mammals
#### 2023.07.20 Ushio
####

# Load libraries
library(tidyverse); packageVersion("tidyverse") # 1.3.2, 2023.7.18
library(cowplot); packageVersion("cowplot") # 1.1.1, 2023.7.18
library(Biostrings); packageVersion("Biostrings") # 2.66.0, 2023.7.18
#library(macam); packageVersion("macam") # 0.1.4, 2023.7.18
theme_set(theme_cowplot())

# Create output directory
set.seed(1234)
(output_folder <- macam::outdir_create())

# Load sequence metadata
md <- read.csv("data_dbseq/marine_mammal_w_nontarget.csv")


# ------------------------------------------------ #
# Load primers
# ------------------------------------------------ #
# Load tested primer names
primer_tested <- list.files("01_usearchOut/1_EachPrimerSet/")
# Drop some primer candidates
## There are dropped due to unexpected amplifications
primer_tested <- primer_tested[!primer_tested %in% c("CheckEditDistance_Cetacea", "CheckEditDistance_MarineMammal")]
## There are dropped due to unexpected amplifications
primer_tested <- primer_tested[primer_tested != "Mu2743F_Mu2950R" & primer_tested != "Mu15826F_Mu15796R"]

# Prepare the final object
var_all_df <- list()


# ------------------------------------------------ #
# Main looop
# ------------------------------------------------ #
for (primer_i in primer_tested) {
  # Load sequences
  primer <- readDNAStringSet(sprintf("data_primerset/each_primerset/%s.fa", primer_i))
  amplicon <- readDNAStringSet(sprintf("01_usearchOut/1_EachPrimerSet/CheckEditDistance_MarineMammal/%s/hits_all.fa", primer_i))
  
  # Extract meta data
  ncbi_id <- names(amplicon) %>% str_split("/") %>% sapply('[', 1) %>% str_sub(start = -11)
  md_sub <- md[match(ncbi_id, md$ncbi_id),]
  md_sub <- md_sub[md_sub$Taxa_ID %>% order,]
  # Sort sequence so that it matches metadata
  amplicon <- amplicon[match(md_sub$ncbi_id, names(amplicon) %>% str_split("/") %>% sapply('[', 1) %>% str_sub(start = -11))]
  
  # Remove primers
  n_fp <- width(primer)[1]
  n_rp <- width(primer)[2]
  amplicon <- subseq(amplicon, start = n_fp+1, end = -n_rp-1)
  writeXStringSet(amplicon, sprintf("%s/amp_tmp.fa", output_folder))
  
  # Alignment by MAFFT
  system(sprintf("mafft %s/amp_tmp.fa > %s/amp_aligned.fa", output_folder, output_folder))
  amplicon_aln <- readDNAStringSet(sprintf("%s/amp_aligned.fa", output_folder))
  system(sprintf("rm %s/amp_tmp.fa %s/amp_aligned.fa", output_folder, output_folder))
  
  # ------------------------------------------------ #
  # Evaluate interspecific variations
  # ------------------------------------------------ #
  # Biostrings package
  amplicon_mt <- stringDist(amplicon_aln, method = "levenshtein") %>% as.matrix
  colnames(amplicon_mt) <- rownames(amplicon_mt) <- md_sub$Genus_Species

  # Prepare the summary object
  var_df <- data.frame(n_diff_label = as.character(0:6), n_diff = 0:6)
  var_df$n_diff_label[7] <- ">5"
  var_df$species <- 0
  var_df$genus <- 0
  var_df$family <- 0
  
  # Inter-species variations
  # Collect and match values
  dist_val <- amplicon_mt[upper.tri(amplicon_mt)]
  for (k in 0:6) var_df[var_df$n_diff == k,"species"] <- sum(dist_val == k, na.rm = T)
  var_df[7,"species"] <- sum(dist_val > 6, na.rm = T)
  
  # Inter-genus variations
  genus_mt <- matrix(NA, ncol = ncol(amplicon_mt), nrow = nrow(amplicon_mt))
  colnames(genus_mt) <- rownames(genus_mt) <- md_sub$Genus_Species
  for (i in 1:nrow(genus_mt)) {
    for (j in 1:ncol(genus_mt)) {
      if ((j > i) & (md_sub$Genus[i] != md_sub$Genus[j])) genus_mt[i,j] <- amplicon_mt[i,j]
    }
  }
  # Collect and match values
  dist_val <- genus_mt[upper.tri(genus_mt)]
  for (k in 0:6) var_df[var_df$n_diff == k,"genus"] <- sum(dist_val == k, na.rm = T)
  var_df[7,"genus"] <- sum(dist_val > 6, na.rm = T)
  
  # Inter-family variations
  fami_mt <- matrix(NA, ncol = ncol(amplicon_mt), nrow = nrow(amplicon_mt))
  colnames(fami_mt) <- rownames(fami_mt) <- md_sub$Genus_Species
  for (i in 1:nrow(fami_mt)) {
    for (j in 1:ncol(fami_mt)) {
      if ((j > i) & (md_sub$Family[i] != md_sub$Family[j])) fami_mt[i,j] <- amplicon_mt[i,j]
    }
  }
  # Collect and match values
  dist_val <- fami_mt[upper.tri(fami_mt)]
  for (k in 0:6) var_df[var_df$n_diff == k,"family"] <- sum(dist_val == k, na.rm = T)
  var_df[7,"family"] <- sum(dist_val > 6, na.rm = T)
  rm(dist_val)
  
  ## Primer name
  var_df$primer <- primer_i

  # Store in the final object
  var_all_df <- c(var_all_df, list(var_df))
  names(var_all_df)[length(var_all_df)] <- primer_i
  
  # Message
  cat(sprintf("Sequences amplified by the primer \`%s\` analyzed`\n\n", primer_i))
}


# ------------------------------------------------ #
# Visualize results
# ------------------------------------------------ #
# Primer order
primer_order <- c("MiMammal", "MarVer1", "MarVer2", "MarVer3", "Ceto2", "Riaz12S",
                  "Mu31F_Dc320R", "Mu31F_MiMammalR", "Dc321F_Dc495R", "Dc494F_Mu643R",
                  "Mu642F_Dc1015R", "Mu642F_Mu1021R", "Dc671F_Dc1015R", "Dc671F_Mu1021R",
                  "Dc1458_1638", "Dc1465F_Dc1646R", "Mu2084F_Dc2438R", "Dc2173_2580",
                  "Mu2187F_Dc2438R", "Mu2187F_Mu2563R", "Dc2385_2580", "Dc2430_2580",
                  "Dc2438F_Mu2563R", "Dc2505_2580", "Mu9459F_Mu9822R", "Mu9459F_Mu9834R")

# Flatten the list
var_all_df2 <- list_rbind(var_all_df) %>% 
  pivot_longer(cols = -c(n_diff, n_diff_label, primer), names_to = "resolution", values_to = "freq")
# Adjust labels
var_all_df2$n_diff_label <- factor(var_all_df2$n_diff_label, levels = c(as.character(0:5), ">5"))
var_all_df2$resolution <- factor(var_all_df2$resolution, levels = c("species", "genus", "family"))
var_all_df2$primer <- factor(var_all_df2$primer, levels = primer_order)

# Visualize
g1 <- var_all_df2 %>%
  ggplot(aes(x = n_diff_label, y = freq + 1, fill = resolution)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("lightblue", "burlywood2", "purple1"), name = "Resolution") +
  facet_wrap(. ~ primer) + panel_border() +
  xlab("Edit distance (bases)") + ylab("N of combinations + 1") +
  theme(strip.text = element_text(size = 10),
        panel.grid.major.y = element_line(linewidth = 0.1, color = "gray70")) +
  scale_y_log10() + panel_border() +
  ggtitle("Edit distance between amplicons of marine mammals")

g2 <- var_all_df2 %>% filter(n_diff < 6) %>% 
  ggplot(aes(x = n_diff_label, y = freq, fill = resolution)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("lightblue", "burlywood2", "purple1"), name = "Resolution") +
  facet_wrap(. ~ primer) + panel_border() +
  xlab("Edit distance (bases)") + ylab("N of combinations") +
  theme(strip.text = element_text(size = 10),
        panel.grid.major.y = element_line(linewidth = 0.1, color = "gray70")) +
  panel_border() + ggtitle("Edit distance between amplicons of marine mammals")



# ------------------------------------------------ #
# Save results
# ------------------------------------------------ #
# Save results
saveRDS(var_all_df, sprintf("%s/var_all_df.obj", output_folder))
ggsave(sprintf("%s/VariabilityMarine.pdf", output_folder),
       plot = g1, width = 12, height = 10)
ggsave(sprintf("%s/VariabilityMarine2.pdf", output_folder),
       plot = g2, width = 12, height = 10)

# Save workspace
#save(list = ls(all.names = TRUE), file = sprintf("%s/%s.RData", output_folder, output_folder))

# Save sessioninfo
macam::save_session_info()
