####
#### Summarize in silico PCR
#### 2023.07.19 Ushio
####

# Load libraries
library(tidyverse); packageVersion("tidyverse") # 1.3.2, 2023.7.18
library(stringr); packageVersion("stringr") # 1.5.0, 2023.7.18
library(macam); packageVersion("macam") # 0.1.4, 2023.7.18
library(Biostrings); packageVersion("Biostrings") # 2.66.0, 2023.7.18
library(scatterpie); packageVersion("scatterpie") # 0.2.1, 2023.7.18
#library(macam); packageVersion("macam") # 0.1.4, 2023.7.18

# Create output directory
set.seed(1234)
dir.create("00_SessionInfo")
(output_folder <- macam::outdir_create())


# ------------------------------------------------ #
# Preparations
# ------------------------------------------------ #
acc_taxdb <- read.csv("02_CompileTaxaOut/acc_tax_db_ed.csv") # <- The original tax_db was manually edited
if (F) {
  # Amend taxdb
  acc_taxdb[which(acc_taxdb$order == "Testudines"), "class"] <- "Reptilia" # Turtle
  acc_taxdb[which(acc_taxdb$order == "Crocodylia"), "class"] <- "Reptilia" # Crocodile
  acc_taxdb[which(acc_taxdb$order == "Ceratodontiformes"), "class"] <- "Ceratodontiformes" # Fish-like
  acc_taxdb[which(acc_taxdb$order == "Coelacanthiformes"), "class"] <- "Coelacanthiformes" # Fish-like
  fish_taxa_names <- c("Hyperoartia", "Myxini", "Chondrichthyes", "Ceratodontiformes", "Coelacanthiformes", "Actinopteri", "Cladistia")
  acc_taxdb$rep_tax <- NA; acc_taxdb$cat <- "vertebrate"
  # Assign rep_tax
  ## rep_tax
  acc_taxdb$rep_tax <- acc_taxdb$class
  acc_taxdb$rep_tax[which(acc_taxdb$class %in% fish_taxa_names)] <- "fish"
  acc_taxdb$rep_tax[which(acc_taxdb$infraorder == "Cetacea")] <- "cetacea"
  acc_taxdb$rep_tax[which(acc_taxdb$class == "Mammalia")] <- "mammal"
  acc_taxdb$rep_tax[which(acc_taxdb$class == "Aves")] <- "bird"
  acc_taxdb$rep_tax[which(acc_taxdb$class == "Amphibia")] <- "amphibian"
  acc_taxdb$rep_tax[which(acc_taxdb$class == "Reptilia")] <- "reptile"
  acc_taxdb$rep_tax[which(acc_taxdb$class == "Lepidosauria")] <- "reptile"
  acc_taxdb$rep_tax[which(acc_taxdb$class == "Leptocardii")] <- "other"
  ## cat
  acc_taxdb$cat[which(acc_taxdb$class %in% fish_taxa_names)] <- "fish"
  acc_taxdb$cat[which(acc_taxdb$infraorder == "Cetacea")] <- "cetacea"
  write.csv(acc_taxdb, "02_CompileTaxaOut/acc_tax_db_ed.csv", row.names = F)
}

# Check acc_taxdb
primer_tested <- list.files("01_usearchOut/1_EachPrimerSet/")
primer_tested <- primer_tested[!primer_tested %in% c("CheckEditDistance_Cetacea", "CheckEditDistance_MarineMammal")]
# Remove two primer sets that seem unsuitable
primer_tested <- primer_tested[primer_tested != "Mu2743F_Mu2950R" & primer_tested != "Mu15826F_Mu15796R"]
radius_fac1 <- 0.00005
radius_fac2 <- 0.0002
# Column names from https://www.drive5.com/usearch/manual/pcrout.html
res_colnames <- c("query_label", "query_start", "query_end", "query_len",
                  "forward_primer", "strand_fprimer", "alignment_fprimer",
                  "reverse_primer", "strand_rprimer", "alignment_rprimer",
                  "amplicon_len", "amplicon", "n_diff_fprimer", "n_diff_rprimer", "total_diff")
# Prepare summary object
match_seq_all <- data.frame()
match_sp_all <- data.frame()
match_seq_detail <- list()
match_sp_detail <- list()
amplicon_len_list <- list()
amplicon_len_others_list <- list()
check_n_amplicon <- data.frame(primer_name = primer_tested)
check_n_amplicon$n_amp_seq <- NA
check_n_amplicon$n_amp_region <- NA


# ------------------------------------------------ #
# Main loop
# ------------------------------------------------ #
for (file_i in primer_tested) {
  # Load result files
  res_all <- read.csv(sprintf("01_usearchOut/1_EachPrimerSet/%s/hits_all_text.txt", file_i),
                      sep = "\t", header = F, col.names = res_colnames)
  
  # Extract Accession numbers
  res_all$accession_id <- res_all$query_label %>% str_replace("UNVERIFIED: ", "") %>% str_split(" ") %>% sapply(`[`, 1)

  # Dereplicate data
  # (i.e., sometimes two regions of one sequence are amplified)
  # In that case, we remove the amplification with more mismatches
  dup_acc_df <- res_all %>% group_by(accession_id) %>% summarize(n = n()) %>% filter(n > 1)
  check_n_amplicon[check_n_amplicon$primer_name == file_i, "n_amp_region"] <- nrow(res_all)
  check_n_amplicon[check_n_amplicon$primer_name == file_i, "n_amp_seq"] <- nrow(res_all) - (sum(dup_acc_df$n) - nrow(dup_acc_df))
  if (nrow(dup_acc_df) > 0) {
    for (i in 1:nrow(dup_acc_df)) {
      # Duplicated IDs
      res_dup_id <- which(res_all$accession_id == dup_acc_df[i,] %>% pull(accession_id))
      # Duplicated data
      res_dup <- res_all %>% filter(accession_id == dup_acc_df[i,] %>% pull(accession_id))
      # Remove dupulicated data (keep minimum total_diff record)
      res_all <- res_all[-res_dup_id[-which.min(res_dup$total_diff)],]
    }
  }

  # Add taxonomy information
  res_all_taxdf <- match(res_all$accession_id, acc_taxdb$accession_id) %>% acc_taxdb[.,]
  res_all$species <- res_all_taxdf$species
  res_all$tax_id <- res_all_taxdf$tax_id
  res_all$class <- res_all_taxdf$class
  res_all$infraorder <- res_all_taxdf$infraorder
  res_all$order <- res_all_taxdf$order
  res_all$family <- res_all_taxdf$family
  res_all$genus <- res_all_taxdf$genus
  
  # Set target-nontarget categories
  res_all$rep_tax <- NA
  res_all$cat <- "vertebrate"
  fish_taxa_names <- c("Hyperoartia", "Myxini", "Chondrichthyes", "Ceratodontiformes", "Coelacanthiformes", "Actinopteri", "Cladistia")
  
  # Assign rep_tax
  ## rep_tax
  res_all$rep_tax <- res_all$class
  res_all$rep_tax[which(res_all$class %in% fish_taxa_names)] <- "fish"
  res_all$rep_tax[which(res_all$infraorder == "Cetacea")] <- "cetacea"
  res_all$rep_tax[which(res_all$class == "Mammalia")] <- "mammal"
  res_all$rep_tax[which(res_all$class == "Aves")] <- "bird"
  res_all$rep_tax[which(res_all$class == "Amphibia")] <- "amphibian"
  res_all$rep_tax[which(res_all$class == "Reptilia")] <- "reptile"
  res_all$rep_tax[which(res_all$class == "Lepidosauria")] <- "reptile"
  res_all$rep_tax[which(res_all$class == "Leptocardii")] <- "other"
  ## cat
  res_all$cat[which(res_all$class %in% fish_taxa_names)] <- "fish"
  res_all$cat[which(res_all$infraorder == "Cetacea")] <- "cetacea"
  
  
  # ------------------------------------------- #
  # Summarize match-mismatch (N of sequences)
  # ------------------------------------------- #
  match_tab <- table(res_all$n_diff_fprimer, res_all$n_diff_rprimer, res_all$cat) %>% data.frame
  
  # Compile Var1 and Var2
  match_tab$Var1 <- match_tab$Var1 %>% as.character %>% as.numeric
  match_tab$Var2 <- match_tab$Var2 %>% as.character %>% as.numeric
  match_tab$Var4 <- paste0(match_tab$Var1, "-", match_tab$Var2)
  match_tab$n_diff <- match_tab$Var1 + match_tab$Var2
  
  # Convert to a wider form
  match_wide <- match_tab %>%
    pivot_wider(names_from = Var3, values_from = Freq) %>% 
    replace_na(list(cetacea = 0, vertebrate = 0, fish = 0)) %>% data.frame
  # Add categories if not exist
  missing_id <- which(!c("cetacea", "vertebrate", "fish") %in% unique(match_tab$Var3))
  if(length(missing_id) == 1) {
    if (missing_id == 1) match_wide$cetacea <- 0
    if (missing_id == 2) match_wide$vertebrate <- 0
    if (missing_id == 3) match_wide$fish <- 0
  } else if (length(missing_id) == 2) {
    if (!1 %in% missing_id) { match_wide$vertebrate <- 0; match_wide$fish <- 0 }
    if (!2 %in% missing_id) { match_wide$cetacea <- 0; match_wide$fish <- 0 }
    if (!3 %in% missing_id) { match_wide$cetacea <- 0; match_wide$vertebrate <- 0 }
  } else if (length(missing_id) == 3) {
    match_wide$cetacea <- 0; match_wide$vertebrate <- 0; match_wide$fish <- 0
  }
  # Calculate proportion
  match_wide$total_hit <- match_wide$cetacea + match_wide$vertebrate + match_wide$fish
  match_wide$vertebrate_p <- match_wide$vertebrate/match_wide$total_hit
  match_wide$cetacea_p <- match_wide$cetacea/match_wide$total_hit
  match_wide$fish_p <- match_wide$fish/match_wide$total_hit
  match_wide$vertebrate_p <- match_wide %>% pull(vertebrate_p) %>% replace_na(0)
  match_wide$cetacea_p <- match_wide %>% pull(cetacea_p) %>% replace_na(0)
  match_wide$fish_p <- match_wide %>% pull(fish_p) %>% replace_na(0)
  
  # Draw scattered pie chart
  g1 <- ggplot() + geom_scatterpie(data = match_wide,
                                   aes(x = Var1, y = Var2, r = radius_fac1*total_hit),
                                   cols = c("cetacea_p", "fish_p", "vertebrate_p")) +
    ggtitle(paste0("The number of sequences amplified by ", file_i)) +
    scale_fill_manual(values = c("red2", "lightblue3", "burlywood2"),
                      name = c("Taxa"),
                      label = c("Cetacea", "Fish", "Other vertebrates")) +
    xlab("N of mismatches to the forward primer") +
    ylab("N of mismatches to the reverse primer") + 
    theme_light() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  ggsave(sprintf("%s/seq_%s.pdf", output_folder, file_i), plot = g1, width = 6.2, height = 5)
  
  
  # ------------------------------------------- #
  # Summarize match-mismatch (N of species)
  # ------------------------------------------- #
  n_sp_all <- res_all %>% group_by(n_diff_fprimer, n_diff_rprimer, species) %>%
    summarize(total_hit = n())
  n_sp_all$class <- acc_taxdb[match(n_sp_all$species, acc_taxdb$species), "class"]
  n_sp_all$order <- acc_taxdb[match(n_sp_all$species, acc_taxdb$species), "order"]
  n_sp_all$infraorder <- acc_taxdb[match(n_sp_all$species, acc_taxdb$species), "infraorder"]
  n_sp_all$family <- acc_taxdb[match(n_sp_all$species, acc_taxdb$species), "family"]
  n_sp_all$rep_tax <- NA
  n_sp_all$cat <- "vertebrate"
  
  # Remove duplicated species
  dup_sp_df <- n_sp_all %>% group_by(species) %>% summarize(n = n()) %>% filter(n > 1)
  if (nrow(dup_sp_df) > 0) {
    for (i in 1:nrow(dup_sp_df)) {
      # Duplicated IDs
      sp_dup_id <- which(n_sp_all$species == dup_sp_df[i,] %>% pull(species))
      # Duplicated data
      sp_dup <- n_sp_all %>% filter(species == dup_sp_df[i,] %>% pull(species))
      # Remove dupulicated data (keep minimum total_diff record)
      n_sp_all <- n_sp_all[-sp_dup_id[-which.min(sp_dup$n_diff_fprimer + sp_dup$n_diff_rprimer)],]
    }
  }
  
  # Assign rep_tax
  ## rep_tax
  n_sp_all$rep_tax <- n_sp_all$class
  n_sp_all$rep_tax[which(n_sp_all$class %in% fish_taxa_names)] <- "fish"
  n_sp_all$rep_tax[which(n_sp_all$infraorder == "Cetacea")] <- "cetacea"
  n_sp_all$rep_tax[which(n_sp_all$class == "Mammalia")] <- "mammal"
  n_sp_all$rep_tax[which(n_sp_all$class == "Aves")] <- "bird"
  n_sp_all$rep_tax[which(n_sp_all$class == "Amphibia")] <- "amphibian"
  n_sp_all$rep_tax[which(n_sp_all$class == "Reptilia")] <- "reptile"
  n_sp_all$rep_tax[which(n_sp_all$class == "Lepidosauria")] <- "reptile"
  n_sp_all$rep_tax[which(n_sp_all$class == "Leptocardii")] <- "other"
  ## cat
  n_sp_all$cat[which(n_sp_all$class %in% fish_taxa_names)] <- "fish"
  n_sp_all$cat[which(n_sp_all$infraorder == "Cetacea")] <- "cetacea"
  #dim(n_sp_all)

  match_sp_tab <- table(n_sp_all$n_diff_fprimer, n_sp_all$n_diff_rprimer, n_sp_all$cat) %>% data.frame
  
  # Compile Var1 and Var2
  match_sp_tab$Var1 <- match_sp_tab$Var1 %>% as.character %>% as.numeric
  match_sp_tab$Var2 <- match_sp_tab$Var2 %>% as.character %>% as.numeric
  match_sp_tab$Var4 <- paste0(match_sp_tab$Var1, "-", match_sp_tab$Var2)
  match_sp_tab$n_diff <- match_sp_tab$Var1 + match_sp_tab$Var2
  
  # Convert to a wider form
  match_sp_wide <- match_sp_tab %>%
    pivot_wider(names_from = Var3, values_from = Freq) %>% 
    replace_na(list(cetacea = 0, vertebrate = 0, fish = 0)) %>% data.frame
  
  # Add categories if not exist
  missing_id2 <- which(!c("cetacea", "vertebrate", "fish") %in% unique(match_sp_tab$Var3))
  if(length(missing_id2) == 1) {
    if (missing_id2 == 1) match_sp_wide$cetacea <- 0
    if (missing_id2 == 2) match_sp_wide$vertebrate <- 0
    if (missing_id2 == 3) match_sp_wide$fish <- 0
  } else if (length(missing_id2) == 2) {
    if (!1 %in% missing_id2) { match_sp_wide$vertebrate <- 0; match_sp_wide$fish <- 0 }
    if (!2 %in% missing_id2) { match_sp_wide$cetacea <- 0; match_sp_wide$fish <- 0 }
    if (!3 %in% missing_id2) { match_sp_wide$cetacea <- 0; match_sp_wide$vertebrate <- 0 }
  } else if (length(missing_id2) == 3) {
    match_sp_wide$cetacea <- 0; match_sp_wide$vertebrate <- 0; match_sp_wide$fish <- 0
  }
  
  # Calculate proportion
  match_sp_wide$total_hit <- match_sp_wide$cetacea + match_sp_wide$vertebrate + match_sp_wide$fish
  match_sp_wide$vertebrate_p <- match_sp_wide$vertebrate/match_sp_wide$total_hit
  match_sp_wide$cetacea_p <- match_sp_wide$cetacea/match_sp_wide$total_hit
  match_sp_wide$fish_p <- match_sp_wide$fish/match_sp_wide$total_hit
  match_sp_wide$vertebrate_p <- match_sp_wide %>% pull(vertebrate_p) %>% replace_na(0)
  match_sp_wide$cetacea_p <- match_sp_wide %>% pull(cetacea_p) %>% replace_na(0)
  match_sp_wide$fish_p <- match_sp_wide %>% pull(fish_p) %>% replace_na(0)
  
  # Draw scattered pie chart
  g2 <- ggplot() + geom_scatterpie(data = match_sp_wide,
                                   aes(x = Var1, y = Var2, r = radius_fac2*total_hit),
                                   cols = c("cetacea_p", "fish_p", "vertebrate_p")) +
    ggtitle(paste0("The number of species amplified by ", file_i)) +
    scale_fill_manual(values = c("red2", "lightblue3", "burlywood2"),
                      name = c("Taxa"),
                      label = c("Cetacea", "Fish", "Other vertebrates")) +
    xlab("N of mismatches to the forward primer") +
    ylab("N of mismatches to the reverse primer") + 
    theme_light() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  ggsave(sprintf("%s/spp_%s.pdf", output_folder, file_i), plot = g2, width = 6.2, height = 5)
  
  
  # ------------------------------------------- #
  # Prepare the summary stats for main figures
  # ------------------------------------------- #
  # Summarize sequence hits
  match_seq <- match_tab %>% group_by(n_diff, Var3) %>%
    summarize(total_hit_seq = sum(Freq)) %>% 
    pivot_wider(names_from = Var3, values_from = total_hit_seq) %>% 
    replace_na(list(cetacea = 0, vertebrate = 0, fish = 0)) %>%
    mutate(primer_name = file_i)
  match_seq_all <- rbind(match_seq_all, match_seq)

  # Summarize species hits
  match_sp <- match_sp_tab %>% group_by(n_diff, Var3) %>%
    summarize(total_hit_sp = sum(Freq)) %>% 
    pivot_wider(names_from = Var3, values_from = total_hit_sp) %>% 
    replace_na(list(cetacea = 0, vertebrate = 0, fish = 0)) %>%
    mutate(primer_name = file_i)
  match_sp_all <- rbind(match_sp_all, match_sp)

  # Keep detailed results
  match_seq_detail <- c(match_seq_detail, list(res_all))
  match_sp_detail <- c(match_sp_detail, list(n_sp_all))
  names(match_seq_detail)[length(match_seq_detail)] <- file_i
  names(match_sp_detail)[length(match_sp_detail)] <- file_i
  
  # Collect amplicon length
  amplicon_len_list <- c(amplicon_len_list, list(c(res_all %>% filter(cat == "cetacea") %>% select(amplicon_len))))
  amplicon_len_others_list <- c(amplicon_len_others_list, list(c(res_all %>% filter(cat != "cetacea") %>% select(amplicon_len))))
  names(amplicon_len_list)[length(amplicon_len_list)] <- file_i
  names(amplicon_len_others_list)[length(amplicon_len_others_list)] <- file_i
}


# Check duplicated amplifications
check_n_amplicon$n_dup <- check_n_amplicon$n_amp_region - check_n_amplicon$n_amp_seq

# Clean up global environments
rm(match_sp); rm(match_seq)
rm(radius_fac1); rm(radius_fac2)
rm(g1); rm(g2); rm(file_i)
rm(sp_dup); rm(sp_dup_id); rm(res_dup); rm(res_dup_id) 
rm(dup_acc_df); rm(dup_sp_df)
rm(match_tab); rm(match_wide)
rm(match_sp_tab); rm(match_sp_wide)
rm(n_sp_all); rm(res_all)


# ------------------------------------------------ #
# Save results
# ------------------------------------------------ #
# Save results
write.csv(check_n_amplicon, sprintf("%s/check_n_amplicon.csv", output_folder), row.names = F)
saveRDS(amplicon_len_list, sprintf("%s/amplicon_len_list.obj", output_folder))
saveRDS(amplicon_len_others_list, sprintf("%s/amplicon_len_others_list.obj", output_folder))
saveRDS(match_seq_detail, sprintf("%s/match_seq_detail.obj", output_folder))
saveRDS(match_sp_detail, sprintf("%s/match_sp_detail.obj", output_folder))
saveRDS(match_seq_all, sprintf("%s/match_seq_all.obj", output_folder))
saveRDS(match_sp_all, sprintf("%s/match_sp_all.obj", output_folder))

# Save sessioninfo
macam::save_session_info()

