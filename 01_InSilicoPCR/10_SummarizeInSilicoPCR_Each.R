####
#### Summarize in silico PCR (each primer analysis)
#### 2025.01.24 Ushio
####

# Load libraries
library(tidyverse); packageVersion("tidyverse") # 1.3.2, 2023.7.18
library(stringr); packageVersion("stringr") # 1.5.0, 2023.7.18
library(Biostrings); packageVersion("Biostrings") # 2.66.0, 2023.7.18
library(macam); packageVersion("macam") # 0.1.9, 2025.1.31

# Create output directory
set.seed(1234)
dir.create("00_SessionInfo")
(output_folder <- macam::outdir_create())


# ------------------------------------------------ #
# Preparations
# ------------------------------------------------ #
acc_taxdb <- read.csv("02_CompileTaxaOut/acc_tax_db_ed.csv")

# Check acc_taxdb
primer_tested <- list.files("01_usearchOut/2_EachPrimer/")
# Column names from https://www.drive5.com/usearch/manual/pcrout.html
res_colnames <- c("query_label", "primer_name", "strand", "n_diff",
                  "start_loc", "end_loc", "alighment")
# Prepare summary object
match_seq_all <- data.frame()
match_sp_all <- data.frame()
match_seq_detail <- list()
match_sp_detail <- list()
check_n_amplicon <- data.frame(primer_name = primer_tested)
check_n_amplicon$n_amp_seq <- NA
check_n_amplicon$n_amp_region <- NA


# ------------------------------------------------ #
# Main loop
# ------------------------------------------------ #
for (file_i in primer_tested) {
  # Load result files
  res_all <- read.csv(sprintf("01_usearchOut/2_EachPrimer/%s/hits_all_text.txt", file_i),
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
      res_dup_id <- which(res_all$accession_id == (dup_acc_df[i,] %>% pull(accession_id)))
      # Duplicated data
      res_dup <- res_all %>% filter(accession_id == (dup_acc_df[i,] %>% pull(accession_id)))
      # Remove dupulicated data (keep minimum total_diff record)
      res_all <- res_all[-res_dup_id[-which.min(res_dup$n_diff)],]
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
  match_tab <- table(res_all$n_diff, res_all$cat) %>% data.frame
  
  # Compile Var1 and Var2
  match_tab$Var1 <- match_tab$Var1 %>% as.character %>% as.numeric

  # Convert to a wider form
  match_wide <- match_tab %>%
    pivot_wider(names_from = Var2, values_from = Freq) %>% 
    replace_na(list(cetacea = 0, vertebrate = 0, fish = 0)) %>% data.frame
  # Add categories if not exist
  missing_id <- which(!c("cetacea", "vertebrate", "fish") %in% unique(match_tab$Var2))
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
  
  
  # ------------------------------------------- #
  # Summarize match-mismatch (N of species)
  # ------------------------------------------- #
  n_sp_all <- res_all %>% group_by(n_diff, species) %>% summarize(total_hit = n())
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
      n_sp_all <- n_sp_all[-sp_dup_id[-which.min(sp_dup$n_diff)],]
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

  # Prepare "match" tab
  match_sp_tab <- table(n_sp_all$n_diff, n_sp_all$cat) %>% data.frame
  
  # Compile Var1 and Var2
  match_sp_tab$Var1 <- match_sp_tab$Var1 %>% as.character %>% as.numeric

  # Convert to a wider form
  match_sp_wide <- match_sp_tab %>%
    pivot_wider(names_from = Var2, values_from = Freq) %>% 
    replace_na(list(cetacea = 0, vertebrate = 0, fish = 0)) %>% data.frame
  
  # Add categories if not exist
  missing_id2 <- which(!c("cetacea", "vertebrate", "fish") %in% unique(match_sp_tab$Var2))
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
  
  
  # ------------------------------------------- #
  # Prepare the summary stats for main figures
  # ------------------------------------------- #
  # Rename columns
  match_tab$n_diff <- match_tab$Var1
  match_sp_tab$n_diff <- match_sp_tab$Var1
  
  # Summarize sequence hits
  match_seq <- match_tab %>% group_by(n_diff, Var2) %>%
    summarize(total_hit_seq = sum(Freq)) %>% 
    pivot_wider(names_from = Var2, values_from = total_hit_seq) %>% 
    replace_na(list(cetacea = 0, vertebrate = 0, fish = 0)) %>%
    mutate(primer_name = file_i)
  match_seq_all <- rbind(match_seq_all, match_seq)

  # Summarize species hits
  match_sp <- match_sp_tab %>% group_by(n_diff, Var2) %>%
    summarize(total_hit_sp = sum(Freq)) %>% 
    pivot_wider(names_from = Var2, values_from = total_hit_sp) %>% 
    replace_na(list(cetacea = 0, vertebrate = 0, fish = 0)) %>%
    mutate(primer_name = file_i)
  match_sp_all <- rbind(match_sp_all, match_sp)

  # Keep detailed results
  match_seq_detail <- c(match_seq_detail, list(res_all))
  match_sp_detail <- c(match_sp_detail, list(n_sp_all))
  names(match_seq_detail)[length(match_seq_detail)] <- file_i
  names(match_sp_detail)[length(match_sp_detail)] <- file_i
}


# Check duplicated amplifications
check_n_amplicon$n_dup <- check_n_amplicon$n_amp_region - check_n_amplicon$n_amp_seq

# Clean up global environments
rm(match_sp); rm(match_seq); rm(file_i)
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
saveRDS(match_seq_detail, sprintf("%s/match_seq_detail.obj", output_folder))
saveRDS(match_sp_detail, sprintf("%s/match_sp_detail.obj", output_folder))
saveRDS(match_seq_all, sprintf("%s/match_seq_all.obj", output_folder))
saveRDS(match_sp_all, sprintf("%s/match_sp_all.obj", output_folder))

# Save workspace
#save(list = ls(all.names = TRUE), file = sprintf("%s/%s.RData", output_folder, output_folder))

# Save sessioninfo
macam::save_session_info()
