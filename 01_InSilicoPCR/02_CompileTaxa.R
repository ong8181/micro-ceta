####
#### Compile taxa information before summarizing in silico PCR
#### 2023.07.19 Ushio
####

# Load libraries
library(tidyverse); packageVersion("tidyverse") # 1.3.2, 2023.7.18
library(stringr); packageVersion("stringr") # 1.5.0, 2023.7.18
library(Biostrings); packageVersion("Biostrings") # 2.66.0, 2023.7.18
library(scatterpie); packageVersion("scatterpie") # 0.2.1, 2023.7.18
library(taxize); packageVersion("taxize") # 0.9.100, 2023.7.19
library(macam); packageVersion("macam") # 0.1.4, 2023.7.18

# Create output directory
set.seed(1234)
dir.create("00_SessionInfo")
(output_folder <- macam::outdir_create())


# ------------------------------------------------ #
# Generate Accession ID - TaxID - Taxa database
# ------------------------------------------------ #

# ------------------------------------------------ #
# Retrieve tax ID
# ------------------------------------------------ #
# Original database
all_seq <- Biostrings::readDNAStringSet("data_dbseq/vertebrate_db.fa")

# Retrieve accession IDs and species names
acc_ids <- names(all_seq) %>% str_replace("UNVERIFIED: ", "") %>% str_split(" ") %>% sapply(`[`, 1)
genus_names <- names(all_seq) %>% str_replace("\\(", "") %>% str_replace("\\)", "") %>%
  str_replace("UNVERIFIED: ", "") %>% str_split(" ") %>% sapply(`[`, 2)
species_names <- names(all_seq) %>% str_replace("\\(", "") %>% str_replace("\\)", "") %>%
  str_replace("UNVERIFIED: ", "") %>% str_split(" ") %>% sapply(`[`, 3)
acc_tax_db <- data.frame(accession_id = acc_ids, species = paste(genus_names, species_names))

# Retrieve TaxID and Taxa information
tax_id_list <- genbank2uid(acc_tax_db$accession_id) # Batch-size = 100 (~ 500 sec?)
names(tax_id_list) <- acc_tax_db$accession_id
tax_id <- tax_id_list %>% sapply('[', 1)
acc_tax_db$tax_id <- tax_id
saveRDS(tax_id_list, sprintf("%s/tax_id_list.obj", output_folder))
saveRDS(tax_id, sprintf("%s/tax_id.obj", output_folder))


# ------------------------------------------------ #
# Retrieve tax information
# ------------------------------------------------ #
# Compile taxa_list
tax_list <- list()
tax_id_split <- split(1:length(tax_id), ceiling(seq_along(1:length(tax_id)) / 50))
for (tax_i in 1:length(tax_id_split)) {
  result <- tryCatch({
    classification(tax_id[tax_id_split[[tax_i]]], db = "ncbi") # Batch-size = 50 (~ 600 sec?)
  }, error = function(e) {
    print(paste("Error occurred: Try to fetch again", e$message))
    classification(tax_id[tax_id_split[[tax_i]]], db = "ncbi") # Batch-size = 50 (~ 600 sec?)
    return(NULL)
  })
  # Combine results
  tax_list <- c(tax_list, list(result))
  # Message
  message(sprintf("\nCycle %s / %s finished!\n", tax_i, length(tax_id_split)))
}

# Unlist tax_list
tax_list2 <- list()
for (i in 1:length(tax_list)) {
  tax_list_tmp <- tax_list[[i]]
  tax_list2 <- c(tax_list2, tax_list_tmp)
}
na_id <- tax_list2 %>% map(function(x) all(is.na(x))) %>% unlist %>% which
names(tax_list2[na_id])

na_tax_tmp1 <- classification(tax_id[na_id][1:26], db = "ncbi") # Batch-size = 50 (~ 600 sec?)
na_tax_tmp2 <- classification(933952, db = "ncbi") # Batch-size = 50 (~ 600 sec?)
na_tax_tmp3 <- classification(tax_id[na_id][28:50], db = "ncbi") # Batch-size = 50 (~ 600 sec?)
na_all <- c(na_tax_tmp1, na_tax_tmp2, na_tax_tmp3)

# Save tax_list2
saveRDS(tax_list2,  sprintf("%s/tax_info.obj", output_folder))

# Add NA list
tax_list2[na_id] <- na_all
acc_tax_db$class <- NA
acc_tax_db$order <- NA
acc_tax_db$infraorder <- NA
acc_tax_db$family <- NA
acc_tax_db$genus <- NA

# Compile tax list
for (i in 1:length(tax_list2)) {
  rank1 <- which(tax_list2[[i]][,"rank"] == "class") %>% tax_list2[[i]][.,"name"] %>% .[1]
  rank2 <- which(tax_list2[[i]][,"rank"] == "order") %>% tax_list2[[i]][.,"name"] %>% .[1]
  rank3 <- which(tax_list2[[i]][,"rank"] == "infraorder") %>% tax_list2[[i]][.,"name"] %>% .[1]
  rank4 <- which(tax_list2[[i]][,"rank"] == "family") %>% tax_list2[[i]][.,"name"] %>% .[1]
  rank5 <- which(tax_list2[[i]][,"rank"] == "genus") %>% tax_list2[[i]][.,"name"] %>% .[1]
  acc_tax_db$class[i] <- rank1
  acc_tax_db$order[i] <- rank2
  acc_tax_db$infraorder[i] <- rank3
  acc_tax_db$family[i] <- rank4
  acc_tax_db$genus[i] <- rank5
}


# ------------------------------------------------ #
# Save results
# ------------------------------------------------ #
write.csv(acc_tax_db, sprintf("%s/acc_tax_db.csv", output_folder), row.names = F)

# Save sessioninfo
#macam::save_session_info()
