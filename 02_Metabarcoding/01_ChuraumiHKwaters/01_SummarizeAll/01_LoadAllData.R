####
#### Summarize pool test results
#### 2025.02.04, R4.3.2
####

# Load library and functions
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")

# Load custom R package
library(macam); packageVersion("macam") # 0.1.9

# Generate output folder
(output_folder <- macam::outdir_create())
dir.create("00_SessionInfo")


# --------------------------------------------- #
# Load data
# --------------------------------------------- #
# Primer set
# ps_all_p1 = all sequence reads, raw data
# ps_rar_p1 = after 99.9% coverage-based rarefaction
# ps_rel_p1 = the relative abundance of ps_rar_p1, cells < 0.5% were replaced with 0
# Primer set 1 (Cet1_Mu31F)
ps_all_p1 <- readRDS("../Cet1_Mu31F/05_SummarizeOut/ps_all.obj")
ps_rar_p1 <- readRDS("../Cet1_Mu31F/06_RarefactionOut/ps_rare.obj")
ps_rel_p1 <- readRDS("../Cet1_Mu31F/06_RarefactionOut/ps_rel.obj")
# Primer set 2 (Cet2_Dc671F)
ps_all_p2 <- readRDS("../Cet2_Dc671F/05_SummarizeOut/ps_all.obj")
ps_rar_p2 <- readRDS("../Cet2_Dc671F/06_RarefactionOut/ps_rare.obj")
ps_rel_p2 <- readRDS("../Cet2_Dc671F/06_RarefactionOut/ps_rel.obj")
# Primer set 3 (Cet3_Mu2084F)
ps_all_p3 <- readRDS("../Cet3_Mu2084F/05_SummarizeOut/ps_all.obj")
ps_rar_p3 <- readRDS("../Cet3_Mu2084F/06_RarefactionOut/ps_rare.obj")
ps_rel_p3 <- readRDS("../Cet3_Mu2084F/06_RarefactionOut/ps_rel.obj")
# Primer set 4 (Cet4_Mu9459F)
ps_all_p4 <- readRDS("../Cet4_Mu9459F/05_SummarizeOut/ps_all.obj")
ps_rar_p4 <- readRDS("../Cet4_Mu9459F/06_RarefactionOut/ps_rare.obj")
ps_rel_p4 <- readRDS("../Cet4_Mu9459F/06_RarefactionOut/ps_rel.obj")
# Primer set 5 (Cet5_MiMammal)
ps_all_p5 <- readRDS("../Cet5_MiMammal/05_SummarizeOut/ps_all.obj")
ps_rar_p5 <- readRDS("../Cet5_MiMammal/06_RarefactionOut/ps_rare.obj")
ps_rel_p5 <- readRDS("../Cet5_MiMammal/06_RarefactionOut/ps_rel.obj")
# Primer set 6 (Cet6_Ceto2)
ps_all_p6 <- readRDS("../Cet6_Ceto2/05_SummarizeOut/ps_all.obj")
ps_rar_p6 <- readRDS("../Cet6_Ceto2/06_RarefactionOut/ps_rare.obj")
ps_rel_p6 <- readRDS("../Cet6_Ceto2/06_RarefactionOut/ps_rel.obj")


# --------------------------------------------- #
# Compile sample and tax names
# --------------------------------------------- #
# Define new names
new_names_main <- factor(c("Lagoon_main", "Lagoon_north", "Lagoon_south", "Reproduction_pool",
                           "Main_pool", "Show_pool", "MedicalTrt_pool", "NC_ChField",
                           "SaiKung_Surface", "SaiKung_Subsurface", "NC_OPField",
                           "NC_DNAext1", "NC_DNAext2", "St_PM15", "NC_PCR1", "NC_PCR2"),
                         levels = c("Lagoon_main", "Lagoon_north", "Lagoon_south", "Reproduction_pool", "Main_pool", "Show_pool", "MedicalTrt_pool", 
                                    "SaiKung_Surface", "SaiKung_Subsurface", "St_PM15",
                                    "NC_ChField", "NC_OPField", "NC_DNAext1", "NC_DNAext2", "NC_PCR1", "NC_PCR2"))
new_names_si <- factor(c("Lagoon_main", "Lagoon_north", "Lagoon_south", "Reproduction_pool",
                           "Main_pool", "Show_pool", "MedicalTrt_pool", "SaiKung_Surface", "SaiKung_Subsurface", "St_PM15", "NC_PCR1", "NC_PCR2"),
                       levels = c("Lagoon_main", "Lagoon_north", "Lagoon_south", "Reproduction_pool", "Main_pool", "Show_pool", "MedicalTrt_pool", 
                                  "SaiKung_Surface", "SaiKung_Subsurface", "St_PM15",
                                  "NC_ChField", "NC_OPField", "NC_DNAext1", "NC_DNAext2", "NC_PCR1", "NC_PCR2"))
new_names1 <- new_names2 <- new_names3 <- new_names4 <- new_names_main
new_names5 <- new_names6 <- new_names_si
names(new_names1) <- sprintf("S0%02d", 1:16)
names(new_names2) <- sprintf("S0%02d", 17:32)
names(new_names3) <- sprintf("S0%02d", 33:48)
names(new_names4) <- sprintf("S0%02d", 49:64)
names(new_names5) <- sprintf("R0%02d", 1:12)
names(new_names6) <- sprintf("R0%02d", 13:24)
# Add new names to the phyloseq objects
## Primer set 1
sample_data(ps_all_p1)$new_name <- new_names1
sample_data(ps_rar_p1)$new_name <- sample_data(ps_rel_p1)$new_name <- new_names1[sample_data(ps_rar_p1)$Sample_Name2]
sample_data(ps_all_p1) <- sample_data(sample_data(ps_all_p1) %>% data.frame %>% select(-index, -I7_Index_ID, -index2, -I5_Index_ID))
sample_data(ps_rar_p1) <- sample_data(sample_data(ps_rar_p1) %>% data.frame %>% select(-index, -I7_Index_ID, -index2, -I5_Index_ID))
sample_data(ps_rel_p1) <- sample_data(sample_data(ps_rel_p1) %>% data.frame %>% select(-index, -I7_Index_ID, -index2, -I5_Index_ID))
## Primer set 2
sample_data(ps_all_p2)$new_name <- new_names2
sample_data(ps_rar_p2)$new_name <- sample_data(ps_rel_p2)$new_name <- new_names2[sample_data(ps_rar_p2)$Sample_Name2]
sample_data(ps_all_p2) <- sample_data(sample_data(ps_all_p2) %>% data.frame %>% select(-index, -I7_Index_ID, -index2, -I5_Index_ID))
sample_data(ps_rar_p2) <- sample_data(sample_data(ps_rar_p2) %>% data.frame %>% select(-index, -I7_Index_ID, -index2, -I5_Index_ID))
sample_data(ps_rel_p2) <- sample_data(sample_data(ps_rel_p2) %>% data.frame %>% select(-index, -I7_Index_ID, -index2, -I5_Index_ID))
## Primer set 3
sample_data(ps_all_p3)$new_name <- new_names3
sample_data(ps_rar_p3)$new_name <- sample_data(ps_rel_p3)$new_name <- new_names3[sample_data(ps_rar_p3)$Sample_Name2]
sample_data(ps_all_p3) <- sample_data(sample_data(ps_all_p3) %>% data.frame %>% select(-index, -I7_Index_ID, -index2, -I5_Index_ID))
sample_data(ps_rar_p3) <- sample_data(sample_data(ps_rar_p3) %>% data.frame %>% select(-index, -I7_Index_ID, -index2, -I5_Index_ID))
sample_data(ps_rel_p3) <- sample_data(sample_data(ps_rel_p3) %>% data.frame %>% select(-index, -I7_Index_ID, -index2, -I5_Index_ID))
## Primer set 4
sample_data(ps_all_p4)$new_name <- new_names4
sample_data(ps_rar_p4)$new_name <- sample_data(ps_rel_p4)$new_name <- new_names4[sample_data(ps_rar_p4)$Sample_Name2]
sample_data(ps_all_p4) <- sample_data(sample_data(ps_all_p4) %>% data.frame %>% select(-index, -I7_Index_ID, -index2, -I5_Index_ID))
sample_data(ps_rar_p4) <- sample_data(sample_data(ps_rar_p4) %>% data.frame %>% select(-index, -I7_Index_ID, -index2, -I5_Index_ID))
sample_data(ps_rel_p4) <- sample_data(sample_data(ps_rel_p4) %>% data.frame %>% select(-index, -I7_Index_ID, -index2, -I5_Index_ID))
## Primer set 5
sample_data(ps_all_p5)$new_name <- new_names5
sample_data(ps_rar_p5)$new_name <- sample_data(ps_rel_p5)$new_name <- new_names5[sample_data(ps_rar_p5)$Sample_Name2]
## Primer set 6
sample_data(ps_all_p6)$new_name <- new_names6
sample_data(ps_rar_p6)$new_name <- sample_data(ps_rel_p6)$new_name <- new_names6[sample_data(ps_rar_p6)$Sample_Name2]

# Rename tax names and remove singleton
## ps_all
ps_all_list <- list(ps_all_p1, ps_all_p2, ps_all_p3, ps_all_p4, ps_all_p5, ps_all_p6)
for (i in 1:6) {
  otu_tmp <- otu_table(ps_all_list[[i]]) %>% as.data.frame
  tax_tmp <- tax_table(ps_all_list[[i]]) %>% as.data.frame
  colnames(otu_tmp) <- sprintf("P%s_%s", i, colnames(otu_table(ps_all_list[[i]])))
  rownames(tax_tmp) <- sprintf("P%s_%s", i, rownames(tax_table(ps_all_list[[i]])))
  # Remove singletons
  otu_tmp[otu_tmp == 1] <- 0
  # Re-create phyloseq object
  ps_all_list[[i]] <- phyloseq(otu_table(otu_tmp, taxa_are_rows = FALSE),
                                   sample_data(ps_all_list[[i]]),
                                   tax_table(as.matrix(tax_tmp)))
}; rm(otu_tmp); rm(tax_tmp); rm(i)
p_all <- ps_all_list; rm(ps_all_list)

## ps_rar
ps_rar_all_list <- list(ps_rar_p1, ps_rar_p2, ps_rar_p3, ps_rar_p4, ps_rar_p5, ps_rar_p6)
for (i in 1:6) {
  otu_tmp <- otu_table(ps_rar_all_list[[i]]) %>% as.data.frame
  tax_tmp <- tax_table(ps_rar_all_list[[i]]) %>% as.data.frame
  colnames(otu_tmp) <- sprintf("P%s_%s", i, colnames(otu_table(ps_rar_all_list[[i]])))
  rownames(tax_tmp) <- sprintf("P%s_%s", i, rownames(tax_table(ps_rar_all_list[[i]])))
  # Re-create phyloseq object
  ps_rar_all_list[[i]] <- phyloseq(otu_table(otu_tmp, taxa_are_rows = FALSE),
                                   sample_data(ps_rar_all_list[[i]]),
                                   tax_table(as.matrix(tax_tmp)))
}; rm(otu_tmp); rm(tax_tmp); rm(i)
prar <- ps_rar_all_list; rm(ps_rar_all_list)

## ps_rel
ps_rel_all_list <- list(ps_rel_p1, ps_rel_p2, ps_rel_p3, ps_rel_p4, ps_rel_p5, ps_rel_p6)
for (i in 1:6) {
  otu_tmp <- otu_table(ps_rel_all_list[[i]]) %>% as.data.frame
  tax_tmp <- tax_table(ps_rel_all_list[[i]]) %>% as.data.frame
  colnames(otu_tmp) <- sprintf("P%s_%s", i, colnames(otu_table(ps_rel_all_list[[i]])))
  rownames(tax_tmp) <- sprintf("P%s_%s", i, rownames(tax_table(ps_rel_all_list[[i]])))
  # Re-create phyloseq object
  ps_rel_all_list[[i]] <- phyloseq(otu_table(otu_tmp, taxa_are_rows = FALSE),
                                   sample_data(ps_rel_all_list[[i]]),
                                   tax_table(as.matrix(tax_tmp)))
}; rm(otu_tmp); rm(tax_tmp); rm(i)
prel <- ps_rel_all_list; rm(ps_rel_all_list)


# --------------------------------------------- #
# Combine all the phyloseq objects
# --------------------------------------------- #
ps_all_raw <- merge_phyloseq(p_all[[1]], p_all[[2]], p_all[[3]], p_all[[4]], p_all[[5]], p_all[[6]])
ps_all_reads <- merge_phyloseq(prar[[1]], prar[[2]], prar[[3]], prar[[4]], prar[[5]], prar[[6]])
ps_all_rel <- merge_phyloseq(prel[[1]], prel[[2]], prel[[3]], prel[[4]], prel[[5]], prel[[6]]) %>% 
  transform_sample_counts(function(x) x/sum(x))

# Compile data
sample_data(ps_all_raw)$Primer <- factor(sample_data(ps_all_raw)$Primer,
                                         levels = c("Mu31F_Dc320R","Dc671F_Dc1015R",
                                                    "Mu2084F_Dc2438R","Mu9459F_Mu9822R",
                                                    "MiMammalF-R","Ceto2F-R"))
sample_data(ps_all_reads)$Primer <- factor(sample_data(ps_all_reads)$Primer,
                                         levels = c("Mu31F_Dc320R","Dc671F_Dc1015R",
                                                    "Mu2084F_Dc2438R","Mu9459F_Mu9822R",
                                                    "MiMammalF-R","Ceto2F-R"))
sample_data(ps_all_rel)$Primer <- factor(sample_data(ps_all_rel)$Primer,
                                     levels = c("Mu31F_Dc320R","Dc671F_Dc1015R",
                                                "Mu2084F_Dc2438R","Mu9459F_Mu9822R",
                                                "MiMammalF-R","Ceto2F-R"))


# --------------------------------------------- #
# Save results and figures
# --------------------------------------------- #
# Save figure
saveRDS(ps_all_raw, "01_LoadAllDataOut/ps_all_raw.obj")
saveRDS(ps_all_reads, "01_LoadAllDataOut/ps_all_rarefy_reads.odj")
saveRDS(ps_all_rel, "01_LoadAllDataOut/ps_all_rel.odj")

# Save tables
write.csv(sample_data(ps_all_raw) %>% data.frame, "01_LoadAllDataOut/ps_all_raw_smd.csv")
write.csv(otu_table(ps_all_raw) %>% data.frame, "01_LoadAllDataOut/ps_all_raw_otu.csv")
write.csv(tax_table(ps_all_raw) %>% data.frame, "01_LoadAllDataOut/ps_all_raw_tax.csv")

# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("00_SessionInfo/", output_folder, "_", substr(Sys.time(), 1, 10), ".txt"))

save(list = ls(all.names = TRUE),
     file = paste0(output_folder, "/", output_folder, ".RData"))

