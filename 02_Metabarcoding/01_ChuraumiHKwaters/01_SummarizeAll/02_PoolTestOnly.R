####
#### Summarize pool test results
#### 2025.02.04 Ushio, R4.3.2
####

# Load library and functions
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(cowplot); packageVersion("cowplot")
library(cols4all); packageVersion("cols4all")
library(ggtext); packageVersion("ggtext")
theme_set(theme_cowplot())

# Load custom R package
library(macam); packageVersion("macam") # 0.1.9

# Generate output folder
(output_folder <- macam::outdir_create())


# --------------------------------------------- #
# Load data
# --------------------------------------------- #
ps_all <- readRDS("01_LoadAllDataOut/ps_all_raw.obj")
ps_rar <- readRDS("01_LoadAllDataOut/ps_all_rarefy_reads.odj")
ps_rel <- readRDS("01_LoadAllDataOut/ps_all_rel.odj")

# Remove samples if the total original reads were less than 10
keep_samples <- names(which(sample_sums(ps_all) >= 10))
ps_rel <- ps_rel %>% prune_samples(keep_samples, .)


# --------------------------------------------- #
# Select pool test samples only
# --------------------------------------------- #
chura_test <- sample_data(ps_rel)$new_name %>% unique() %>% .[1:7]
ps_chr <- prune_samples(sample_data(ps_rel)$new_name %in% chura_test, ps_rel)


# --------------------------------------------- #
# Select Cetacean ASVs only and visualize
# --------------------------------------------- #
ps_chr <- ps_chr %>% subset_taxa(infraorder == "Cetacea") %>%
  transform_sample_counts(function(x) x/sum(x)) %>% 
  prune_taxa(taxa_sums(.) > 0, .)

# Bundle taxa
tax_rank_list <- (tax_table(ps_chr) %>% colnames)[2:17]
target_taxa <- "species"
ps_sub <- macamseq::taxa_name_bundle(ps_chr, target_taxa, top_taxa_n = 10, taxa_rank_list = tax_rank_list)
ps_mt1 <- speedyseq::psmelt(ps_sub)
ps_mt2 <- ps_mt1 %>% group_by_at(c("new_name", "Primer", "bundled_tax")) %>%
  summarize(sequence_reads = sum(Abundance))


f1 <- ggplot(ps_mt2, aes_(x = as.name("new_name"), y = as.name("sequence_reads"),
                         fill = as.name("bundled_tax"))) +
  geom_bar(stat = "identity", colour = NA) +
  facet_wrap(.~ Primer) +
  scale_fill_manual(values = c(c4a("brewer.paired"))) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 8)) + 
  labs(title = "99.9% coverage-based rarefaction",
       y = "Relative abundance", x = NULL) +
  NULL


# --------------------------------------------- #
# Save cetacean ASVs for manual editing
# --------------------------------------------- #
write.csv(otu_table(ps_chr) %>% as.data.frame, "02_PoolTestOnlyOut/otu_table.csv")
write.csv(data.frame(sample_data(ps_chr)), "02_PoolTestOnlyOut/sample_data.csv")
write.csv(tax_table(ps_chr) %>% as.data.frame, "02_PoolTestOnlyOut/tax_table.csv")

# Load manually edited file
sample_sheet_ed <- read.csv("02_PoolTestOnlyOut/data_manual_ed/sample_data_ed.csv", row.names = 1)
otu_sheet_ed <- read.csv("02_PoolTestOnlyOut/data_manual_ed/otu_table_ed.csv", row.names = 1)
tax_sheet_ed <- read.csv("02_PoolTestOnlyOut/data_manual_ed/tax_table_ed.csv", row.names = 1)

# Check structure
dim(sample_sheet_ed); dim(sample_data(ps_chr))
colnames(sample_sheet_ed); colnames(sample_data(ps_chr))
all(rownames(sample_sheet_ed) == rownames(sample_data(ps_chr)))
dim(otu_sheet_ed); dim(otu_table(ps_chr))
all(otu_sheet_ed == otu_table(ps_chr))
dim(tax_sheet_ed); dim(tax_table(ps_chr))
all(tax_sheet_ed$seq == (data.frame(tax_table(ps_chr)) %>% pull(seq)))

# Import to phyloseq
ps_chr_ed <- phyloseq(otu_table(otu_sheet_ed, taxa_are_rows = FALSE),
                   sample_data(sample_sheet_ed),
                   tax_table(as.matrix(tax_sheet_ed))) %>% 
  transform_sample_counts(function(x) x/sum(x))
sample_data(ps_chr_ed)$Primer <- factor(sample_data(ps_chr_ed)$Primer,
                                        levels = c("Mu31F_Dc320R","Dc671F_Dc1015R",
                                                   "Mu2084F_Dc2438R","Mu9459F_Mu9822R",
                                                   "MiMammalF-R","Ceto2F-R"))
# Edit taxa
ps_mt3 <- speedyseq::psmelt(ps_chr_ed)
spp_names <- c("_Pseudorca crassidens_", "_Tursiops truncatus_", "_Tursiops aduncus_", "_Steno bredanensis_", "_Feresa attenuata_", 
               "Unidentified Delphinidae", "Unidentified Cetacea", "Unreared species")
ps_mt3$manual_species3 <- as.character(ps_mt3$manual_species2)
ps_mt3$manual_species3[ps_mt3$manual_species3 == "Pseudorca crassidens"] <- spp_names[1]
ps_mt3$manual_species3[ps_mt3$manual_species3 == "Tursiops truncatus"] <- spp_names[2]
ps_mt3$manual_species3[ps_mt3$manual_species3 == "Tursiops aduncus"] <- spp_names[3]
ps_mt3$manual_species3[ps_mt3$manual_species3 == "Steno bredanensis"] <- spp_names[4]
ps_mt3$manual_species3[ps_mt3$manual_species3 == "Feresa attenuata"] <- spp_names[5]
ps_mt3$manual_species3 <- factor(ps_mt3$manual_species3, levels = spp_names)

pool_spp_names <- c("Lagoon main (_P. crassidens_, _T. truncatus_)", "Lagoon north (_T. aduncus_)",
                    "Lagoon south (_P. crassidens_, _T. truncatus_)", "Main pool (_P. crassidens_, _T. truncatus_, _T. aduncus_)",
                    "Medical treatment (_S. bredanensis_)", "Reproduction (_P. crassidens_, _T. truncatus_ x _aduncus_)",
                    "Show pool (_T. aduncus_, _F. attenuata_)")
ps_mt3$new_name2 <- ps_mt3$new_name
ps_mt3$new_name2[ps_mt3$new_name2 == "Show_pool"] <- pool_spp_names[7]
ps_mt3$new_name2[ps_mt3$new_name2 == "Reproduction_pool"] <- pool_spp_names[6]
ps_mt3$new_name2[ps_mt3$new_name2 == "MedicalTrt_pool"] <- pool_spp_names[5]
ps_mt3$new_name2[ps_mt3$new_name2 == "Main_pool"] <- pool_spp_names[4]
ps_mt3$new_name2[ps_mt3$new_name2 == "Lagoon_south"] <- pool_spp_names[3]
ps_mt3$new_name2[ps_mt3$new_name2 == "Lagoon_north"] <- pool_spp_names[2]
ps_mt3$new_name2[ps_mt3$new_name2 == "Lagoon_main"] <- pool_spp_names[1]
ps_mt3$new_name2 <- factor(ps_mt3$new_name2, levels = rev(pool_spp_names))

ps_mt4 <- ps_mt3 %>% group_by_at(c("new_name2", "Primer", "manual_species3")) %>%
  summarize(sequence_reads = sum(Abundance))

f2 <- ggplot(ps_mt4, aes_(x = as.name("new_name2"), y = as.name("sequence_reads"),
                          fill = as.name("manual_species3"))) +
  geom_bar(stat = "identity", colour = NA) +
  facet_wrap(.~ Primer) +
  scale_fill_manual(name = "Detectd eDNA", values = c(c4a("brewer.paired")[1:5], "gray80", "gray40", "red3")) +
  coord_flip() +
  theme(legend.position = "bottom", axis.text.y = ggtext::element_markdown(),
        legend.text = ggtext::element_markdown()) + 
  labs(title = "99.9% coverage-based rarefaction and remove singlton and contaminations",
       y = "Relative abundance", x = NULL) +
  NULL


# --------------------------------------------- #
# Save results and figures
# --------------------------------------------- #
# Save figure
ggsave(sprintf("%s/Pool_Cetacean.pdf", output_folder), f1, width = 10, height = 6)
ggsave(sprintf("%s/Pool_Cetacean_ed.pdf", output_folder), f2, width = 14, height = 8)

# Save session info
writeLines(capture.output(sessionInfo()),
           paste0("00_SessionInfo/", output_folder, "_", substr(Sys.time(), 1, 10), ".txt"))

save(list = ls(all.names = TRUE),
     file = paste0(output_folder, "/", output_folder, ".RData"))

