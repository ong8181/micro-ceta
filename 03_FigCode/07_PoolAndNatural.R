####
#### Summarize all data
#### 2025.02.20, Ushio
####

# Load library and functions
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(cowplot); packageVersion("cowplot")
library(cols4all); packageVersion("cols4all")
library(ggtext); packageVersion("ggtext")
theme_set(theme_cowplot())


# --------------------------------------------- #
# Load data
# --------------------------------------------- #
ps_all <- readRDS("../02_Metabarcoding/01_ChuraumiHKwaters/01_SummarizeAll/01_LoadAllDataOut/ps_all_raw.obj")
ps_rar <- readRDS("../02_Metabarcoding/01_ChuraumiHKwaters/01_SummarizeAll/01_LoadAllDataOut/ps_all_rarefy_reads.odj")
ps_rel <- readRDS("../02_Metabarcoding/01_ChuraumiHKwaters/01_SummarizeAll/01_LoadAllDataOut/ps_all_rel.odj")

# Remove samples if the total original reads were less than 10
keep_samples <- names(which(sample_sums(ps_all) >= 10))
ps_rel <- ps_rel %>% prune_samples(keep_samples, .)


# --------------------------------------------- #
# Visualize overall composition
# --------------------------------------------- #
# Combine tax names
tax_rank_list <- (tax_table(ps_rel) %>% colnames)[2:17]
target_taxa <- "infraorder"
ps_sub <- macamseq::taxa_name_bundle(ps_rel, target_taxa, top_taxa_n = 10, taxa_rank_list = tax_rank_list)
ps_mt1 <- speedyseq::psmelt(ps_sub)

# Edit taxa names
ps_mt1$bundled_tax[ps_mt1$superkingdom == "Bacteria"] <- "Bacteria"
ps_mt1$bundled_tax[ps_mt1$species == "Homo sapiens"] <- "_Homo sapiens_"
ps_mt1$bundled_tax[ps_mt1$family == "Hominidae"] <- "_Homo sapiens_"
ps_mt1$bundled_tax[ps_mt1$assignSp_Species == "sapiens"] <- "_Homo sapiens_"
ps_mt1$bundled_tax[ps_mt1$species == "Canis lupus"] <- "_Canis lupus_"
ps_mt1$bundled_tax[ps_mt1$species == "Sus scrofa"] <- "_Sus scrofa_"
ps_mt1$bundled_tax[ps_mt1$species == "Bos taurus"] <- "_Bos taurus_"
ps_mt1$bundled_tax[ps_mt1$bundled_tax == "Undetermined" & ps_mt1$superkingdom == ""] <- "Undetermined"
ps_mt1$bundled_tax[ps_mt1$bundled_tax == "Undetermined" & ps_mt1$class == "Mammalia"] <- "Others"
ps_mt1$bundled_tax[ps_mt1$bundled_tax == "Undetermined" & ps_mt1$superkingdom == "Eukaryota"] <- "Other Eukaryota"
# Edit order
ps_mt1$bundled_tax <- factor(ps_mt1$bundled_tax,
                             levels = c("Cetacea","Bacteria","_Homo sapiens_","_Sus scrofa_","_Bos taurus_",
                                        "_Canis lupus_", "Other Eukaryota","Others","Undetermined"))

ps_mt2 <- ps_mt1 %>% group_by_at(c("new_name", "Primer", "bundled_tax")) %>%
  summarize(sequence_reads = sum(Abundance))

pool_names <- c("Lagoon main", "Lagoon north", "Lagoon south", "Main pool",
                    "Medical treatment pool", "Reproduction pool",
                    "Show pool", "Natural sample 1", "Natural sample 2", "Natural sample 3")
ps_mt2$new_name2 <- as.character(ps_mt2$new_name)
ps_mt2$new_name2[ps_mt2$new_name2 == "St_PM15"] <- pool_names[10]
ps_mt2$new_name2[ps_mt2$new_name2 == "SaiKung_Subsurface"] <- pool_names[9]
ps_mt2$new_name2[ps_mt2$new_name2 == "SaiKung_Surface"] <- pool_names[8]
ps_mt2$new_name2[ps_mt2$new_name2 == "Show_pool"] <- pool_names[7]
ps_mt2$new_name2[ps_mt2$new_name2 == "Reproduction_pool"] <- pool_names[6]
ps_mt2$new_name2[ps_mt2$new_name2 == "MedicalTrt_pool"] <- pool_names[5]
ps_mt2$new_name2[ps_mt2$new_name2 == "Main_pool"] <- pool_names[4]
ps_mt2$new_name2[ps_mt2$new_name2 == "Lagoon_south"] <- pool_names[3]
ps_mt2$new_name2[ps_mt2$new_name2 == "Lagoon_north"] <- pool_names[2]
ps_mt2$new_name2[ps_mt2$new_name2 == "Lagoon_main"] <- pool_names[1]
ps_mt2$new_name2 <- factor(ps_mt2$new_name2, levels = rev(pool_names))

# Rename primers
ps_mt2$Primer <- ps_mt2$Primer %>% str_replace_all("_", "/")
ps_mt2$Primer <- ps_mt2$Primer %>% str_replace_all("F-R", "")
ps_mt2$Primer <- factor(ps_mt2$Primer,
                        levels = c("Mu31F/Dc320R","Dc671F/Dc1015R",
                                   "Mu2084F/Dc2438R","Mu9459F/Mu9822R",
                                   "MiMammal","Ceto2"))

# Visualize data
f1 <- ps_mt2 %>% 
  ggplot(aes_(x = as.name("new_name2"), y = as.name("sequence_reads"),
              fill = as.name("bundled_tax"))) +
  geom_bar(stat = "identity", colour = NA) +
  facet_wrap(.~ Primer) +
  scale_fill_manual(values = c(c4a("brewer.paired")[1:8],"gray80","gray50")) +
  coord_flip() +
  theme(legend.position = "bottom", axis.text.y = ggtext::element_markdown(),
        legend.text = ggtext::element_markdown()) + 
  labs(title = "99.9% coverage-based rarefaction and remove singleton",
       y = "Relative abundance", x = NULL) +
  NULL

f2 <- ps_mt2 %>% 
  filter(Primer != "MiMammal" & Primer != "Ceto2") %>% 
  ggplot(aes_(x = as.name("new_name2"), y = as.name("sequence_reads"),
                         fill = as.name("bundled_tax"))) +
  geom_bar(stat = "identity", colour = NA) +
  facet_wrap(.~ Primer) +
  scale_fill_manual(values = c(c4a("brewer.paired")[1:8],"gray80","gray50")) +
  coord_flip() +
  theme(legend.position = "bottom", axis.text.y = ggtext::element_markdown(),
        legend.text = ggtext::element_markdown()) + 
  labs(title = "99.9% coverage-based rarefaction and remove singleton",
       y = "Relative abundance", x = NULL) +
  NULL

f3 <- ps_mt2 %>% 
  filter(Primer == "MiMammal" | Primer == "Ceto2") %>% 
  ggplot(aes_(x = as.name("new_name2"), y = as.name("sequence_reads"),
                          fill = as.name("bundled_tax"))) +
  geom_bar(stat = "identity", colour = NA) +
  facet_wrap(.~ Primer) +
  scale_fill_manual(values = c(c4a("brewer.paired")[1:8],"gray80","gray50")) +
  coord_flip() +
  theme(legend.position = "bottom", axis.text.y = ggtext::element_markdown(),
        legend.text = ggtext::element_markdown()) + 
  labs(title = "99.9% coverage-based rarefaction and remove singleton",
       y = "Relative abundance", x = NULL) +
  NULL


# --------------------------------------------- #
# Save results and figures
# --------------------------------------------- #
# Save figure
saveRDS(list(f1, f2, f3), "data_robj/PoolNaturalTest.obj")

# Save sessioninfo
macam::save_session_info()

