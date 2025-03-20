####
#### Visualize in silico PCR
#### 2025.01.24 Ushio
####

# Load libraries
library(tidyverse); packageVersion("tidyverse") # 1.3.2, 2023.7.18
library(cowplot); packageVersion("cowplot") # 1.1.1, 2023.7.18
library(viridis); packageVersion("viridis") # 0.6.3, 2023.7.19
library(cols4all); packageVersion("cols4all") # 0.4, 2023.7.19
theme_set(theme_cowplot())

# Create output directory
set.seed(1234)
(output_folder <- macam::outdir_create())


# ------------------------------------------------ #
# Load the previous results
# ------------------------------------------------ #
match_seq_all <- readRDS("10_SummarizeInSilicoPCR_EachOut/match_seq_all.obj")
match_sp_all <- readRDS("10_SummarizeInSilicoPCR_EachOut/match_sp_all.obj")
acc_taxdb <- read.csv("02_CompileTaxaOut/acc_tax_db_ed.csv")

# Modify the primer order
primer_order <- c("MiMammalF","MiMammalR","MarVer1F","MarVer1R","MarVer2F","MarVer2R",
                  "MarVer3F","MarVer3R","Ceto2F","Ceto2R","Riaz12SF","Riaz12SR","Mu31F",
                  "Dc320R","Dc321F","Dc494F","Dc495R","Mu642F","Mu643R","Dc671F",
                  "Dc1015R","Mu1021R","Dc1458F","Dc1465F","Dc1638R","Dc1646R","Mu2084F",
                  "Dc2173F","Mu2187F","Dc2385F","Dc2430F","Dc2438F","Dc2438R","Dc2505F",
                  "Mu2563R","Dc2580R","Mu9459F","Mu9822R","Mu9834R")
match_seq_all$primer_name <- factor(match_seq_all$primer_name, levels = primer_order)
match_sp_all$primer_name <- factor(match_sp_all$primer_name, levels = primer_order)

# Compile for visualization
match_seq_all$n_diff <- factor(match_seq_all$n_diff, levels = 4:0)
match_sp_all$n_diff <- factor(match_sp_all$n_diff, levels = 4:0)


# ------------------------------------------------ #
# Visualize results
# ------------------------------------------------ #
# Collect sep and species number
n_seq_cet <- acc_taxdb %>% filter(cat == "cetacea") %>% nrow
n_spp_cet <- acc_taxdb %>% filter(cat == "cetacea") %>% pull(species) %>% unique %>% length
n_seq_ver <- acc_taxdb %>% filter(cat == "vertebrate") %>% nrow
n_spp_ver <- acc_taxdb %>% filter(cat == "vertebrate") %>% pull(tax_id) %>% unique %>% length
n_seq_fis <- acc_taxdb %>% filter(cat == "fish") %>% nrow
n_spp_fis <- acc_taxdb %>% filter(cat == "fish") %>% pull(tax_id) %>% unique %>% length


# ------------------------------------------------ #
## N of sequences
# ------------------------------------------------ #
g1a <- match_seq_all %>% pivot_longer(cols = -c(n_diff, primer_name)) %>% 
  filter(name == "cetacea") %>% 
  ggplot(aes(x = primer_name, y = value, fill = n_diff)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = n_seq_cet, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = rev(c4a("tol.sunset", 5)), name = "primer\nmismatch") +
  xlab(NULL) + ylab("N of sequence amplified") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  panel_border() + ggtitle("N of sequences, Cetacea")

g1b <- match_seq_all %>% pivot_longer(cols = -c(n_diff, primer_name)) %>% 
  filter(name == "vertebrate") %>% 
  ggplot(aes(x = primer_name, y = value, fill = n_diff)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = n_seq_ver, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = rev(c4a("tol.sunset", 5)), name = "primer\nmismatch") +
  xlab(NULL) + ylab("N of sequence amplified") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  panel_border() + ggtitle("N of sequences, Non-fish vertebrates")

g1c <- match_seq_all %>% pivot_longer(cols = -c(n_diff, primer_name)) %>% 
  filter(name == "fish") %>% 
  ggplot(aes(x = primer_name, y = value, fill = n_diff)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = n_seq_fis, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = rev(c4a("tol.sunset", 5)), name = "primer\nmismatch") +
  xlab(NULL) + ylab("N of sequence amplified") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  panel_border() + ggtitle("N of sequences, Fish")


# ------------------------------------------------ #
## N of species
# ------------------------------------------------ #
g2a <- match_sp_all %>% pivot_longer(cols = -c(n_diff, primer_name)) %>% 
  filter(name == "cetacea") %>% 
  ggplot(aes(x = primer_name, y = value, fill = n_diff)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = n_spp_cet, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = rev(c4a("tol.sunset", 5)), name = "primer\nmismatch") +
  xlab(NULL) + ylab("N of species amplified") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  panel_border() + ggtitle("N of species, Cetacea")

g2b <- match_sp_all %>% pivot_longer(cols = -c(n_diff, primer_name)) %>% 
  filter(name == "vertebrate") %>% 
  ggplot(aes(x = primer_name, y = value, fill = n_diff)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = n_spp_ver, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = rev(c4a("tol.sunset", 5)), name = "primer\nmismatch") +
  xlab(NULL) + ylab("N of species amplified") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  panel_border() + ggtitle("N of species, Non-fish vertebrates")

g2c <- match_sp_all %>% pivot_longer(cols = -c(n_diff, primer_name)) %>% 
  filter(name == "fish") %>% 
  ggplot(aes(x = primer_name, y = value, fill = n_diff)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = n_spp_fis, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = rev(c4a("tol.sunset", 5)), name = "primer\nmismatch") +
  xlab(NULL) + ylab("N of species amplified") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  panel_border() + ggtitle("N of species, Fish")
 

# ------------------------------------------------ #
# Save results
# ------------------------------------------------ #
# Save figures
ggsave(file = sprintf("%s/nseq_amplified_allprimers.pdf", output_folder),
       plot = plot_grid(g1a, g1c, g1b, ncol = 1, align = "hv", axis = "lrbt"),
       width = 10, height = 14)
ggsave(file = sprintf("%s/nspp_amplified_allprimers.pdf", output_folder),
       plot = plot_grid(g2a, g2c, g2b, ncol = 1, align = "hv", axis = "lrbt"),
       width = 10, height = 14)

# Save workspace
#save(list = ls(all.names = TRUE), file = sprintf("%s/%s.RData", outdir, outdir))

# Save sessioninfo
macam::save_session_info()

