####
#### Visualize in silico PCR
#### 2025.02.20, Ushio
####

# Load libraries
library(tidyverse); packageVersion("tidyverse")
library(cowplot); packageVersion("cowplot")
library(viridis); packageVersion("viridis")
library(cols4all); packageVersion("cols4all")
library(ggtext); packageVersion("ggtext")
theme_set(theme_cowplot())

# Create output directory
set.seed(1234)


# ------------------------------------------------ #
# Load the previous results
# ------------------------------------------------ #
amplicon_len_list <- readRDS("../01_InSilicoPCR/03_SummarizeInSilicoPCROut/amplicon_len_list.obj")
amplicon_len_others_list <- readRDS("../01_InSilicoPCR/03_SummarizeInSilicoPCROut/amplicon_len_others_list.obj")
match_seq_all <- readRDS("../01_InSilicoPCR/03_SummarizeInSilicoPCROut/match_seq_all.obj")
match_sp_all <- readRDS("../01_InSilicoPCR/03_SummarizeInSilicoPCROut/match_sp_all.obj")
acc_taxdb <- read.csv("../01_InSilicoPCR/02_CompileTaxaOut/acc_tax_db_ed.csv")

# Modify the primer order
## Replace "_" with "/"
match_seq_all$primer_name <- match_seq_all$primer_name %>% str_replace_all("_", "/")
match_sp_all$primer_name <- match_sp_all$primer_name %>% str_replace_all("_", "/")
match_seq_all$primer_name <- match_seq_all$primer_name %>% str_replace_all("Mu31F/Dc320R", "**µCeta | Mu31F/Dc320R**")
match_sp_all$primer_name <- match_sp_all$primer_name %>% str_replace_all("Mu31F/Dc320R", "**µCeta | Mu31F/Dc320R**")
match_seq_all$primer_name <- match_seq_all$primer_name %>% str_replace_all("Dc671F/Dc1015R", "**Dc671F/Dc1015R**")
match_sp_all$primer_name <- match_sp_all$primer_name %>% str_replace_all("Dc671F/Dc1015R", "**Dc671F/Dc1015R**")
match_seq_all$primer_name <- match_seq_all$primer_name %>% str_replace_all("Mu2084F/Dc2438R", "**Mu2084F/Dc2438R**")
match_sp_all$primer_name <- match_sp_all$primer_name %>% str_replace_all("Mu2084F/Dc2438R", "**Mu2084F/Dc2438R**")
match_seq_all$primer_name <- match_seq_all$primer_name %>% str_replace_all("Mu9459F/Mu9822R", "**Mu9459F/Mu9822R**")
match_sp_all$primer_name <- match_sp_all$primer_name %>% str_replace_all("Mu9459F/Mu9822R", "**Mu9459F/Mu9822R**")

## Define the primer order
primer_order <- c("MiMammal", "MarVer1", "MarVer2", "MarVer3", "Ceto2", "Riaz12S",
                  "**µCeta | Mu31F/Dc320R**", "Mu31F/MiMammalR", "Dc321F/Dc495R", "Dc494F/Mu643R",
                  "Mu642F/Dc1015R", "Mu642F/Mu1021R", "**Dc671F/Dc1015R**", "Dc671F/Mu1021R",
                  "Dc1458/1638", "Dc1465F/Dc1646R", "**Mu2084F/Dc2438R**", "Dc2173/2580",
                  "Mu2187F/Dc2438R", "Mu2187F/Mu2563R", "Dc2385/2580", "Dc2430/2580",
                  "Dc2438F/Mu2563R", "Dc2505/2580", "**Mu9459F/Mu9822R**", "Mu9459F/Mu9834R")
match_seq_all$primer_name <- factor(match_seq_all$primer_name, levels = primer_order)
match_sp_all$primer_name <- factor(match_sp_all$primer_name, levels = primer_order)

# Compile for visualization
match_seq_all$n_diff <- factor(match_seq_all$n_diff, levels = 6:0)
match_sp_all$n_diff <- factor(match_sp_all$n_diff, levels = 6:0)


# ------------------------------------------------ #
# Visualize results
# ------------------------------------------------ #
# Collect sep and species number
n_seq_cet <- acc_taxdb %>% filter(cat == "cetacea") %>% nrow
n_spp_cet <- acc_taxdb %>% filter(cat == "cetacea") %>% pull(species) %>% unique %>% length
n_seq_ver <- acc_taxdb %>% filter(cat == "vertebrate") %>% nrow
n_spp_ver <- acc_taxdb %>% filter(cat == "vertebrate") %>% pull(species) %>% unique %>% length
n_seq_fis <- acc_taxdb %>% filter(cat == "fish") %>% nrow
n_spp_fis <- acc_taxdb %>% filter(cat == "fish") %>% pull(species) %>% unique %>% length

acc_taxdb %>% filter(rep_tax == "mammal") %>% pull(tax_id) %>% unique %>% length
acc_taxdb %>% filter(rep_tax == "bird") %>% pull(tax_id) %>% unique %>% length
acc_taxdb %>% filter(rep_tax == "fish") %>% pull(tax_id) %>% unique %>% length
acc_taxdb %>% filter(rep_tax == "reptile") %>% pull(tax_id) %>% unique %>% length
acc_taxdb %>% filter(rep_tax == "amphibian") %>% pull(tax_id) %>% unique %>% length
acc_taxdb %>% filter(rep_tax != "mammal" &
                       rep_tax != "bird" &
                       rep_tax != "fish" &
                       rep_tax != "reptile" &
                       rep_tax != "amphibian") %>% pull(tax_id) %>% unique %>% length
length(acc_taxdb %>% pull(tax_id) %>% unique)


# ------------------------------------------------ #
## N of sequences
# ------------------------------------------------ #
g1 <- match_seq_all %>% pivot_longer(cols = -c(n_diff, primer_name)) %>% 
  filter(name == "cetacea") %>% 
  ggplot(aes(x = primer_name, y = value, fill = n_diff)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = n_seq_cet, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = rev(c4a("tol.sunset")), name = "primer\nmismatch") +
  xlab(NULL) + ylab("N of sequence amplified") +
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5)) +
  panel_border() + ggtitle("N of sequences amplified, Cetacea")

g2 <- match_seq_all %>% pivot_longer(cols = -c(n_diff, primer_name)) %>% 
  filter(name == "fish") %>% 
  ggplot(aes(x = primer_name, y = value, fill = n_diff)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = n_seq_fis, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = rev(c4a("tol.sunset")), name = "primer\nmismatch") +
  xlab(NULL) + ylab("N of sequence amplified") +
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5)) +
  panel_border() + ggtitle("N of sequences amplified, Fish")

g3 <- match_seq_all %>% pivot_longer(cols = -c(n_diff, primer_name)) %>% 
  filter(name == "vertebrate") %>% 
  ggplot(aes(x = primer_name, y = value, fill = n_diff)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = n_seq_ver, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = rev(c4a("tol.sunset")), name = "primer\nmismatch") +
  xlab(NULL) + ylab("N of sequence amplified") +
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5)) +
  panel_border() + ggtitle("N of sequences amplified, Non-fish vertebrates")


# ------------------------------------------------ #
# N of species
# ------------------------------------------------ #
g4 <- match_sp_all %>% pivot_longer(cols = -c(n_diff, primer_name)) %>% 
  filter(name == "cetacea") %>% 
  ggplot(aes(x = primer_name, y = value, fill = n_diff)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = n_spp_cet, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = rev(c4a("tol.sunset")), name = "primer\nmismatch") +
  xlab(NULL) + ylab("N of species amplified") +
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5)) +
  panel_border() + ggtitle("N of species amplified, Cetacea")

g5 <- match_sp_all %>% pivot_longer(cols = -c(n_diff, primer_name)) %>% 
  filter(name == "fish") %>% 
  ggplot(aes(x = primer_name, y = value, fill = n_diff)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = n_spp_fis, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = rev(c4a("tol.sunset")), name = "primer\nmismatch") +
  xlab(NULL) + ylab("N of species amplified") +
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5)) +
  panel_border() + ggtitle("N of species amplified, Fish")

g6 <- match_sp_all %>% pivot_longer(cols = -c(n_diff, primer_name)) %>% 
  filter(name == "vertebrate") %>% 
  ggplot(aes(x = primer_name, y = value, fill = n_diff)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = n_spp_ver, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = rev(c4a("tol.sunset")), name = "primer\nmismatch") +
  xlab(NULL) + ylab("N of species amplified") +
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5)) +
  panel_border() + ggtitle("N of species amplified, Non-fish vertebrates")


# ------------------------------------------------ #
# Amplicon length (including primers)
# ------------------------------------------------ #
amplicon_long <- data.frame()
for (i in 1:length(amplicon_len_list)) {
  amplicon_tmp <- data.frame(amplicon_length = amplicon_len_list[[i]],
                             primer = names(amplicon_len_list)[i])
  amplicon_long <- rbind(amplicon_long, amplicon_tmp)
}
amplicon_others_long <- data.frame()
for (i in 1:length(amplicon_len_others_list)) {
  amplicon_others_tmp <- data.frame(amplicon_length = amplicon_len_others_list[[i]],
                             primer = names(amplicon_len_others_list)[i])
  amplicon_others_long <- rbind(amplicon_others_long, amplicon_others_tmp)
}
## Delete temporal object
rm(amplicon_tmp); rm(amplicon_others_tmp)
## Replace "_" with "/"
amplicon_long$primer <- amplicon_long$primer %>% str_replace_all("_", "/")
amplicon_long$primer <- amplicon_long$primer %>% str_replace_all("Mu31F/Dc320R", "**µCeta | Mu31F/Dc320R**")
amplicon_long$primer <- amplicon_long$primer %>% str_replace_all("Dc671F/Dc1015R", "**Dc671F/Dc1015R**")
amplicon_long$primer <- amplicon_long$primer %>% str_replace_all("Mu2084F/Dc2438R", "**Mu2084F/Dc2438R**")
amplicon_long$primer <- amplicon_long$primer %>% str_replace_all("Mu9459F/Mu9822R", "**Mu9459F/Mu9822R**")
amplicon_others_long$primer <- amplicon_others_long$primer %>% str_replace_all("_", "/")
amplicon_others_long$primer <- amplicon_others_long$primer %>% str_replace_all("Mu31F/Dc320R", "**µCeta | Mu31F/Dc320R**")
amplicon_others_long$primer <- amplicon_others_long$primer %>% str_replace_all("Dc671F/Dc1015R", "**Dc671F/Dc1015R**")
amplicon_others_long$primer <- amplicon_others_long$primer %>% str_replace_all("Mu2084F/Dc2438R", "**Mu2084F/Dc2438R**")
amplicon_others_long$primer <- amplicon_others_long$primer %>% str_replace_all("Mu9459F/Mu9822R", "**Mu9459F/Mu9822R**")
## Re-order primers
amplicon_long$primer <- factor(amplicon_long$primer, levels = primer_order)
amplicon_others_long$primer <- factor(amplicon_others_long$primer, levels = primer_order)
## Visualize
b1 <- amplicon_long %>%
  ggplot(aes(x = primer, y = amplicon_len)) +
  geom_jitter(alpha = 0.1, height = 0, width = 0.2) +
  geom_boxplot(outlier.size = 0, outlier.shape = NA, size = 0.2, alpha = 0.7) +
  xlab(NULL) + ylab("Amplicon length (bp)") + ggtitle("Simulated amplicon length (Cetacea)") +
  scale_y_continuous(breaks = seq(0, 2000, by = 100)) + coord_cartesian(ylim = c(0,1000)) +
  theme(axis.text.x = element_markdown(vjust = 0.5, hjust = 1, angle = 90),
        panel.grid = element_line(linewidth = 0.05))
b2 <- amplicon_others_long %>%
  ggplot(aes(x = primer, y = amplicon_len)) +
  geom_jitter(alpha = 0.1, height = 0, width = 0.2) +
  geom_boxplot(outlier.size = 0, outlier.shape = NA, size = 0.2, alpha = 0.9) +
  xlab(NULL) + ylab("Amplicon length (bp)") + ggtitle("Simulated amplicon length (Others)") +
  scale_y_continuous(breaks = seq(0, 2000, by = 100)) + coord_cartesian(ylim = c(0,1000)) +
  theme(axis.text.x = element_markdown(vjust = 0.5, hjust = 1, angle = 90),
        panel.grid = element_line(linewidth = 0.05))


# ------------------------------------------------ #
# Save results
# ------------------------------------------------ #
# Save figures
saveRDS(list(g1, g2, g3), "data_robj/InSilicoPCR_nseq.obj")
saveRDS(list(g4, g5, g6), "data_robj/InSilicoPCR_nspp.obj")
saveRDS(list(b1, b2), "data_robj/InSilicoPCR_amplicon_length.obj")

# Amplicon length
ggsave("data_img/InSilicoPCR_amplicon_length_cetacea.jpg", plot = b1, width = 12, height = 7)
ggsave("data_img/InSilicoPCR_amplicon_length_others.jpg", plot = b2, width = 12, height = 7)

# Save sessioninfo
macam::save_session_info(create_session_info_dir = TRUE)
