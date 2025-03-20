####
#### Visualize in silico PCR
#### 2023.07.21 Ushio
####

# Load libraries
library(tidyverse); packageVersion("tidyverse") # 1.3.2, 2023.7.18
library(cowplot); packageVersion("cowplot") # 1.1.1, 2023.7.18
library(viridis); packageVersion("viridis") # 0.6.3, 2023.7.19
library(cols4all); packageVersion("cols4all") # 0.4, 2023.7.19
#library(macam); packageVersion("macam") # 0.1.4, 2023.7.18
theme_set(theme_cowplot())

# Create output directory
set.seed(1234)
(output_folder <- macam::outdir_create())


# ------------------------------------------------ #
# Load the previous results
# ------------------------------------------------ #
amplicon_len_list <- readRDS("03_SummarizeInSilicoPCROut/amplicon_len_list.obj")
amplicon_len_others_list <- readRDS("03_SummarizeInSilicoPCROut/amplicon_len_others_list.obj")
match_seq_all <- readRDS("03_SummarizeInSilicoPCROut/match_seq_all.obj")
match_sp_all <- readRDS("03_SummarizeInSilicoPCROut/match_sp_all.obj")
acc_taxdb <- read.csv("02_CompileTaxaOut/acc_tax_db_ed.csv")

# Modify the primer order
primer_order <- c("MiMammal", "MarVer1", "MarVer2", "MarVer3", "Ceto2", "Riaz12S",
                  "Mu31F_Dc320R", "Mu31F_MiMammalR", "Dc321F_Dc495R", "Dc494F_Mu643R",
                  "Mu642F_Dc1015R", "Mu642F_Mu1021R", "Dc671F_Dc1015R", "Dc671F_Mu1021R",
                  "Dc1458_1638", "Dc1465F_Dc1646R", "Mu2084F_Dc2438R", "Dc2173_2580",
                  "Mu2187F_Dc2438R", "Mu2187F_Mu2563R", "Dc2385_2580", "Dc2430_2580",
                  "Dc2438F_Mu2563R", "Dc2505_2580", "Mu9459F_Mu9822R", "Mu9459F_Mu9834R")
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
n_spp_ver <- acc_taxdb %>% filter(cat == "vertebrate") %>% pull(tax_id) %>% unique %>% length
n_seq_fis <- acc_taxdb %>% filter(cat == "fish") %>% nrow
n_spp_fis <- acc_taxdb %>% filter(cat == "fish") %>% pull(tax_id) %>% unique %>% length


# ------------------------------------------------ #
## N of sequences
# ------------------------------------------------ #
g1a <- match_seq_all %>% pivot_longer(cols = -c(n_diff, primer_name)) %>% 
  filter(name == "cetacea") %>% 
  ggplot(aes(x = n_diff, y = value)) +
  geom_bar(stat = "identity", fill = "red2") +
  geom_hline(yintercept = n_seq_cet, linetype = 2, linewidth = 0.3) +
  scale_x_discrete(limits = rev) + facet_wrap(. ~ primer_name) +
  xlab("N of total mismatches") + ylab("N of sequence amplified") +
  theme(strip.text = element_text(size = 6)) +
  panel_border() + ggtitle("N of sequences amplified, Cetacea")

g1b <- match_seq_all %>% pivot_longer(cols = -c(n_diff, primer_name)) %>% 
  filter(name == "cetacea") %>% 
  ggplot(aes(x = primer_name, y = value, fill = n_diff)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = n_seq_cet, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = rev(c4a("tol.sunset")), name = "primer\nmismatch") +
  xlab(NULL) + ylab("N of sequence amplified") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  panel_border() + ggtitle("N of sequences amplified, Cetacea")

g2a <- match_seq_all %>% pivot_longer(cols = -c(n_diff, primer_name)) %>% 
  filter(name == "vertebrate" | name == "fish") %>% 
  ggplot(aes(x = n_diff, y = value)) +
  geom_bar(stat = "identity", aes(fill = name)) +
  geom_hline(yintercept = n_seq_ver + n_seq_fis, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = c("lightblue", "burlywood2"), name = "Taxa") +
  scale_x_discrete(limits = rev) + facet_wrap(. ~ primer_name) +
  xlab("N of total mismatches") + ylab("N of sequence amplified") +
  theme(strip.text = element_text(size = 6)) +
  panel_border() + ggtitle("N of sequences amplified, Other vertebrates")

g2b <- match_seq_all %>% pivot_longer(cols = -c(n_diff, primer_name)) %>% 
  filter(name == "vertebrate") %>% 
  ggplot(aes(x = primer_name, y = value, fill = n_diff)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = n_seq_ver, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = rev(c4a("tol.sunset")), name = "primer\nmismatch") +
  xlab(NULL) + ylab("N of sequence amplified") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  panel_border() + ggtitle("N of sequences amplified, Other vertebrates")

g2c <- match_seq_all %>% pivot_longer(cols = -c(n_diff, primer_name)) %>% 
  filter(name == "fish") %>% 
  ggplot(aes(x = primer_name, y = value, fill = n_diff)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = n_seq_fis, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = rev(c4a("tol.sunset")), name = "primer\nmismatch") +
  xlab(NULL) + ylab("N of sequence amplified") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  panel_border() + ggtitle("N of sequences amplified, Fish")


# ------------------------------------------------ #
## N of species
# ------------------------------------------------ #
g3a <- match_sp_all %>% pivot_longer(cols = -c(n_diff, primer_name)) %>% 
  filter(name == "cetacea") %>% 
  ggplot(aes(x = n_diff, y = value)) +
  geom_bar(stat = "identity", fill = "red2") +
  geom_hline(yintercept = n_spp_cet, linetype = 2, linewidth = 0.3) +
  scale_x_discrete(limits = rev) + facet_wrap(. ~ primer_name) +
  xlab("N of total mismatches") + ylab("N of species amplified") +
  theme(strip.text = element_text(size = 6)) +
  panel_border() + ggtitle("N of species amplified, Cetacea")

g3b <- match_sp_all %>% pivot_longer(cols = -c(n_diff, primer_name)) %>% 
  filter(name == "cetacea") %>% 
  ggplot(aes(x = primer_name, y = value, fill = n_diff)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = n_spp_cet, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = rev(c4a("tol.sunset")), name = "primer\nmismatch") +
  xlab(NULL) + ylab("N of species amplified") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  panel_border() + ggtitle("N of species amplified, Cetacea")

g4a <- match_sp_all %>% pivot_longer(cols = -c(n_diff, primer_name)) %>% 
  filter(name == "vertebrate" | name == "fish") %>% 
  ggplot(aes(x = n_diff, y = value)) +
  geom_bar(stat = "identity", aes(fill = name)) +
  geom_hline(yintercept = n_spp_ver + n_spp_fis, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = c("lightblue", "burlywood2"), name = "Taxa") +
  scale_x_discrete(limits = rev) + facet_wrap(. ~ primer_name) +
  xlab("N of total mismatches") + ylab("N of species amplified") +
  theme(strip.text = element_text(size = 6)) +
  panel_border() + ggtitle("N of species amplified, Other vertebrates")

g4b <- match_sp_all %>% pivot_longer(cols = -c(n_diff, primer_name)) %>% 
  filter(name == "vertebrate") %>% 
  ggplot(aes(x = primer_name, y = value, fill = n_diff)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = n_spp_ver, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = rev(c4a("tol.sunset")), name = "primer\nmismatch") +
  xlab(NULL) + ylab("N of species amplified") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  panel_border() + ggtitle("N of species amplified, Other vertebrates")

g4c <- match_sp_all %>% pivot_longer(cols = -c(n_diff, primer_name)) %>% 
  filter(name == "fish") %>% 
  ggplot(aes(x = primer_name, y = value, fill = n_diff)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = n_spp_fis, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = rev(c4a("tol.sunset")), name = "primer\nmismatch") +
  xlab(NULL) + ylab("N of species amplified") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  panel_border() + ggtitle("N of species amplified, Fish")


# ------------------------------------------------ #
## Amplicon length (including primers)
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
amplicon_long$primer <- factor(amplicon_long$primer, levels = primer_order)
amplicon_others_long$primer <- factor(amplicon_others_long$primer, levels = primer_order)
## Visualize
b1 <- amplicon_long %>%
  ggplot(aes(x = primer, y = amplicon_len)) +
  geom_boxplot(outlier.size = 0, outlier.shape = NA, size = 0.2) +
  geom_jitter(alpha = 0.1, position=position_jitter(0.2)) +
  xlab(NULL) + ylab("Amplicon length (bp)") + ggtitle("Simulated amplicon length (Cetacea)") +
  scale_y_continuous(breaks = seq(0, 2000, by = 100)) + coord_cartesian(ylim = c(0,1000)) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90),
        panel.grid = element_line(linewidth = 0.05))

b2 <- amplicon_others_long %>%
  ggplot(aes(x = primer, y = amplicon_len)) +
  geom_boxplot(outlier.size = 0, outlier.shape = NA, size = 0.2) +
  geom_jitter(alpha = 0.1, position=position_jitter(0.2)) +
  xlab(NULL) + ylab("Amplicon length (bp)") + ggtitle("Simulated amplicon length (Others)") +
  scale_y_continuous(breaks = seq(0, 2000, by = 100)) + coord_cartesian(ylim = c(0,1000)) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 1, angle = 90),
        panel.grid = element_line(linewidth = 0.05))


# ------------------------------------------------ #
# Save results
# ------------------------------------------------ #
# Save figures
ggsave(file = sprintf("%s/facet_nseq_amplified.pdf", output_folder),
       plot = plot_grid(g1a, g2a, ncol = 1, align = "hv", axis = "lrbt"),
       width = 10, height = 14)
ggsave(file = sprintf("%s/facet_nspp_amplified.pdf", output_folder),
       plot = plot_grid(g3a, g4a, ncol = 1, align = "hv", axis = "lrbt"),
       width = 10, height = 14)
ggsave(file = sprintf("%s/bar_nseq_amplified.pdf", output_folder),
       plot = plot_grid(g1b, g2b, g2c, ncol = 1, align = "hv", axis = "lrbt"),
       width = 10, height = 14)
ggsave(file = sprintf("%s/bar_nspp_amplified.pdf", output_folder),
       plot = plot_grid(g3b, g4b, g4c, ncol = 1, align = "hv", axis = "lrbt"),
       width = 10, height = 14)


# Amplicon length
ggsave(file = sprintf("%s/amplicon_length.jpg", output_folder),
       plot = b1, width = 12, height = 7)
ggsave(file = sprintf("%s/amplicon_length_others.jpg", output_folder),
       plot = b2, width = 12, height = 7)
ggsave(file = sprintf("%s/amplicon_length_all_ylim1000.jpg", output_folder),
       plot = plot_grid(b1, b2, ncol = 1, align = "hv", axis = "lrbt"),
       width = 12, height = 12)

# Save sessioninfo
macam::save_session_info()

