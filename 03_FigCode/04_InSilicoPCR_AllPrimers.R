####
#### Visualize other taxa information
#### 2025.02.20, Ushio
####


# Load libraries
library(tidyverse); packageVersion("tidyverse")
library(cowplot); packageVersion("cowplot")
library(cols4all); packageVersion("cols4all")
library(ggtext); packageVersion("ggtext")
theme_set(theme_cowplot())

# Create output directory
set.seed(1234)


# ------------------------------------------------ #
# Load the previous results
# ------------------------------------------------ #
match_seq_all <- readRDS("../01_InSilicoPCR/10_SummarizeInSilicoPCR_EachOut/match_seq_all.obj")
match_sp_all <- readRDS("../01_InSilicoPCR/10_SummarizeInSilicoPCR_EachOut/match_sp_all.obj")
acc_taxdb <- read.csv("../01_InSilicoPCR/02_CompileTaxaOut/acc_tax_db_ed.csv")

# Modify the primer order
## Rename primers
### match_seq_all
match_seq_all$primer_name[match_seq_all$primer_name == "Dc2173F"] <- "Dc2173"
match_seq_all$primer_name[match_seq_all$primer_name == "Dc2385F"] <- "Dc2385"
match_seq_all$primer_name[match_seq_all$primer_name == "Dc2430F"] <- "Dc2430"
match_seq_all$primer_name[match_seq_all$primer_name == "Dc2580R"] <- "Dc2580"
match_seq_all$primer_name <- match_seq_all$primer_name %>% str_replace_all("Mu31F", "**µCeta-F | Mu31F**")
match_seq_all$primer_name <- match_seq_all$primer_name %>% str_replace_all("Dc320R", "**µCeta-R | Dc320R**")
match_seq_all$primer_name <- match_seq_all$primer_name %>% str_replace_all("Dc671F", "**Dc671F**")
match_seq_all$primer_name <- match_seq_all$primer_name %>% str_replace_all("Dc1015R", "**Dc1015R**")
match_seq_all$primer_name <- match_seq_all$primer_name %>% str_replace_all("Mu2084F", "**Mu2084F**")
match_seq_all$primer_name <- match_seq_all$primer_name %>% str_replace_all("Dc2438R", "**Dc2438R**")
match_seq_all$primer_name <- match_seq_all$primer_name %>% str_replace_all("Mu9459F", "**Mu9459F**")
match_seq_all$primer_name <- match_seq_all$primer_name %>% str_replace_all("Mu9822R", "**Mu9822R**")
### match_sp_all
match_sp_all$primer_name[match_sp_all$primer_name == "Dc2173F"] <- "Dc2173"
match_sp_all$primer_name[match_sp_all$primer_name == "Dc2385F"] <- "Dc2385"
match_sp_all$primer_name[match_sp_all$primer_name == "Dc2430F"] <- "Dc2430"
match_sp_all$primer_name[match_sp_all$primer_name == "Dc2580R"] <- "Dc2580"
match_sp_all$primer_name <- match_sp_all$primer_name %>% str_replace_all("Mu31F", "**µCeta-F | Mu31F**")
match_sp_all$primer_name <- match_sp_all$primer_name %>% str_replace_all("Dc320R", "**µCeta-R | Dc320R**")
match_sp_all$primer_name <- match_sp_all$primer_name %>% str_replace_all("Dc671F", "**Dc671F**")
match_sp_all$primer_name <- match_sp_all$primer_name %>% str_replace_all("Dc1015R", "**Dc1015R**")
match_sp_all$primer_name <- match_sp_all$primer_name %>% str_replace_all("Mu2084F", "**Mu2084F**")
match_sp_all$primer_name <- match_sp_all$primer_name %>% str_replace_all("Dc2438R", "**Dc2438R**")
match_sp_all$primer_name <- match_sp_all$primer_name %>% str_replace_all("Mu9459F", "**Mu9459F**")
match_sp_all$primer_name <- match_sp_all$primer_name %>% str_replace_all("Mu9822R", "**Mu9822R**")
## Define the primer order,
primer_order <- c("MiMammalF","MiMammalR","MarVer1F","MarVer1R","MarVer2F","MarVer2R",
                  "MarVer3F","MarVer3R","Ceto2F","Ceto2R","Riaz12SF","Riaz12SR","**µCeta-F | Mu31F**",
                  "**µCeta-R | Dc320R**","Dc321F","Dc494F","Dc495R","Mu642F","Mu643R","**Dc671F**",
                  "**Dc1015R**","Mu1021R","Dc1458F","Dc1465F","Dc1638R","Dc1646R","**Mu2084F**",
                  "Dc2173","Mu2187F","Dc2385","Dc2430","Dc2438F","**Dc2438R**","Dc2505F",
                  "Mu2563R","Dc2580","**Mu9459F**","**Mu9822R**","Mu9834R")
match_seq_all$primer_name <- factor(match_seq_all$primer_name, levels = primer_order)
match_sp_all$primer_name <- factor(match_sp_all$primer_name, levels = primer_order)

# Compile for visualization
match_seq_all$n_diff <- factor(match_seq_all$n_diff, levels = 4:0)
match_sp_all$n_diff <- factor(match_sp_all$n_diff, levels = 4:0)


# ------------------------------------------------ #
# Compile data results
# ------------------------------------------------ #
# Collect sep and species number
n_seq_cet <- acc_taxdb %>% filter(cat == "cetacea") %>% nrow
n_spp_cet <- acc_taxdb %>% filter(cat == "cetacea") %>% pull(species) %>% unique %>% length
n_seq_ver <- acc_taxdb %>% filter(cat == "vertebrate") %>% nrow
n_spp_ver <- acc_taxdb %>% filter(cat == "vertebrate") %>% pull(species) %>% unique %>% length
n_seq_fis <- acc_taxdb %>% filter(cat == "fish") %>% nrow
n_spp_fis <- acc_taxdb %>% filter(cat == "fish") %>% pull(species) %>% unique %>% length


# ------------------------------------------------ #
## N of sequences
# ------------------------------------------------ #
g1 <- match_seq_all %>% pivot_longer(cols = -c(n_diff, primer_name)) %>% 
  filter(name == "cetacea") %>% 
  ggplot(aes(x = primer_name, y = value, fill = n_diff)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = n_seq_cet, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = rev(c4a("tol.sunset", 5)), name = "primer\nmismatch") +
  xlab(NULL) + ylab("N of sequence matched") +
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5)) +
  panel_border() + ggtitle("N of sequences, Cetacea")

g2 <- match_seq_all %>% pivot_longer(cols = -c(n_diff, primer_name)) %>% 
  filter(name == "vertebrate") %>% 
  ggplot(aes(x = primer_name, y = value, fill = n_diff)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = n_seq_ver, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = rev(c4a("tol.sunset", 5)), name = "primer\nmismatch") +
  xlab(NULL) + ylab("N of sequence matched") +
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5)) +
  panel_border() + ggtitle("N of sequences, Non-fish vertebrates")

g3 <- match_seq_all %>% pivot_longer(cols = -c(n_diff, primer_name)) %>% 
  filter(name == "fish") %>% 
  ggplot(aes(x = primer_name, y = value, fill = n_diff)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = n_seq_fis, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = rev(c4a("tol.sunset", 5)), name = "primer\nmismatch") +
  xlab(NULL) + ylab("N of sequence matched") +
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5)) +
  panel_border() + ggtitle("N of sequences, Fish")


# ------------------------------------------------ #
## N of species
# ------------------------------------------------ #
g4 <- match_sp_all %>% pivot_longer(cols = -c(n_diff, primer_name)) %>% 
  filter(name == "cetacea") %>% 
  ggplot(aes(x = primer_name, y = value, fill = n_diff)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = n_spp_cet, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = rev(c4a("tol.sunset", 5)), name = "primer\nmismatch") +
  xlab(NULL) + ylab("N of species matched") +
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5)) +
  panel_border() + ggtitle("N of species, Cetacea")

g5 <- match_sp_all %>% pivot_longer(cols = -c(n_diff, primer_name)) %>% 
  filter(name == "vertebrate") %>% 
  ggplot(aes(x = primer_name, y = value, fill = n_diff)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = n_spp_ver, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = rev(c4a("tol.sunset", 5)), name = "primer\nmismatch") +
  xlab(NULL) + ylab("N of species matched") +
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5)) +
  panel_border() + ggtitle("N of species, Non-fish vertebrates")

g6 <- match_sp_all %>% pivot_longer(cols = -c(n_diff, primer_name)) %>% 
  filter(name == "fish") %>% 
  ggplot(aes(x = primer_name, y = value, fill = n_diff)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = n_spp_fis, linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = rev(c4a("tol.sunset", 5)), name = "primer\nmismatch") +
  xlab(NULL) + ylab("N of species matched") +
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5)) +
  panel_border() + ggtitle("N of species, Fish")


# ------------------------------------------------ #
# Save results
# ------------------------------------------------ #
# Save figures
saveRDS(list(g4, g6, g5), "data_robj/InSilicoPCR_nspp_AllPrimers.obj")

# Save sessioninfo
macam::save_session_info()
