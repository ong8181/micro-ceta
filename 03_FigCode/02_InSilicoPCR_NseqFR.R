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
set.seed(1234)


# ------------------------------------------------ #
# Load the previous results
# ------------------------------------------------ #
df_seq <- readRDS("../01_InSilicoPCR/03_SummarizeInSilicoPCROut/match_seq_detail.obj")
acc_taxdb <- read.csv("../01_InSilicoPCR/02_CompileTaxaOut/acc_tax_db_ed.csv")
## Define the primer order
primer_order <- c("MiMammal", "MarVer1", "MarVer2", "MarVer3", "Ceto2", "Riaz12S",
                  "**µCeta | Mu31F/Dc320R**", "Mu31F/MiMammalR", "Dc321F/Dc495R", "Dc494F/Mu643R",
                  "Mu642F/Dc1015R", "Mu642F/Mu1021R", "**Dc671F/Dc1015R**", "Dc671F/Mu1021R",
                  "Dc1458/1638", "Dc1465F/Dc1646R", "**Mu2084F/Dc2438R**", "Dc2173/2580",
                  "Mu2187F/Dc2438R", "Mu2187F/Mu2563R", "Dc2385/2580", "Dc2430/2580",
                  "Dc2438F/Mu2563R", "Dc2505/2580", "**Mu9459F/Mu9822R**", "Mu9459F/Mu9834R")



# ------------------------------------------------ #
# Compile results
# ------------------------------------------------ #
# N seq
com_spp <- c("Canis lupus", "Bos taurus", "Gallus gallus", "Bos grunniens", "Equus caballus", "Capra hircus", "Ovis aries")
## Major groups
n_seq_hum <- acc_taxdb %>% filter(species == "Homo sapiens") %>% nrow
n_seq_com <- acc_taxdb %>% filter(species %in% com_spp) %>% nrow
n_seq_mam <- acc_taxdb %>% filter(rep_tax == "mammal" & species != "Homo sapiens" & !(species %in% com_spp)) %>% nrow
n_seq_cet <- acc_taxdb %>% filter(cat == "cetacea") %>% nrow
n_seq_bid <- acc_taxdb %>% filter(rep_tax == "bird") %>% nrow
n_seq_rep <- acc_taxdb %>% filter(rep_tax == "reptile") %>% nrow
n_seq_amp <- acc_taxdb %>% filter(rep_tax == "amphibian") %>% nrow
n_seq_fis <- acc_taxdb %>% filter(rep_tax == "fish") %>% nrow
## Marine mammals
n_seq_mms <- acc_taxdb %>% filter(family == "Odobenidae" | family == "Otariidae" | family == "Phocidae" | genus == "Lontra" | genus == "Enhydra") %>% nrow


# ----------------------------------------------- #
# Compile data
# ----------------------------------------------- #
# Define wrapper function
## For data
data_compiler <- function(df) {
  # Replace "_" with "/"
  df$primer_name <- df$primer_name %>% str_replace_all("_", "/")
  df$primer_name <- df$primer_name %>% str_replace_all("Mu31F/Dc320R", "**µCeta | Mu31F/Dc320R**")
  df$primer_name <- df$primer_name %>% str_replace_all("Dc671F/Dc1015R", "**Dc671F/Dc1015R**")
  df$primer_name <- df$primer_name %>% str_replace_all("Mu2084F/Dc2438R", "**Mu2084F/Dc2438R**")
  df$primer_name <- df$primer_name %>% str_replace_all("Mu9459F/Mu9822R", "**Mu9459F/Mu9822R**")
  # Data processing
  df$primer_name <- factor(df$primer_name, levels = primer_order)
  df$Var1 <- df$Var1 %>% as.character %>% as.numeric
  df$Var2 <- df$Var2 %>% as.character %>% as.numeric
  df$n_diff <- factor(df$Var1 + df$Var2, levels = 6:0)
  return(df %>% pivot_longer(cols = -c(n_diff, primer_name, Var1, Var2)))
}


# ----------------------------------------------- #
# Data
# ----------------------------------------------- #
# Extract subset data (Human)
d01a <- data.frame()
for (i in 1:length(df_seq)) {
  df_tmp <- df_seq[[i]] %>% filter(species == "Homo sapiens")
  if (nrow(df_tmp) > 0) {
    tab_seq_tmp <- table(df_tmp$n_diff_fprimer, df_tmp$n_diff_rprimer) %>% data.frame
    tab_seq_tmp$primer_name <- names(df_seq)[[i]]
    d01a <- rbind(d01a, tab_seq_tmp)
  }
}
# Extract subset data (Common mammals)
d02a <- data.frame()
for (i in 1:length(df_seq)) {
  df_tmp <- df_seq[[i]] %>% filter(species %in% com_spp)
  if (nrow(df_tmp) > 0) {
    tab_seq_tmp <- table(df_tmp$n_diff_fprimer, df_tmp$n_diff_rprimer) %>% data.frame
    tab_seq_tmp$primer_name <- names(df_seq)[[i]]
    d02a <- rbind(d02a, tab_seq_tmp)
  }
}
# Extract subset data (Other mammals)
d03a <- data.frame()
for (i in 1:length(df_seq)) {
  df_tmp <- df_seq[[i]] %>% filter(rep_tax == "mammal" & species != "Homo sapiens" & !(species %in% com_spp))
  if (nrow(df_tmp) > 0) {
    tab_seq_tmp <- table(df_tmp$n_diff_fprimer, df_tmp$n_diff_rprimer) %>% data.frame
    tab_seq_tmp$primer_name <- names(df_seq)[[i]]
    d03a <- rbind(d03a, tab_seq_tmp)
  }
}
# Extract subset data (Birds)
d04a <- data.frame()
for (i in 1:length(df_seq)) {
  df_tmp <- df_seq[[i]] %>% filter(rep_tax == "bird")
  if (nrow(df_tmp) > 0) {
    tab_seq_tmp <- table(df_tmp$n_diff_fprimer, df_tmp$n_diff_rprimer) %>% data.frame
    tab_seq_tmp$primer_name <- names(df_seq)[[i]]
    d04a <- rbind(d04a, tab_seq_tmp)
  }
}
# Extract subset data (Reptiles)
d05a <- data.frame()
for (i in 1:length(df_seq)) {
  df_tmp <- df_seq[[i]] %>% filter(rep_tax == "reptile")
  if (nrow(df_tmp) > 0) {
    tab_seq_tmp <- table(df_tmp$n_diff_fprimer, df_tmp$n_diff_rprimer) %>% data.frame
    tab_seq_tmp$primer_name <- names(df_seq)[[i]]
    d05a <- rbind(d05a, tab_seq_tmp)
  }
}
# Extract subset data (Amphibian)
d06a <- data.frame()
for (i in 1:length(df_seq)) {
  df_tmp <- df_seq[[i]] %>% filter(rep_tax == "amphibian")
  if (nrow(df_tmp) > 0) {
    tab_seq_tmp <- table(df_tmp$n_diff_fprimer, df_tmp$n_diff_rprimer) %>% data.frame
    tab_seq_tmp$primer_name <- names(df_seq)[[i]]
    d06a <- rbind(d06a, tab_seq_tmp)
  }
}
# Extract subset data (Fish)
d07a <- data.frame()
for (i in 1:length(df_seq)) {
  df_tmp <- df_seq[[i]] %>% filter(rep_tax == "fish")
  if (nrow(df_tmp) > 0) {
    tab_seq_tmp <- table(df_tmp$n_diff_fprimer, df_tmp$n_diff_rprimer) %>% data.frame
    tab_seq_tmp$primer_name <- names(df_seq)[[i]]
    d07a <- rbind(d07a, tab_seq_tmp)
  }
}
# Extract subset data (Cetacea)
d08a <- data.frame()
for (i in 1:length(df_seq)) {
  df_tmp <- df_seq[[i]] %>% filter(cat == "cetacea")
  if (nrow(df_tmp) > 0) {
    tab_seq_tmp <- table(df_tmp$n_diff_fprimer, df_tmp$n_diff_rprimer) %>% data.frame
    tab_seq_tmp$primer_name <- names(df_seq)[[i]]
    d08a <- rbind(d08a, tab_seq_tmp)
  }
}
# Extract subset data (Marine mammals)
d09a <- data.frame()
for (i in 1:length(df_seq)) {
  df_tmp <- df_seq[[i]] %>% filter(family == "Odobenidae" | family == "Otariidae" | family == "Phocidae" | genus == "Lontra" | genus == "Enhydra")
  if (nrow(df_tmp) > 0) {
    tab_seq_tmp <- table(df_tmp$n_diff_fprimer, df_tmp$n_diff_rprimer) %>% data.frame
    tab_seq_tmp$primer_name <- names(df_seq)[[i]]
    d09a <- rbind(d09a, tab_seq_tmp)
  }
}

# Combine all data
d_all <- rbind(data_compiler(d01a) %>% mutate(cat = "human"),
               data_compiler(d02a) %>% mutate(cat = "common_mammal"),
               data_compiler(d03a) %>% mutate(cat = "other_mammal"),
               data_compiler(d04a) %>% mutate(cat = "bird"),
               data_compiler(d05a) %>% mutate(cat = "reptile"),
               data_compiler(d06a) %>% mutate(cat = "amphibian"),
               data_compiler(d07a) %>% mutate(cat = "fish"),
               data_compiler(d08a) %>% mutate(cat = "cetacea"),
               data_compiler(d09a) %>% mutate(cat = "marine_mammal"))

# Modify levels
tax_levels <- c("cetacea", "marine_mammal", "bird", "human", "other_mammal",
                "reptile", "common_mammal", "fish", "amphibian")
d_all$cat <- factor(d_all$cat, levels = tax_levels)

new_labels <- c("cetacea" = "Cetacea", "marine_mammal" = "Marine mammals", "bird" = "Bird", 
                "human" = "Human", "other_mammal" = "Other mammals", "reptile" = "Reptile", 
                "common_mammal" = "Common mammals", "fish" = "Fish", "amphibian" = "Amphibian")
## The number of sequences for each group
hline_data <- data.frame(cat = factor(tax_levels, levels = tax_levels),
                         yintercept = c(n_seq_cet, n_seq_mms, n_seq_bid, n_seq_hum, n_seq_mam, n_seq_rep, n_seq_com, n_seq_fis, n_seq_amp))
## New labels
new_labels <- c("cetacea" = "Cetacea", "marine_mammal" = "Marine mammals", "bird" = "Bird", 
                "human" = "Human", "other_mammal" = "Other mammals", "reptile" = "Reptile", 
                "common_mammal" = "Common mammals", "fish" = "Fish", "amphibian" = "Amphibian")


# Visualize
g1 <- d_all %>% 
  ggplot(aes(x = primer_name, y = value, fill = n_diff)) +
  geom_bar(stat = "identity") +
  geom_hline(data = hline_data, aes(yintercept = yintercept), linetype = 2, linewidth = 0.3) +
  facet_wrap(. ~ cat, scale = "free_y", labeller = as_labeller(new_labels)) +
  scale_fill_manual(values = rev(c4a("tol.sunset")), name = "primer\nmismatch") +
  xlab(NULL) + panel_border() +
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major = element_line(linewidth = 0.1, color = "gray70")) +
  ylab("N of sequence amplified") + ggtitle("N of sequences amplified")

# Forward + Reverse primers
d_all$Var1 <- factor(d_all$Var1, levels = 3:0)
d_all$Var2 <- factor(d_all$Var2, levels = 3:0)

# Forward primer only
g2 <- d_all %>% 
  ggplot(aes(x = primer_name, y = value, fill = Var1)) +
  geom_bar(stat = "identity") +
  geom_hline(data = hline_data, aes(yintercept = yintercept), linetype = 2, linewidth = 0.3) +
  facet_wrap(. ~ cat, scale = "free_y", labeller = as_labeller(new_labels)) +
  scale_fill_manual(values = rev(c4a("tol.sunset", 4)), name = "primer\nmismatch") +
  xlab(NULL) + panel_border() +
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major = element_line(linewidth = 0.1, color = "gray70")) +
  ylab("N of sequence matched with the forward primer") + ggtitle("Forward primers: N of sequences matched")

# Reverse primer only
g3 <- d_all %>% 
  ggplot(aes(x = primer_name, y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  geom_hline(data = hline_data, aes(yintercept = yintercept), linetype = 2, linewidth = 0.3) +
  facet_wrap(. ~ cat, scale = "free_y", labeller = as_labeller(new_labels)) +
  scale_fill_manual(values = rev(c4a("tol.sunset", 4)), name = "primer\nmismatch") +
  xlab(NULL) + panel_border() +
  theme(axis.text.x = element_markdown(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major = element_line(linewidth = 0.1, color = "gray70")) +
  ylab("N of sequence matched with the reverse primer") + ggtitle("Reverse primers: N of sequences matched")


# ------------------------------------------------ #
# Save results
# ------------------------------------------------ #
# Save figures
saveRDS(list(g1, g2, g3), "data_robj/InSilicoPCR_nseq_FR.obj")

# Save sessioninfo
macam::save_session_info()
