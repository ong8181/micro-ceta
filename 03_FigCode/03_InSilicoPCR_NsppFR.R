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
df_sp <- readRDS("../01_InSilicoPCR/03_SummarizeInSilicoPCROut/match_sp_detail.obj")
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
# N spp
com_spp <- c("Canis lupus", "Bos taurus", "Gallus gallus", "Bos grunniens", "Equus caballus", "Capra hircus", "Ovis aries")
## Major groups
n_spp_mam <- acc_taxdb %>% filter(rep_tax == "mammal" & species != "Homo sapiens" & !(species %in% com_spp)) %>%
  pull(species) %>% unique %>% length
n_spp_cet <- acc_taxdb %>% filter(cat == "cetacea") %>% pull(species) %>% unique %>% length
n_spp_bid <- acc_taxdb %>% filter(rep_tax == "bird") %>% pull(species) %>% unique %>% length
n_spp_rep <- acc_taxdb %>% filter(rep_tax == "reptile") %>% pull(species) %>% unique %>% length
n_spp_amp <- acc_taxdb %>% filter(rep_tax == "amphibian") %>% pull(species) %>% unique %>% length
n_spp_fis <- acc_taxdb %>% filter(rep_tax == "fish") %>% pull(species) %>% unique %>% length
## Marine mammals
n_spp_mms <- acc_taxdb %>% filter(family == "Odobenidae" | family == "Otariidae" | family == "Phocidae" | genus == "Lontra" | genus == "Enhydra") %>%
  pull(species) %>% unique %>% length


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
d01b <- data.frame()
for (i in 1:length(df_sp)) {
  df_tmp <- df_sp[[i]] %>% filter(species == "Homo sapiens")
  if (nrow(df_tmp) > 0) {
    tab_seq_tmp <- table(df_tmp$n_diff_fprimer, df_tmp$n_diff_rprimer) %>% data.frame
    tab_seq_tmp$primer_name <- names(df_sp)[[i]]
    d01b <- rbind(d01b, tab_seq_tmp)
  }
}
# Extract subset data (Common mammals)
d02b <- data.frame()
for (i in 1:length(df_sp)) {
  df_tmp <- df_sp[[i]] %>% filter(species %in% com_spp)
  if (nrow(df_tmp) > 0) {
    tab_seq_tmp <- table(df_tmp$n_diff_fprimer, df_tmp$n_diff_rprimer) %>% data.frame
    tab_seq_tmp$primer_name <- names(df_sp)[[i]]
    d02b <- rbind(d02b, tab_seq_tmp)
  }
}
# Extract subset data (Other mammals)
d03b <- data.frame()
for (i in 1:length(df_sp)) {
  df_tmp <- df_sp[[i]] %>% filter(rep_tax == "mammal" & species != "Homo sapiens" & !(species %in% com_spp))
  if (nrow(df_tmp) > 0) {
    tab_seq_tmp <- table(df_tmp$n_diff_fprimer, df_tmp$n_diff_rprimer) %>% data.frame
    tab_seq_tmp$primer_name <- names(df_sp)[[i]]
    d03b <- rbind(d03b, tab_seq_tmp)
  }
}
# Extract subset data (Birds)
d04b <- data.frame()
for (i in 1:length(df_sp)) {
  df_tmp <- df_sp[[i]] %>% filter(rep_tax == "bird")
  if (nrow(df_tmp) > 0) {
    tab_seq_tmp <- table(df_tmp$n_diff_fprimer, df_tmp$n_diff_rprimer) %>% data.frame
    tab_seq_tmp$primer_name <- names(df_sp)[[i]]
    d04b <- rbind(d04b, tab_seq_tmp)
  }
}
# Extract subset data (Reptiles)
d05b <- data.frame()
for (i in 1:length(df_sp)) {
  df_tmp <- df_sp[[i]] %>% filter(rep_tax == "reptile")
  if (nrow(df_tmp) > 0) {
    tab_seq_tmp <- table(df_tmp$n_diff_fprimer, df_tmp$n_diff_rprimer) %>% data.frame
    tab_seq_tmp$primer_name <- names(df_sp)[[i]]
    d05b <- rbind(d05b, tab_seq_tmp)
  }
}
# Extract subset data (Amphibian)
d06b <- data.frame()
for (i in 1:length(df_sp)) {
  df_tmp <- df_sp[[i]] %>% filter(rep_tax == "amphibian")
  if (nrow(df_tmp) > 0) {
    tab_seq_tmp <- table(df_tmp$n_diff_fprimer, df_tmp$n_diff_rprimer) %>% data.frame
    tab_seq_tmp$primer_name <- names(df_sp)[[i]]
    d06b <- rbind(d06b, tab_seq_tmp)
  }
}
# Extract subset data (Fish)
d07b <- data.frame()
for (i in 1:length(df_sp)) {
  df_tmp <- df_sp[[i]] %>% filter(rep_tax == "fish")
  if (nrow(df_tmp) > 0) {
    tab_seq_tmp <- table(df_tmp$n_diff_fprimer, df_tmp$n_diff_rprimer) %>% data.frame
    tab_seq_tmp$primer_name <- names(df_sp)[[i]]
    d07b <- rbind(d07b, tab_seq_tmp)
  }
}
# Extract subset data (Cetacea)
d08b <- data.frame()
for (i in 1:length(df_sp)) {
  df_tmp <- df_sp[[i]] %>% filter(cat == "cetacea")
  if (nrow(df_tmp) > 0) {
    tab_seq_tmp <- table(df_tmp$n_diff_fprimer, df_tmp$n_diff_rprimer) %>% data.frame
    tab_seq_tmp$primer_name <- names(df_sp)[[i]]
    d08b <- rbind(d08b, tab_seq_tmp)
  }
}
# Extract subset data (Marine mammals)
d09b <- data.frame()
for (i in 1:length(df_sp)) {
  df_tmp <- df_sp[[i]]
  df_tmp$genus <- str_split(df_tmp$species, " ") %>% sapply(`[`, 1)
  df_tmp <- df_tmp %>% filter(family == "Odobenidae" | family == "Otariidae" | family == "Phocidae" | genus == "Lontra" | genus == "Enhydra")
  if (nrow(df_tmp) > 0) {
    tab_seq_tmp <- table(df_tmp$n_diff_fprimer, df_tmp$n_diff_rprimer) %>% data.frame
    tab_seq_tmp$primer_name <- names(df_sp)[[i]]
    d09b <- rbind(d09b, tab_seq_tmp)
  }
}

# Combine all data
d_all <- rbind(data_compiler(d01b) %>% mutate(cat = "human"),
               data_compiler(d02b) %>% mutate(cat = "common_mammal"),
               data_compiler(d03b) %>% mutate(cat = "other_mammal"),
               data_compiler(d04b) %>% mutate(cat = "bird"),
               data_compiler(d05b) %>% mutate(cat = "reptile"),
               data_compiler(d06b) %>% mutate(cat = "amphibian"),
               data_compiler(d07b) %>% mutate(cat = "fish"),
               data_compiler(d08b) %>% mutate(cat = "cetacea"),
               data_compiler(d09b) %>% mutate(cat = "marine_mammal"))


# Modify levels
tax_levels <- c("cetacea", "marine_mammal", "bird", "human", "other_mammal",
                "reptile", "common_mammal", "fish", "amphibian")
d_all$cat <- factor(d_all$cat, levels = tax_levels)
## The number of sequences for each group
hline_data <- data.frame(cat = factor(tax_levels, levels = tax_levels),
  yintercept = c(n_spp_cet, n_spp_mms, n_spp_bid, 1, n_spp_mam, n_spp_rep, length(com_spp), n_spp_fis, n_spp_amp))
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
  ylab("N of species amplified") + ggtitle("N of species amplified")

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
  ylab("N of species matched with the forward primer") + ggtitle("Forward primers: N of species matched")

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
  ylab("N of species matched with the reverse primer") + ggtitle("Reverse primers: N of species matched")




# ------------------------------------------------ #
# Save results
# ------------------------------------------------ #
# Save figures
saveRDS(list(g1, g2, g3), "data_robj/InSilicoPCR_nspp_FR.obj")

# Save sessioninfo
macam::save_session_info()
