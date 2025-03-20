####
#### Visualize other taxa information
#### 2023.07.19 Ushio
####

# Load libraries
library(tidyverse); packageVersion("tidyverse") # 1.3.2, 2023.7.18
library(cowplot); packageVersion("cowplot") # 1.1.1, 2023.7.18
library(cols4all); packageVersion("cols4all") # 0.4, 2023.7.19
#library(macam); packageVersion("macam") # 0.1.4, 2023.7.18
theme_set(theme_cowplot())

# Create output directory
set.seed(1234)
(output_folder <- macam::outdir_create())


# ------------------------------------------------ #
# Load the previous results
# ------------------------------------------------ #
df_seq <- readRDS("03_SummarizeInSilicoPCROut/match_seq_detail.obj")
df_sp <- readRDS("03_SummarizeInSilicoPCROut/match_sp_detail.obj")
acc_taxdb <- read.csv("02_CompileTaxaOut/acc_tax_db_ed.csv")
primer_order <- c("MiMammal", "MarVer1", "MarVer2", "MarVer3", "Ceto2", "Riaz12S",
                  "Mu31F_Dc320R", "Mu31F_MiMammalR", "Dc321F_Dc495R", "Dc494F_Mu643R",
                  "Mu642F_Dc1015R", "Mu642F_Mu1021R", "Dc671F_Dc1015R", "Dc671F_Mu1021R",
                  "Dc1458_1638", "Dc1465F_Dc1646R", "Mu2084F_Dc2438R", "Dc2173_2580",
                  "Mu2187F_Dc2438R", "Mu2187F_Mu2563R", "Dc2385_2580", "Dc2430_2580",
                  "Dc2438F_Mu2563R", "Dc2505_2580", "Mu9459F_Mu9822R", "Mu9459F_Mu9834R")


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
d_all$cat <- factor(d_all$cat, levels = c("cetacea", "marine_mammal", "bird", "human", "other_mammal",
                                          "reptile", "common_mammal", "fish", "amphibian"))

# Visualize
gg_obj <- d_all %>% 
  ggplot(aes(x = primer_name, y = value, fill = n_diff)) +
  geom_bar(stat = "identity") +
  facet_wrap(. ~ cat, scale = "free_y") +
  scale_fill_manual(values = rev(c4a("tol.sunset")), name = "primer\nmismatch") +
  xlab(NULL) + panel_border() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major = element_line(linewidth = 0.1, color = "gray70")) +
  ylab("N of sequence amplified") + ggtitle("N of sequences amplified")

# Forward + Reverse primers
d_all$Var1 <- factor(d_all$Var1, levels = 3:0)
d_all$Var2 <- factor(d_all$Var2, levels = 3:0)

# Forward primer only
g2 <- d_all %>% 
  ggplot(aes(x = primer_name, y = value, fill = Var1)) +
  geom_bar(stat = "identity") +
  facet_wrap(. ~ cat, scale = "free_y") +
  scale_fill_manual(values = rev(c4a("tol.sunset", 4)), name = "primer\nmismatch") +
  xlab(NULL) + panel_border() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major = element_line(linewidth = 0.1, color = "gray70")) +
  ylab("N of sequence matched by the forward primer") + ggtitle("Forward primers: N of sequences matched")

# Reverse primer only
g3 <- d_all %>% 
  ggplot(aes(x = primer_name, y = value, fill = Var2)) +
  geom_bar(stat = "identity") +
  facet_wrap(. ~ cat, scale = "free_y") +
  scale_fill_manual(values = rev(c4a("tol.sunset", 4)), name = "primer\nmismatch") +
  xlab(NULL) + panel_border() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major = element_line(linewidth = 0.1, color = "gray70")) +
  ylab("N of sequence matched by the reverse primer") + ggtitle("Reverse primers: N of sequences matched")


# ------------------------------------------------ #
# Save results
# ------------------------------------------------ #
# Save figures
ggsave(file = sprintf("%s/nseq_amplified_all.pdf", output_folder),
       plot = gg_obj,
       width = 16, height = 10)
ggsave(file = sprintf("%s/nseq_amplified_all_frPrimers.pdf", output_folder),
       plot = plot_grid(g2, g3, ncol = 1, align = "hv", axis = "lrtb"),
       width = 16, height = 20)

# Save sessioninfo
macam::save_session_info()
