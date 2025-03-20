####
#### Evaluate interspecific variations
#### 2025.02.20, Ushio
####


# Load libraries
library(tidyverse); packageVersion("tidyverse")
library(cowplot); packageVersion("cowplot")
library(ggtext); packageVersion("ggtext")
theme_set(theme_cowplot())
set.seed(1234)

# Load sequence metadata
var_all_df <- readRDS("../01_InSilicoPCR/07_EvaluateDistanceCetaOut/var_all_df.obj")


# ------------------------------------------------ #
# Visualize results
# ------------------------------------------------ #
# Define the primer order
primer_order <- c("MiMammal", "MarVer1", "MarVer2", "MarVer3", "Ceto2", "Riaz12S",
                  "**µCeta | Mu31F/Dc320R**", "Mu31F/MiMammalR", "Dc321F/Dc495R", "Dc494F/Mu643R",
                  "Mu642F/Dc1015R", "Mu642F/Mu1021R", "**Dc671F/Dc1015R**", "Dc671F/Mu1021R",
                  "Dc1458/1638", "Dc1465F/Dc1646R", "**Mu2084F/Dc2438R**", "Dc2173/2580",
                  "Mu2187F/Dc2438R", "Mu2187F/Mu2563R", "Dc2385/2580", "Dc2430/2580",
                  "Dc2438F/Mu2563R", "Dc2505/2580", "**Mu9459F/Mu9822R**", "Mu9459F/Mu9834R")

# Flatten the list
var_all_df2 <- list_rbind(var_all_df) %>% 
  pivot_longer(cols = -c(n_diff, n_diff_label, primer), names_to = "resolution", values_to = "freq")
# Adjust labels
var_all_df2$n_diff_label <- factor(var_all_df2$n_diff_label, levels = c(as.character(0:5), ">5"))
var_all_df2$resolution <- factor(var_all_df2$resolution, levels = c("species", "genus", "family"))
# Replace "_" in the primer names with "/"
var_all_df2$primer <- var_all_df2$primer %>% str_replace_all("_", "/")
var_all_df2$primer <- var_all_df2$primer %>% str_replace_all("Mu31F/Dc320R", "**µCeta | Mu31F/Dc320R**")
var_all_df2$primer <- var_all_df2$primer %>% str_replace_all("Dc671F/Dc1015R", "**Dc671F/Dc1015R**")
var_all_df2$primer <- var_all_df2$primer %>% str_replace_all("Mu2084F/Dc2438R", "**Mu2084F/Dc2438R**")
var_all_df2$primer <- var_all_df2$primer %>% str_replace_all("Mu9459F/Mu9822R", "**Mu9459F/Mu9822R**")
var_all_df2$primer <- factor(var_all_df2$primer, levels = primer_order)

# Visualize
g1 <- var_all_df2 %>%
  ggplot(aes(x = n_diff_label, y = freq + 1, fill = resolution)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("lightblue", "burlywood2", "purple1"), name = "Resolution") +
  facet_wrap(. ~ primer) + panel_border() +
  xlab("Edit distance (bases)") + ylab("N of combinations + 1") +
  theme(strip.text = element_markdown(size = 8),
        panel.grid.major.y = element_line(linewidth = 0.1, color = "gray70")) +
  scale_y_log10() + panel_border() +
  ggtitle("Edit distance between amplicons of Cetacea")

g2 <- var_all_df2 %>% filter(n_diff < 6) %>% 
  ggplot(aes(x = n_diff_label, y = freq, fill = resolution)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("lightblue", "burlywood2", "purple1"), name = "Resolution") +
  facet_wrap(. ~ primer) + panel_border() +
  xlab("Edit distance (bases)") + ylab("N of combinations") +
  theme(strip.text = element_markdown(size = 8),
        panel.grid.major.y = element_line(linewidth = 0.1, color = "gray70")) +
  panel_border() + ggtitle("Edit distance between amplicons of Cetacea")

g3 <- var_all_df2 %>% filter(n_diff < 6) %>% 
  ggplot(aes(x = n_diff_label, y = freq, fill = resolution)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("lightblue", "burlywood2", "purple1"), name = "Resolution") +
  facet_wrap(. ~ primer) + panel_border() +
  xlab("Edit distance (bases)") + ylab("N of combinations") +
  coord_cartesian(ylim = c(0,100)) +
  theme(strip.text = element_markdown(size = 8),
        panel.grid.major.y = element_line(linewidth = 0.1, color = "gray70")) +
  panel_border() + ggtitle("Edit distance between amplicons of Cetacea")


# ------------------------------------------------ #
# Save results
# ------------------------------------------------ #
# Save results
saveRDS(list(g1, g2, g3), "data_robj/InSilicoPCR_EditDist.obj")

# Save sessioninfo
macam::save_session_info()
