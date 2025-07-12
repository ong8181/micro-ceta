####
#### µCeta paper
#### Compile and format figures
#### 2025.03.07, Ushio
#### 2025.07.14, revised, Ushio
####

# Load libraries
library(tidyverse); packageVersion("tidyverse")
library(patchwork); packageVersion("patchwork")
library(cowplot); packageVersion("cowplot")
library(cols4all); packageVersion("cols4all")
library(ggforce); packageVersion("ggforce")
library(ggtext); packageVersion("ggtext")
library(png); packageVersion("png")

# Create output directory
dir.create("00_SessionInfo")
dir.create("formatted_figs")


# ------------------------------------------------ #
# Load all figure data
# ------------------------------------------------ #
# In silico PCR
p1 <- readRDS("../03_FigCode/data_robj/InSilicoPCR_nseq.obj")
p2 <- readRDS("../03_FigCode/data_robj/InSilicoPCR_nspp.obj")
p3 <- readRDS("../03_FigCode/data_robj/InSilicoPCR_nseq_FR.obj")
p4 <- readRDS("../03_FigCode/data_robj/InSilicoPCR_nspp_FR.obj")
p5 <- readRDS("../03_FigCode/data_robj/InSilicoPCR_EditDist.obj")
p6 <- readRDS("../03_FigCode/data_robj/InSilicoPCR_amplicon_length.obj")
p7 <- readRDS("../03_FigCode/data_robj/InSilicoPCR_nspp_AllPrimers.obj")
# Pool and natural test
r1 <- readRDS("../03_FigCode/data_robj/PoolTest.obj")
r2 <- readRDS("../03_FigCode/data_robj/PoolNaturalTest.obj")
# Survey image
s1 <- readPNG("../03_FigCode/data_img/SurveyPics_1.png")


# ------------------------------------------------ #
# Compile figures
# ------------------------------------------------ #
# Set figure theme
theme_set(theme_bw())

# Figure: In silico PCR
## Amplified species and sequences
Fig_pcr_nseq <-
  ((p1[[1]] + theme(axis.text.x = element_blank(),
                    plot.tag = element_text(face = "bold"),
                    legend.position = "none")) /
     (p1[[2]] + theme(axis.text.x = element_blank(),
                      plot.tag = element_text(face = "bold"),
                      legend.position = "right")) /
     (p1[[3]] + theme(plot.tag = element_text(face = "bold"),
                      legend.position = "none"))) +
  plot_annotation(tag_levels = list(c("a","b","c")))
Fig_pcr_nspp <-
  ((p2[[1]] + theme(axis.text.x = element_blank(),
                    plot.tag = element_text(face = "bold"),
                    legend.position = "none")) /
     (p2[[2]] + theme(axis.text.x = element_blank(),
                      plot.tag = element_text(face = "bold"),
                      legend.position = "right")) /
     (p2[[3]] + theme(plot.tag = element_text(face = "bold"),
                      legend.position = "none"))) +
  plot_annotation(tag_levels = list(c("a","b","c")))

## Amplified species and sequences, detailed patterns (F- and R-primers)
Fig_pcr_detailAll <- 
  ((p4[[1]] + theme(plot.tag = element_text(face = "bold"))) /
     (p3[[1]] + theme(plot.tag = element_text(face = "bold")))) +
  plot_annotation(tag_levels = list(c("a","b")))
Fig_pcr_nseqFR <-
  ((p3[[2]] + theme(plot.tag = element_text(face = "bold"))) /
     (p3[[3]] + theme(plot.tag = element_text(face = "bold")))) +
  plot_annotation(tag_levels = list(c("a","b")))
Fig_pcr_nsppFR <-
  ((p4[[2]] + theme(plot.tag = element_text(face = "bold"))) /
    (p4[[3]] + theme(plot.tag = element_text(face = "bold")))) +
  plot_annotation(tag_levels = list(c("a","b")))

## Match-mismatch (All primers)
Fig_pcr_allspp <- 
  ((p7[[1]] + theme(axis.text.x = element_blank(),
                    plot.tag = element_text(face = "bold"),
                    legend.position = "none")) /
     (p7[[2]] + theme(axis.text.x = element_blank(),
                      plot.tag = element_text(face = "bold"),
                      legend.position = "right")) /
     (p7[[3]] + theme(plot.tag = element_text(face = "bold"),
                      legend.position = "none"))) +
  plot_annotation(tag_levels = list(c("a","b","c")))

## Edit distance between amplicons
Fig_pcr_dist <- p5[[1]] + labs(title = NULL)

## Amplicon length
Fig_pcr_amplength <-
  (p6[[1]] + coord_cartesian(ylim = c(0,600)) +
     theme(plot.tag = element_text(face = "bold"),
                   axis.text.x = element_blank(),
                   panel.grid = element_line(linewidth = 0.2))) /
  (p6[[2]] + coord_cartesian(ylim = c(0,600)) +
     theme(plot.tag = element_text(face = "bold"),
                   panel.grid = element_line(linewidth = 0.2))) +
  plot_annotation(tag_levels = list(c("a","b")))

# Pool and natural water samples
## Manual colors
pool_cols <- c("#64709E","#88CCEE","#DDCC77","#339643","#CC6677","gray80", "gray40", "red3")
nat_cols <- c(c4a("carto.safe", n = 6), "gray80", "gray40", "red3")
#colorspace::choose_color()
## Pool samples
Fig_pool_cetacea1 <- r1[[2]] +
  scale_fill_manual(values = pool_cols, name = NULL) + 
  facet_wrap(~ Primer, nrow = 2) +
  labs(title = NULL) + theme(legend.position = "top")
Fig_pool_cetacea2 <- r1[[3]] +
  scale_fill_manual(values = pool_cols, name = NULL) + 
  labs(title = NULL) + theme(legend.position = "top")
## Pool and natural samples
new_levels <- c(c("Natural sample 1", "Natural sample 2", "Natural sample 3"),
                                    levels(r2[[1]]$data$new_name2)[4:10])
r2[[1]]$data$new_name2 <- factor(r2[[1]]$data$new_name2, levels = new_levels)
Fig_pool_natural1 <- r2[[1]] +
  scale_fill_manual(values = nat_cols, name = NULL) + 
  facet_wrap(~ Primer, nrow = 2) +
  labs(title = NULL) + theme(legend.position = "top")

# Enlarged figures
Fig_pool_cetacea1_zoom <- Fig_pool_cetacea1 + 
  facet_wrap(.~ Primer, nrow = 1) +
  coord_flip(ylim = c(0, 0.1)) +
  scale_y_continuous(breaks = seq(0, 0.1, by = 0.05),
                     labels = c(0, 0.05, 0.1)) +
  NULL

Fig_pool_cetacea <- (Fig_pool_cetacea1) /
  (Fig_pool_cetacea1_zoom + theme(legend.position = "none") +
     ggtitle("Enlarged figures (0-0.1)")) +
  plot_annotation(tag_levels = list(c("a","b"))) +
  plot_layout(height = c(2, 1))


# ------------------------------------------------ #
# Save compiled figures
# ------------------------------------------------ #
# Main figures
ggsave("formatted_figs/Fig2_InSilicoPCR_nspp.pdf", Fig_pcr_nspp, width = 10, height = 10)
ggsave("formatted_figs/Fig3_InSilicoPCR_dist.pdf", Fig_pcr_dist, width = 11, height = 9)
ggsave("formatted_figs/Fig5_PoolTest.pdf", Fig_pool_cetacea, width = 10, height = 10)

# Supplementary figures
ggsave("formatted_figs/FigS1_AmpliconLength.jpg", Fig_pcr_amplength, width = 8, height = 8, dpi = 300)
ggsave("formatted_figs/FigS2_InSilicoPCR_AllPrimers.pdf", Fig_pcr_allspp, width = 10, height = 10)
ggsave("formatted_figs/FigS3_InSilicoPCR_nseq.pdf", Fig_pcr_nseq, width = 10, height = 10)
ggsave("formatted_figs/FigS4_InSilicoPCR_detailAll.pdf", Fig_pcr_detailAll, width = 12, height = 15)
ggsave("formatted_figs/FigS6_PoolTest.pdf", Fig_pool_cetacea2, width = 10, height = 4)
ggsave("formatted_figs/FigS7_PoolNatTest.pdf", Fig_pool_natural1, width = 10, height = 7)


# ------------------------------------------------ #
# Save compiled data
# ------------------------------------------------ #
# Save workspace and session information
macam::save_session_info()
