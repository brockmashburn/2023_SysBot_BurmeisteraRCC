
#install.packages("FactoMineR")
#install.packages("factoextra")
#install.packages("corrplot")
#install.packages("RColorBrewer")
#install.packages("ggpubr")
#install.packages("dlookr")
#install.packages("writexl")
#install.packages("ggstatsplot")
#install.packages("mclust")

################################################################################################
##method http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/118-principal-component-analysis-in-r-prcomp-vs-princomp/

library(factoextra)
library(FactoMineR)
library(mclust)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(svglite)

################################################################################
###RECURVED COROLLA CLADE DATA ANALYSIS
################################################################################

# establish working directory
setwd("/Users/brockmashburn/Documents/Work/Burmeistera/_Projects/4_RecurvedCorollaClade/Analysis")
# data for all my samples
data.all <-read.csv("Data_RCC_complete.csv")
data.all

################################################################################
## Set Theme
################################################################################

# make a theme: taken from Donoghue et al. 2022 - https://doi.org/10.1038/s41559-022-01823-x
# source code: https://github.com/eaton-lab/Oreinotinus-phylogeny/blob/main/Analyses/Morphology-ecomorphs-and-convergence/plot_igesge_mix.R
theme_clim <- function(){
  theme_bw() +
    theme(axis.text = element_text(size = 16),
          axis.text.x = element_text(face = "italic"),
          # text = element_text(family = "Arial"),
          axis.title = element_text(size = 18),
          axis.line.x = element_line(color = "black"), 
          axis.line.y = element_line(color = "black"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),                                          
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),  
          plot.margin = unit(c(1, 1, 1, 1), units = , "cm"),
          plot.title = element_text(size = 18, vjust = 1, hjust = 0),
                     legend.text = element_text(size = 12),          
                     legend.title = element_blank(),                              
                     legend.position = c(0.95, 0.15), 
                     legend.key = element_blank(),
                     legend.background = element_rect(color = "black", 
                                                      fill = "transparent", 
                                                      size = 2, linetype = "blank"),
          strip.text = element_text(size = 10, color = "black", face = "bold.italic"),
          strip.background = element_rect(color = "white", fill = "white", size = 1))
}

# load  function for geom_flat_violin
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

################################################################################
##PHASE 1: Analysis of all Recurved Corolla Clade (RCC) samples
################################################################################

###########
# Phase 1A: RCC Univariate Analysis
###########

# Characters to test:
# leaf length, leaf width, pedicel length, hyp length, hyp width, 
# calyx lobe length, calyx lobe width, corolla ventral opening, androecium length, exsertion length

### Character 1. Largest Leaf Length
names(data.all)
# select only columns of interest and remove missing samples
leaf.len <- na.omit(data.all[,c("sample_num", "pile_id1", "large_leaf_len")])

# make the plot
char.plot1 <- 
    ggplot(data = leaf.len, aes(pile_id1, large_leaf_len, fill = pile_id1)) +
    geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.7) +
    geom_point(aes(y = large_leaf_len, color = pile_id1), show.legend = FALSE,
               position = position_jitter(width = 0.15), size = 2, alpha = 0.3) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
    labs(y = "Largest Leaf Length (mm)", x = NULL) +
    guides(fill = "none", color = "none") +
    # if scale starting from 0 is desired
    #scale_y_continuous(limits = c(0, 200)) +
    scale_fill_manual(values = c("#1B9E77", "#7570B3", "#D95F02")) +
    scale_color_manual(values = c("#1B9E77", "#7570B3", "#D95F02")) +
    stat_compare_means(comparisons = list(c("B. crispiloba","B. sodiroana"),
                                          c("B. crispiloba","B. succulenta"),
                                          c("B. sodiroana","B. succulenta")),
                       #method = "wilcox.test", paired = FALSE, # these are default
                       size = 6, label = "p.signif") +
    theme_clim()
char.plot1

### Character 2. Largest Leaf Width
names(data.all)
# select only columns of interest and remove missing samples
leaf.wid <- na.omit(data.all[,c("sample_num", "pile_id1", "large_leaf_wid")])

# make the plot
char.plot2 <- 
    ggplot(data = leaf.wid, aes(pile_id1, large_leaf_wid, fill = pile_id1)) +
    geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.7) +
    geom_point(aes(y = large_leaf_wid, color = pile_id1), show.legend = FALSE,
               position = position_jitter(width = 0.15), size = 2, alpha = 0.3) +
    geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
    labs(y = "Largest Leaf Width (mm)", x = NULL) +
    guides(fill = "none", color = "none") +
    # if scale starting from 0 is desired
    #scale_y_continuous(limits = c(0, 200)) +
    scale_fill_manual(values = c("#1B9E77", "#7570B3", "#D95F02")) +
    scale_color_manual(values = c("#1B9E77", "#7570B3", "#D95F02")) +
    stat_compare_means(comparisons = list(c("B. crispiloba","B. sodiroana"),
                                          c("B. crispiloba","B. succulenta"),
                                          c("B. sodiroana","B. succulenta")),
                       #method = "wilcox.test", paired = FALSE, # these are default
                       size = 6, label = "p.signif") +
  theme_clim()
char.plot2

### Character 3. Pedicel Length
names(data.all)
# select only columns of interest and remove missing samples
pedicel.len <- na.omit(data.all[,c("sample_num", "pile_id1", "pedicel_len")])

# make the plot
char.plot3 <- 
  ggplot(data = pedicel.len, aes(pile_id1, pedicel_len, fill = pile_id1)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.7) +
  geom_point(aes(y = pedicel_len, color = pile_id1), show.legend = FALSE,
             position = position_jitter(width = 0.15), size = 2, alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  labs(y = "Pedicel Length (mm)", x = NULL) +
  guides(fill = "none", color = "none") +
  # if scale starting from 0 is desired
  #scale_y_continuous(limits = c(0, 200)) +
  scale_fill_manual(values = c("#1B9E77", "#7570B3", "#D95F02")) +
  scale_color_manual(values = c("#1B9E77", "#7570B3", "#D95F02")) +
  stat_compare_means(comparisons = list(c("B. crispiloba","B. sodiroana"),
                                        c("B. crispiloba","B. succulenta"),
                                        c("B. sodiroana","B. succulenta")),
                     #method = "wilcox.test", paired = FALSE, # these are default
                     size = 6, label = "p.signif") +
  theme_clim()
char.plot3

### Character 4. Hypanthium Length
names(data.all)
# select only columns of interest and remove missing samples
hyp.len <- na.omit(data.all[,c("sample_num", "pile_id1", "hyp_len")])

# make a plot
char.plot4 <- 
  ggplot(data = hyp.len, aes(pile_id1, hyp_len, fill = pile_id1)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.7) +
  geom_point(aes(y = hyp_len, color = pile_id1), show.legend = FALSE,
             position = position_jitter(width = 0.15), size = 2, alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  labs(y = "Hypanthium Length (mm)", x = NULL) +
  guides(fill = "none", color = "none") +
  # if scale starting from 0 is desired
  #scale_y_continuous(limits = c(0, 200)) +
  scale_fill_manual(values = c("#1B9E77", "#7570B3", "#D95F02")) +
  scale_color_manual(values = c("#1B9E77", "#7570B3", "#D95F02")) +
  stat_compare_means(comparisons = list(c("B. crispiloba","B. sodiroana"),
                                        c("B. crispiloba","B. succulenta"),
                                        c("B. sodiroana","B. succulenta")),
                     #method = "wilcox.test", paired = FALSE, # these are default
                     size = 6, label = "p.signif") +
  theme_clim()
char.plot4

### Character 5. Hypanthium Width
names(data.all)
# select only columns of interest and remove missing samples
hyp.wid <- na.omit(data.all[,c("sample_num", "pile_id1", "hyp_wid")])

# make the plot
char.plot5 <- 
  ggplot(data = hyp.wid, aes(pile_id1, hyp_wid, fill = pile_id1)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.7) +
  geom_point(aes(y = hyp_wid, color = pile_id1), show.legend = FALSE,
             position = position_jitter(width = 0.15), size = 2, alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  labs(y = "Hypanthium Width (mm)", x = NULL) +
  guides(fill = "none", color = "none") +
  # if scale starting from 0 is desired
  #scale_y_continuous(limits = c(0, 200)) +
  scale_fill_manual(values = c("#1B9E77", "#7570B3", "#D95F02")) +
  scale_color_manual(values = c("#1B9E77", "#7570B3", "#D95F02")) +
  stat_compare_means(comparisons = list(c("B. crispiloba","B. sodiroana"),
                                        c("B. crispiloba","B. succulenta"),
                                        c("B. sodiroana","B. succulenta")),
                     #method = "wilcox.test", paired = FALSE, # these are default
                     size = 6, label = "p.signif") +
  theme_clim()
char.plot5

### Character 6. Calyx Lobe Length
names(data.all)
# select only columns of interest and remove missing samples
calyxlobe.len <- na.omit(data.all[,c("sample_num", "pile_id1", "calyxlobe_len")])

# make the plot
char.plot6 <- 
  ggplot(data = calyxlobe.len, aes(pile_id1, calyxlobe_len, fill = pile_id1)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.7) +
  geom_point(aes(y = calyxlobe_len, color = pile_id1), show.legend = FALSE,
             position = position_jitter(width = 0.15), size = 2, alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  labs(y = "Calyx Lobe Length (mm)", x = NULL) +
  guides(fill = "none", color = "none") +
  # if scale starting from 0 is desired
  #scale_y_continuous(limits = c(0, 200)) +
  scale_fill_manual(values = c("#1B9E77", "#7570B3", "#D95F02")) +
  scale_color_manual(values = c("#1B9E77", "#7570B3", "#D95F02")) +
  stat_compare_means(comparisons = list(c("B. crispiloba","B. sodiroana"),
                                        c("B. crispiloba","B. succulenta"),
                                        c("B. sodiroana","B. succulenta")),
                     #method = "wilcox.test", paired = FALSE, # these are default
                     size = 6, label = "p.signif") +
  theme_clim()
char.plot6

### Character 7. Calyx Lobe Width
names(data.all)
# select only columns of interest and remove missing samples
calyxlobe.wid <- na.omit(data.all[,c("sample_num", "pile_id1", "calyxlobe_wid")])

# make the plot
char.plot7 <- 
  ggplot(data = calyxlobe.wid, aes(pile_id1, calyxlobe_wid, fill = pile_id1)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.7) +
  geom_point(aes(y = calyxlobe_wid, color = pile_id1), show.legend = FALSE,
             position = position_jitter(width = 0.15), size = 2, alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  labs(y = "Calyx Lobe Width (mm)", x = NULL) +
  guides(fill = "none", color = "none") +
  # if scale starting from 0 is desired
  #scale_y_continuous(limits = c(0, 200)) +
  scale_fill_manual(values = c("#1B9E77", "#7570B3", "#D95F02")) +
  scale_color_manual(values = c("#1B9E77", "#7570B3", "#D95F02")) +
  stat_compare_means(comparisons = list(c("B. crispiloba","B. sodiroana"),
                                        c("B. crispiloba","B. succulenta"),
                                        c("B. sodiroana","B. succulenta")),
                     #method = "wilcox.test", paired = FALSE, # these are default
                     size = 6, label = "p.signif") +
  theme_clim()
char.plot7

### Character 8. Corolla Ventral Opening
names(data.all)
# select only columns of interest and remove missing samples
cortube.vent.len <- na.omit(data.all[,c("sample_num", "pile_id1", "cortube_vent_len")])

# make the plot
char.plot8 <- 
  ggplot(data = cortube.vent.len, aes(pile_id1, cortube_vent_len, fill = pile_id1)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.7) +
  geom_point(aes(y = cortube_vent_len, color = pile_id1), show.legend = FALSE,
             position = position_jitter(width = 0.15), size = 2, alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  labs(y = "Length to Ventral Corolla Opening (mm)", x = NULL) +
  guides(fill = "none", color = "none") +
  # if scale starting from 0 is desired
  #scale_y_continuous(limits = c(0, 200)) +
  scale_fill_manual(values = c("#1B9E77", "#7570B3", "#D95F02")) +
  scale_color_manual(values = c("#1B9E77", "#7570B3", "#D95F02")) +
  stat_compare_means(comparisons = list(c("B. crispiloba","B. sodiroana"),
                                        c("B. crispiloba","B. succulenta"),
                                        c("B. sodiroana","B. succulenta")),
                     #method = "wilcox.test", paired = FALSE, # these are default
                     size = 6, label = "p.signif") +
  theme_clim()
char.plot8

### Character 9. Androecium Length
names(data.all)
# select only columns of interest and remove missing samples
andr.len <- na.omit(data.all[,c("sample_num", "pile_id1", "andr_len")])

# make the plot
char.plot9 <- 
  ggplot(data = andr.len, aes(pile_id1, andr_len, fill = pile_id1)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.7) +
  geom_point(aes(y = andr_len, color = pile_id1), show.legend = FALSE,
             position = position_jitter(width = 0.15), size = 2, alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  labs(y = "Androecium Length (mm)", x = NULL) +
  guides(fill = "none", color = "none") +
  # if scale starting from 0 is desired
  #scale_y_continuous(limits = c(0, 200)) +
  scale_fill_manual(values = c("#1B9E77", "#7570B3", "#D95F02")) +
  scale_color_manual(values = c("#1B9E77", "#7570B3", "#D95F02")) +
  stat_compare_means(comparisons = list(c("B. crispiloba","B. sodiroana"),
                                        c("B. crispiloba","B. succulenta"),
                                        c("B. sodiroana","B. succulenta")),
                     #method = "wilcox.test", paired = FALSE, # these are default
                     size = 6, label = "p.signif") +
  theme_clim()
char.plot9

### Character 10. Exsertion Length
names(data.all)
# select only columns of interest and remove missing samples
exsertion.len <- na.omit(data.all[,c("sample_num", "pile_id1", "exsertion_len")])

# make the plot
char.plot10 <- 
  ggplot(data = exsertion.len, aes(pile_id1, exsertion_len, fill = pile_id1)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.7) +
  geom_point(aes(y = exsertion_len, color = pile_id1), show.legend = FALSE,
             position = position_jitter(width = 0.15), size = 2, alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  labs(y = "Exsertion Length (mm)", x = NULL) +
  guides(fill = "none", color = "none") +
  # if scale starting from 0 is desired
  #scale_y_continuous(limits = c(0, 200)) +
  scale_fill_manual(values = c("#1B9E77", "#7570B3", "#D95F02")) +
  scale_color_manual(values = c("#1B9E77", "#7570B3", "#D95F02")) +
  stat_compare_means(comparisons = list(c("B. crispiloba","B. sodiroana"),
                                        c("B. crispiloba","B. succulenta"),
                                        c("B. sodiroana","B. succulenta")),
                     #method = "wilcox.test", paired = FALSE, # these are default
                     size = 6, label = "p.signif") +
  theme_clim()
char.plot10

### Figure 2
## export individual plots to make Figure 2 in inkscape
#calyx lobe length
ggsave(char.plot6, device = "svg", filename = "F2-CalyxLobeLength.svg",
       height = 6, width = 6, dpi = 600)
#androecium length
ggsave(char.plot9, device = "svg", filename = "F2-AndroeciumLength.svg",
       height = 6, width = 6, dpi = 600)
#exsertion length
ggsave(char.plot10, device = "svg", filename = "F2-ExsertionLength.svg",
       height = 6, width = 6, dpi = 600)

### Figure S1
## make a figure of the remaining characters compared
# Make a panel
multi_plot1 <- ggarrange(char.plot1,char.plot2,char.plot3,
                         char.plot4,char.plot5,char.plot7,
                         char.plot8,
                         labels = c("A","B","C","D","E","F","G"), #labels given each panel 
                         font.label = list(size = 30, color = "black", family = NULL),
                         ncol = 3, nrow = 3, #adjust plot space 
                         common.legend = FALSE) #does the plot have a common legend
#add titles and labels to the multi-panel graph
multi_plot1 <- annotate_figure(multi_plot1,
                               top = text_grob("Figure S1. Phase 1 Univariate Analysis: Additional Traits Measured", 
                                               color = "black", face = "bold", size = 34))
# save Figure S1
ggsave(multi_plot1, device = "tiff", filename = "FigureS1.tiff", 
       height = 20, width = 18, dpi = 350, bg = "white")
ggsave(multi_plot1, device = "pdf", filename = "FigureS1.pdf",
       height = 20, width = 18, dpi = 350, bg = "white")

###########
# Phase 1B: RCC Multivariate Analysis
###########

# make a data frame to hold clustering outputs
RCC.cluster.res <- data.all[,c(1:6)]

# measurement headers
names(data.all)

# dataset of measurements
data.RCC <- na.omit(data.all[,c(1,2,4:6,9,14,15,20:25,27,29,32)])
data.RCC

names(data.RCC)

#run PCA
RCC.pca <- PCA(data.RCC[ , -c(1:6,16)], ncp = 6, graph = FALSE)
RCC.pca
RCC.pca$eig
RCC.pca$var$cos2

RCC.PCA.1 <- fviz_pca_ind(RCC.pca, 
                        geom.ind = c("point","text"), pointsize = 2.5,
                        mean.point = FALSE, repel = TRUE,
                        col.ind = data.RCC$pile_id)
RCC.PCA.1 <- ggpubr::ggpar(RCC.PCA.1,
                   title = "Phase 1 PCA: Sample Numbers",
                   font.main = c(18, "bold", "black"),
                   font.x = 12, font.y = 12,
                   legend = c(0.125, 0.12),
                   font.legend = c(12, "italic"),
                   ggtheme = theme_bw(),
                   palette = c("#1B9E77", "#7570B3", "#D95F02")) +
  theme(legend.title = element_blank())
RCC.PCA.1

#visualize the variables
RCC.PCA.2 <- fviz_pca_biplot(RCC.pca, fill.ind = "white", label = "var", col.var = "black",
                     col.ind = data.RCC$pile_id, alpha.ind = 1, 
                     mean.point = FALSE,
                     pointsize = 3, repel = TRUE)
RCC.PCA.2 <- ggpubr::ggpar(RCC.PCA.2, title = "PCA of the RCC - contributions",
                   font.x = 12, font.y = 12,
                   legend = "none",
                   ggtheme = theme_bw(), 
                   palette = c("#1B9E77", "#7570B3", "#D95F02")) +
  theme(legend.title = element_blank())
RCC.PCA.2

#heirarchical clustering unsupervised
hc.RCC.pca <- HCPC(RCC.pca, nb.clust = -1)
RCC.PCA.3 <- fviz_cluster(hc.RCC.pca, repel = TRUE, pointsize = 3, 
                  show.clust.cent = FALSE, geom = "point",
                  ggtheme = theme_bw(), palette = "jco")
RCC.PCA.3 <- ggpubr::ggpar(RCC.PCA.3, title = "PCA of RCC - clustering",
                   font.x = 12, font.y = 12,
                   legend = "none",
                   )
RCC.PCA.3
hc.RCC.pca$desc.var$quanti
hc.RCC.pca$data.clust$clust

# heirarchical clustering unsupervised - with sample numbers
RCC.PCA.4 <- fviz_cluster(hc.RCC.pca, repel = TRUE, pointsize = 2.5, 
                          show.clust.cent = FALSE, geom = c("point","text"),
                          ggtheme = theme_bw(), palette = "jco",
                          max.overlaps = 20)
RCC.PCA.4 <- ggpubr::ggpar(RCC.PCA.4, 
                           title = "Phase 1 PCA: Clusters with Sample Numbers",
                           font.main = c(18, "bold", "black"),
                           legend = c(0.07, 0.14),
                           font.legend = c(12),
                           font.x = 12, font.y = 12,
                           ggtheme = theme_bw(),
                           )
RCC.PCA.4

# save RCC PCA clustering results to data frame
RCC.cluster.res$RCC_pca_clust <- hc.RCC.pca$data.clust$clust

###Factor Analysis of Mixed Data on RCC
names(data.RCC)
RCC.famd <- FAMD(data.RCC[ , -c(1:5)], ncp = 6, graph = FALSE)
RCC.famd$eig

# FAMD plot with labels
RCC.FAMD.1 <- fviz_famd_ind(RCC.famd, habillage = as.factor(data.RCC$pile_id1),
                          mean.point = FALSE,
                          invisible = "quali.var",
                          geom = c("point", "text"), pointsize = 2.5, repel = TRUE)
RCC.FAMD.1 <- ggpubr::ggpar(RCC.FAMD.1,
                   title = "Phase 1 FAMD: Sample Numbers",
                   font.main = c(18, "bold", "black"),
                   font.x = 12, font.y = 12,
                   legend = c(0.87, 0.1),
                   font.legend = c(12, "italic"),
                   ggtheme = theme_bw(), 
                   palette = c("#1B9E77", "#7570B3", "#D95F02")) + 
  theme(legend.title = element_blank())
RCC.FAMD.1

# look at the contribution of variables
RCC.FAMD.2 <- fviz_famd_ind(RCC.famd,
                   habillage = as.factor(data.RCC$pile_id1),
                   col.quali.var = "black", 
                   #goem = c("point", "text"),
                   label = "quali.var",
                   mean.point = FALSE,
                   alpha.ind = 1, repel = TRUE)
RCC.FAMD.2 <- ggpubr::ggpar(RCC.FAMD.2, title = "FAMD of RCC - contributions",
                   font.x = 12, font.y = 12,
                   legend = "none", orientation = "reverse",
                   ggtheme = theme_bw(), 
                   palette = c("#1B9E77", "#7570B3", "#D95F02")) + 
  theme(legend.title = element_blank())
RCC.FAMD.2

# heirarchical clustering unsupervised
hc.RCC.famd <- HCPC(RCC.famd, nb.clust = -1)
RCC.FAMD.3 <- fviz_cluster(hc.RCC.famd, repel = TRUE, geom = "point",
                  show.clust.cent = FALSE, pointsize = 3,
                  ggtheme = theme_bw(), palette = "jco")
RCC.FAMD.3 <- ggpubr::ggpar(RCC.FAMD.3, title = "FAMD of RCC - Clusters",
                   font.x = 12, font.y = 12,
                   orientation = "reverse",
                   legend = "none")
RCC.FAMD.3
hc.RCC.famd$desc.var$quanti
hc.RCC.famd$data.clust$clust

# heirarchical clustering unsupervised with sample labels
RCC.FAMD.4 <- fviz_cluster(hc.RCC.famd, repel = TRUE, geom = c("point","text"),
                           show.clust.cent = FALSE, pointsize = 2.5,
                           ggtheme = theme_bw(), palette = "jco")
RCC.FAMD.4 <- ggpubr::ggpar(RCC.FAMD.4, title = "Phase 1 FAMD: Clusters with Sample Numbers",
                            font.title = c(18, "bold", "black"), 
                            font.x = 12, font.y = 12,
                            legend = c(0.92, 0.86),
                            font.legend = c(12),
                            orientation = "reverse")
RCC.FAMD.4

# save RCC FAMD clustering results to data frame
RCC.cluster.res$RCC_famd_clust <- hc.RCC.famd$data.clust$clust

##### Model Based Clustering
names(data.RCC)
data.RCC.norm <- scale(data.RCC[, -c(1:6,16)])
mc.RCC <- Mclust(data.RCC.norm)
summary(mc.RCC)
mc.RCC$classification
# what about 3 clusters?
mc.RCC.3 <- Mclust(data.RCC.norm, G = 3)
fviz_mclust(mc.RCC.3, "classification", 
            ggtheme = theme_bw(), palette = "jco")
summary(mc.RCC.3)

# save RCC model based clustering results to data frame
RCC.cluster.res$RCC_mbc_k2 <- mc.RCC$classification
RCC.cluster.res$RCC_mbc_k3 <- mc.RCC.3$classification

# plot the BIC scores
RCC.MBC.1 <- fviz_mclust(mc.RCC, "BIC", 
                 ggtheme = theme_bw(), palette = "jco")
RCC.MBC.1 <- ggpubr::ggpar(RCC.MBC.1, title = "Phase 1 Model-based Clust.: Model Selection",
                   font.main = c(18, "bold", "black"),
                   font.x = 12, font.y = 12)
RCC.MBC.1

# plot 2 clusters in clean format
RCC.MBC.2 <- fviz_mclust(mc.RCC, "classification",
                 show.clust.cent = FALSE,
                 geom = c("point"), pointsize = 3,
                 repel= TRUE, ellipse.type = "convex",
                 ggtheme = theme_bw(), palette = "jco")
RCC.MBC.2 <- ggpubr::ggpar(RCC.MBC.2, title = "MClust of RCC - K = 2", 
                   font.x = 12, font.y = 12,
                   legend = "none",
                   ggtheme = theme_bw()) + 
  theme(plot.subtitle = element_blank()) +
  scale_x_reverse()
RCC.MBC.2

# plot 2 clusters with sample numbers included
RCC.MBC.3 <- fviz_mclust(mc.RCC, "classification",
                         show.clust.cent = FALSE,
                         geom = c("point", "text"), pointsize = 2.5,
                         repel= TRUE, ellipse.type = "convex",
                         ggtheme = theme_bw(), palette = "jco")
RCC.MBC.3 <- ggpubr::ggpar(RCC.MBC.3, title = "Phase 1 Model-based Clust.: n = 2, Sample Numbers",
                           font.main = c(18, "bold", "black"),
                           font.x = 12, font.y = 12,
                           legend = c(0.08, 0.1),
                           font.legend = 12,
                           ggtheme = theme_bw()) + 
  theme(plot.subtitle = element_blank()) +
  scale_x_reverse()
RCC.MBC.3

# plot 3 clusters in clean format
RCC.MBC.4 <- fviz_mclust(mc.RCC.3, "classification",
                 show.clust.cent = FALSE,
                 geom = c("point"), pointsize = 3,
                 repel= TRUE, ellipse.type = "convex",
                 ggtheme = theme_bw(), palette = "jco")
RCC.MBC.4 <- ggpubr::ggpar(RCC.MBC.4, title = "Mclust of RCC - n = 3",
                   show.clust.cent = FALSE,
                   font.x = 12, font.y = 12,
                   legend = "none",
                   ggtheme = theme_bw()) +
  theme(plot.subtitle = element_blank()) +
  scale_x_reverse()
RCC.MBC.4

### save the clustering results as a csv
#write.csv(RCC.cluster.res, "RCC_cluster_results.csv")

### Figure 3
#RCC - PCA with variables plotted
ggsave(RCC.PCA.2, device = "svg", filename = "F3-A.svg",
       height = 6, width = 6, dpi = 600)
#RCC - PCA with clusters
ggsave(RCC.PCA.3, device = "svg", filename = "F3-B.svg",
       height = 6, width = 6, dpi = 600)
##RCC - FAMD with variables plotted
ggsave(RCC.FAMD.2, device = "svg", filename = "F3-C.svg",
       height = 6, width = 6, dpi = 600)
##RCC - FAMD with clusters
ggsave(RCC.FAMD.3, device = "svg", filename = "F3-D.svg",
       height = 6, width = 6, dpi = 600)
##RCC - MClust with 2 clusters
ggsave(RCC.MBC.2, device = "svg", filename = "F3-E.svg",
       height = 6, width = 6, dpi = 600)
##RCC - MClust with 3 clusters
ggsave(RCC.MBC.4, device = "svg", filename = "F3-F.svg",
       height = 6, width = 6, dpi = 600)

### Figure S2
# put clustering results together
multi_plot2 <- ggarrange(RCC.PCA.1,RCC.PCA.4,
                         RCC.FAMD.1,RCC.FAMD.4,
                         RCC.MBC.1,RCC.MBC.3,
                         labels = c("A","B","C","D","E","F"), #labels given each panel
                         font.label = list(size = 30, color = "black", family = NULL),
                         ncol = 2, nrow = 3, #adjust plot space 
                         common.legend = FALSE) #does the plot have a common legend
#add titles and labels to the multi-panel graph
multi_plot2 <- annotate_figure(multi_plot2,
                               top = text_grob("Figure S2. Phase 1 Multivariate Analysis", 
                                               color = "black", face = "bold", size = 26))
# save Figure S2
ggsave(multi_plot2, device = "tiff", filename = "FigureS2.tiff", 
       height = 19, width = 14, dpi = 350, bg = "white")
ggsave(multi_plot2, device = "pdf", filename = "FigureS2.pdf", 
       units = "mm", scale = 2, height = 240, width = 177, dpi = 350, bg = "white")

################################################################################
##PHASE 2: Analysis of B. succulenta s.l. (SUC)
################################################################################
# list samples by species name
data.RCC$pile_id1
# keep just the succulenta samples
data.SUC <- data.RCC[c(69:95),]
names(data.SUC)

# make a data frame to hold clustering outputs
SUC.cluster.res <- data.SUC[,c(1:5)]

# change the names in the country column
for (i in 1:nrow(data.SUC)) {
  if(data.SUC$country[i] == "Venezuela") {
    data.SUC$country[i] <- "Col. & Ven."
  } else if(data.SUC$country[i] == "Colombia") {
    data.SUC$country[i] <- "Col. & Ven."
  }
}
# confirm
data.SUC$country

###########
# Phase 2A: SUC Univariate Analysis
###########

### Character 1. Largest Leaf Length - not significant
names(data.SUC)
# select only columns of interest and remove missing samples
suc.leaf.len <- na.omit(data.SUC[,c("sample_num", "country", "large_leaf_len")])

# make the plot
suc.char.plot1 <- 
  ggplot(data = suc.leaf.len, aes(country, large_leaf_len, fill = country)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.7) +
  geom_point(aes(y = large_leaf_len, color = country), show.legend = FALSE,
             position = position_jitter(width = 0.15), size = 2, alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  labs(y = "Largest Leaf Length (mm)", x = NULL) +
  guides(fill = "none", color = "none") +
  # if scale starting from 0 is desired
  #scale_y_continuous(limits = c(0, 200)) +
  scale_fill_manual(values = c("#1B9E77", "#0072B2")) +
  scale_color_manual(values = c("#1B9E77", "#0072B2")) +
  stat_compare_means(comparisons = list(c("Ecuador","Col. & Ven.")),
                     #method = "wilcox.test", paired = FALSE, # these are default
                     size = 6, label = "p.signif") +
  theme_clim()
suc.char.plot1

### Character 2. Largest Leaf Width - not significant
names(data.SUC)
# select only columns of interest and remove missing samples
suc.leaf.wid <- na.omit(data.SUC[,c("sample_num", "country", "large_leaf_wid")])

# make the plot
suc.char.plot2 <- 
  ggplot(data = suc.leaf.wid, aes(country, large_leaf_wid, fill = country)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.7) +
  geom_point(aes(y = large_leaf_wid, color = country), show.legend = FALSE,
             position = position_jitter(width = 0.15), size = 2, alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  labs(y = "Largest Leaf Width (mm)", x = NULL) +
  guides(fill = "none", color = "none") +
  # if scale starting from 0 is desired
  #scale_y_continuous(limits = c(0, 200)) +
  scale_fill_manual(values = c("#1B9E77", "#0072B2")) +
  scale_color_manual(values = c("#1B9E77", "#0072B2")) +
  stat_compare_means(comparisons = list(c("Ecuador","Col. & Ven.")),
                     #method = "wilcox.test", paired = FALSE, # these are default
                     size = 6, label = "p.signif") +
  theme_clim()
suc.char.plot2

### Character 3. Pedicel Length - significant
names(data.SUC)
# select only columns of interest and remove missing samples
pedicel.len <- na.omit(data.SUC[,c("sample_num", "country", "pedicel_len")])

# make the plot
suc.char.plot3 <- 
  ggplot(data = pedicel.len, aes(country, pedicel_len, fill = country)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.7) +
  geom_point(aes(y = pedicel_len, color = country), show.legend = FALSE,
             position = position_jitter(width = 0.15), size = 2, alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  labs(y = "Pedicel Length (mm)", x = NULL) +
  guides(fill = "none", color = "none") +
  # if scale starting from 0 is desired
  #scale_y_continuous(limits = c(0, 200)) +
  scale_fill_manual(values = c("#1B9E77", "#0072B2")) +
  scale_color_manual(values = c("#1B9E77", "#0072B2")) +
  stat_compare_means(comparisons = list(c("Ecuador","Col. & Ven.")),
                     #method = "wilcox.test", paired = FALSE, # these are default
                     size = 6, label = "p.signif") +
  theme_clim()
suc.char.plot3

### Character 4. Hypanthium Length - significant
names(data.SUC)
# select only columns of interest and remove missing samples
suc.hyp.len <- na.omit(data.SUC[,c("sample_num", "country", "hyp_len")])

# make the plot
suc.char.plot4 <- 
  ggplot(data = suc.hyp.len, aes(country, hyp_len, fill = country)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.7) +
  geom_point(aes(y = hyp_len, color = country), show.legend = FALSE,
             position = position_jitter(width = 0.15), size = 2, alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  labs(y = "Hypanthium Length (mm)", x = NULL) +
  guides(fill = "none", color = "none") +
  # if scale starting from 0 is desired
  #scale_y_continuous(limits = c(0, 200)) +
  scale_fill_manual(values = c("#1B9E77", "#0072B2")) +
  scale_color_manual(values = c("#1B9E77", "#0072B2")) +
  stat_compare_means(comparisons = list(c("Ecuador","Col. & Ven.")),
                     #method = "wilcox.test", paired = FALSE, # these are default
                     size = 6, label = "p.signif") +
  theme_clim()
suc.char.plot4

### Character 5. Hypanthium Width - not significant
names(data.SUC)
# select only columns of interest and remove missing samples
suc.hyp.wid <- na.omit(data.SUC[,c("sample_num", "country", "hyp_wid")])

# make the plot
suc.char.plot5 <- 
  ggplot(data = suc.hyp.wid, aes(country, hyp_wid, fill = country)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.7) +
  geom_point(aes(y = hyp_wid, color = country), show.legend = FALSE,
             position = position_jitter(width = 0.15), size = 2, alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  labs(y = "Hypanthium Width (mm)", x = NULL) +
  guides(fill = "none", color = "none") +
  # if scale starting from 0 is desired
  #scale_y_continuous(limits = c(0, 200)) +
  scale_fill_manual(values = c("#1B9E77", "#0072B2")) +
  scale_color_manual(values = c("#1B9E77", "#0072B2")) +
  stat_compare_means(comparisons = list(c("Ecuador","Col. & Ven.")),
                     #method = "wilcox.test", paired = FALSE, # these are default
                     size = 6, label = "p.signif") +
  theme_clim()
suc.char.plot5

### Character 6. Calyx Lobe Length - not significant
names(data.SUC)
# select only columns of interest and remove missing samples
suc.calyxlobe.len <- na.omit(data.SUC[,c("sample_num", "country", "calyxlobe_len")])

# make the plot
suc.char.plot6 <- 
  ggplot(data = suc.calyxlobe.len, aes(country, calyxlobe_len, fill = country)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.7) +
  geom_point(aes(y = calyxlobe_len, color = country), show.legend = FALSE,
             position = position_jitter(width = 0.15), size = 2, alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  labs(y = "Calyx Lobe Length (mm)", x = NULL) +
  guides(fill = "none", color = "none") +
  # if scale starting from 0 is desired
  #scale_y_continuous(limits = c(0, 200)) +
  scale_fill_manual(values = c("#1B9E77", "#0072B2")) +
  scale_color_manual(values = c("#1B9E77", "#0072B2")) +
  stat_compare_means(comparisons = list(c("Ecuador","Col. & Ven.")),
                     #method = "wilcox.test", paired = FALSE, # these are default
                     size = 6, label = "p.signif") +
  theme_clim()
suc.char.plot6

### Character 7. Calyx Lobe Width - significant
names(data.SUC)
# select only columns of interest and remove missing samples
suc.calyxlobe.wid <- na.omit(data.SUC[,c("sample_num", "country", "calyxlobe_wid")])

# make the plot
suc.char.plot7 <- 
  ggplot(data = suc.calyxlobe.wid, aes(country, calyxlobe_wid, fill = country)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.7) +
  geom_point(aes(y = calyxlobe_wid, color = country), show.legend = FALSE,
             position = position_jitter(width = 0.15), size = 2, alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  labs(y = "Calyx Lobe Width (mm)", x = NULL) +
  guides(fill = "none", color = "none") +
  # if scale starting from 0 is desired
  #scale_y_continuous(limits = c(0, 200)) +
  scale_fill_manual(values = c("#1B9E77", "#0072B2")) +
  scale_color_manual(values = c("#1B9E77", "#0072B2")) +
  stat_compare_means(comparisons = list(c("Ecuador","Col. & Ven.")),
                     #method = "wilcox.test", paired = FALSE, # these are default
                     size = 6, label = "p.signif") +
  theme_clim()
suc.char.plot7

### Character 8. Ventral Opening - not significant
names(data.SUC)
# select only columns of interest and remove missing samples
suc.ventral.len <- na.omit(data.SUC[,c("sample_num", "country", "cortube_vent_len")])

# make the plot
suc.char.plot8 <- 
  ggplot(data = suc.ventral.len, aes(country, cortube_vent_len, fill = country)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.7) +
  geom_point(aes(y = cortube_vent_len, color = country), show.legend = FALSE,
             position = position_jitter(width = 0.15), size = 2, alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  labs(y = "Length to Ventral Corolla Opening (mm)", x = NULL) +
  guides(fill = "none", color = "none") +
  # if scale starting from 0 is desired
  #scale_y_continuous(limits = c(0, 200)) +
  scale_fill_manual(values = c("#1B9E77", "#0072B2")) +
  scale_color_manual(values = c("#1B9E77", "#0072B2")) +
  stat_compare_means(comparisons = list(c("Ecuador","Col. & Ven.")),
                     #method = "wilcox.test", paired = FALSE, # these are default
                     size = 6, label = "p.signif") +
  theme_clim()
suc.char.plot8

### Character 9. Androecium Length
names(data.SUC)
# select only columns of interest and remove missing samples
suc.andr.len <- na.omit(data.SUC[,c("sample_num", "country", "andr_len")])

# make the plot
suc.char.plot9 <- 
  ggplot(data = suc.andr.len, aes(country, andr_len, fill = country)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.7) +
  geom_point(aes(y = andr_len, color = country), show.legend = FALSE,
             position = position_jitter(width = 0.15), size = 2, alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  labs(y = "Androecium Length (mm)", x = NULL) +
  guides(fill = "none", color = "none") +
  # if scale starting from 0 is desired
  #scale_y_continuous(limits = c(0, 200)) +
  scale_fill_manual(values = c("#1B9E77", "#0072B2")) +
  scale_color_manual(values = c("#1B9E77", "#0072B2")) +
  stat_compare_means(comparisons = list(c("Ecuador","Col. & Ven.")),
                     #method = "wilcox.test", paired = FALSE, # these are default
                     size = 6, label = "p.signif") +
  theme_clim()
suc.char.plot9

### Character 10. Exsertion Length
names(data.SUC)
# select only columns of interest and remove missing samples
suc.exsertion.len <- na.omit(data.SUC[,c("sample_num", "country", "exsertion_len")])

# make the plot
suc.char.plot10 <- 
  ggplot(data = suc.exsertion.len, aes(country, exsertion_len, fill = country)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.7) +
  geom_point(aes(y = exsertion_len, color = country), show.legend = FALSE,
             position = position_jitter(width = 0.15), size = 2, alpha = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  labs(y = "Exsertion Length (mm)", x = NULL) +
  guides(fill = "none", color = "none") +
  # if scale starting from 0 is desired
  #scale_y_continuous(limits = c(0, 200)) +
  scale_fill_manual(values = c("#1B9E77", "#0072B2")) +
  scale_color_manual(values = c("#1B9E77", "#0072B2")) +
  stat_compare_means(comparisons = list(c("Ecuador","Col. & Ven.")),
                     #method = "wilcox.test", paired = FALSE, # these are default
                     size = 6, label = "p.signif") +
  theme_clim()
suc.char.plot10

### Figure 4
## export individual plots for final figures
#pedicel length
ggsave(suc.char.plot3, device = "svg", filename = "F4-PedicelLength.svg",
       height = 6, width = 6, dpi = 600)
#hypanthium length
ggsave(suc.char.plot4, device = "svg", filename = "F4-HypanthiumLength.svg",
       height = 6, width = 6, dpi = 600)
#calyx lobe length
ggsave(suc.char.plot6, device = "svg", filename = "F4-CalyxLobeLength.svg",
       height = 6, width = 6, dpi = 600)
#calyx lobe width
ggsave(suc.char.plot7, device = "svg", filename = "F4-CalyxLobeWidth.svg",
       height = 6, width = 6, dpi = 600)
#androecium length
ggsave(suc.char.plot9, device = "svg", filename = "F4-AndroeciumLength.svg",
       height = 6, width = 6, dpi = 600)
#exsertion length
ggsave(suc.char.plot10, device = "svg", filename = "F4-ExsertionLength.svg",
       height = 6, width = 6, dpi = 600)

### Figure S3
## make a figure of the remaining characters compared
# Make a panel
multi_plot3 <- ggarrange(suc.char.plot1,suc.char.plot2,
                         suc.char.plot5,suc.char.plot8,
                         labels = c("A","B","C","D"), #labels given each panel 
                         font.label = list(size = 26, color = "black", family = NULL),
                         ncol = 2, nrow = 2, #adjust plot space 
                         common.legend = FALSE) #does the plot have a common legend
#add titles and labels to the multi-panel graph
multi_plot3 <- annotate_figure(multi_plot3,
                               top = text_grob("Figure S3. Phase 2 Univariate Analysis: Additional Traits Measured", 
                                               color = "black", face = "bold", size = 26))
# save Figure S3
ggsave(multi_plot3, device = "tiff", filename = "FigureS3.tiff", 
       height = 12, width = 12, dpi = 350, bg = "white")
ggsave(multi_plot3, device = "pdf", filename = "FigureS3.pdf",
       height = 12, width = 12, dpi = 350, bg = "white")

###########
# Phase 2B: SUC Multivariate Analysis
###########

names(data.SUC)

##run PCA
SUC.pca <- PCA(data.SUC[ , -c(1:6,16)], ncp = 6, graph = FALSE)
SUC.pca
SUC.pca$eig
SUC.pca$var$cos2
SUC.pca$var$contrib

# PCA with sample labels
SUC.PCA.1 <- fviz_pca_ind(SUC.pca, 
                        geom.ind = c("point", "text"), pointsize = 2.5,
                        mean.point = FALSE, repel = TRUE,
                        col.ind = data.SUC$country)
SUC.PCA.1 <- ggpubr::ggpar(SUC.PCA.1,
                              title = "Phase 2 PCA: Sample Numbers",
                              font.main = c(18, "bold", "black"),
                              font.x = 12, font.y = 12,
                              legend = c(0.88, 0.1),
                              font.legend = 12,
                              ggtheme = theme_bw(),
                              palette = c("#1B9E77", "#0072B2")) +
  theme(legend.title = element_blank())
SUC.PCA.1

#visualize contribution of variables to distributions
SUC.PCA.2 <- fviz_pca_biplot(SUC.pca, geom.ind = "point", fill.ind = "white", 
                         label = "var", col.var = "black",
                         mean.point = FALSE,
                         col.ind = data.SUC$country, alpha.ind = 0.8, 
                         pointsize = 3, repel = TRUE)
SUC.PCA.2 <- ggpubr::ggpar(SUC.PCA.2, title = "PCA of B. succulenta s.l.: Quant. Trait Contributions",
                       font.x = 12, font.y = 12,
                       legend = "none",
                       ggtheme = theme_bw(), 
                       palette = c("#1B9E77", "#0072B2")) +
  theme(legend.title = element_blank())
SUC.PCA.2

#heirarchical clustering unsupervised
hc.SUC.pca <- HCPC(SUC.pca, nb.clust = -1)
SUC.PCA.3 <- fviz_cluster(hc.SUC.pca, repel = TRUE, pointsize = 3,
                      show.clust.cent = FALSE, geom = "point",
                      ggtheme = theme_bw(), palette = "jco")
SUC.PCA.3 <- ggpubr::ggpar(SUC.PCA.3, title = "PCA of B. succulenta s.l.: Clustering Results",
                       font.x = 12, font.y = 12,
                       legend = "none"
                       )
SUC.PCA.3
hc.SUC.pca$desc.var$quanti
hc.SUC.pca$data.clust$clust

#heirarchical clustering unsupervised - with sample numbers
SUC.PCA.4 <- fviz_cluster(hc.SUC.pca, repel = TRUE, pointsize = 2.5,
                          show.clust.cent = FALSE, geom = c("point","text"),
                          ggtheme = theme_bw(), palette = "jco")
SUC.PCA.4 <- ggpubr::ggpar(SUC.PCA.4, title = "Phase 2 PCA: Clusters with Sample Numbers",
                           font.main = c(18, "bold", "black"),
                           font.x = 12, font.y = 12,
                           legend = c(0.92, 0.16),
                           font.legend = 12,
                           ggtheme = theme_bw()
                           )
SUC.PCA.4


# save SUC PCA clustering results to data frame
SUC.cluster.res$SUC_pca_clust <- hc.SUC.pca$data.clust$clust

###Factor Analysis of Mixed Data on SUC
names(data.SUC)
SUC.famd <- FAMD(data.SUC[ , -c(1:5)], ncp = 6, graph = FALSE)

# FAMD plot with labels
SUC.FAMD.1 <- fviz_famd_ind(SUC.famd, col.ind = data.SUC$country,
                          mean.point = FALSE, geom = c("point","text"),
                          pointsize = 2.5, 
                          invisible = "quali.var",
                          repel = TRUE)
SUC.FAMD.1 <- ggpubr::ggpar(SUC.FAMD.1,
                   title = "Phase 2 FAMD: Sample Numbers",
                   font.main = c(18, "bold", "black"),
                   font.x = 12, font.y = 12,
                   legend = c(0.88, 0.1),
                   font.legend = 12,
                   ggtheme = theme_bw(), 
                   palette = c("#1B9E77", "#0072B2")) + 
  theme(legend.title = element_blank())
SUC.FAMD.1

# visualize the contribution of variabls
SUC.FAMD.2 <- fviz_famd_ind(SUC.famd, 
                       habillage = as.factor(data.SUC$country),
                       col.quali.var = "black",
                       #goem = c("point", "text"),
                       label = "quali.var",
                       mean.point = FALSE,
                       alpha.ind = 1, repel = TRUE)
SUC.FAMD.2 <- ggpubr::ggpar(SUC.FAMD.2, title = "FAMD of B. succulenta s.l.: Qual. Trait Contributions",
                       font.x = 12, font.y = 12,
                       legend = "none",
                       ggtheme = theme_bw(),
                       palette = c("#1B9E77", "#0072B2")) + 
  theme(legend.title = element_blank())
SUC.FAMD.2

#heirarchical clustering unsupervised
hc.SUC.famd <- HCPC(SUC.famd, nb.clust = -1)
SUC.FAMD.3 <- fviz_cluster(hc.SUC.famd, repel = TRUE, geom = "point",
                      show.clust.cent = FALSE, pointsize = 3,
                      ggtheme = theme_bw(), palette = "jco")
SUC.FAMD.3 <- ggpubr::ggpar(SUC.FAMD.3, title = "FAMD of B. succulenta s.l.: Clustering Results",
                       font.x = 12, font.y = 12,
                       legend = "none"
                       )
SUC.FAMD.3
hc.SUC.famd$desc.var$quanti
hc.SUC.famd$data.clust$clust

#heirarchical clustering unsupervised - with sample labels
SUC.FAMD.4 <- fviz_cluster(hc.SUC.famd, repel = TRUE, geom = c("point","text"),
                           show.clust.cent = FALSE, pointsize = 2.5,
                           ggtheme = theme_bw(), palette = "jco")
SUC.FAMD.4 <- ggpubr::ggpar(SUC.FAMD.4, title = "Phase 2 FAMD: Clusters with Sample Numbers",
                            font.title = c(18, "bold", "black"), 
                            font.x = 12, font.y = 12,
                            legend = c(0.92, 0.14),
                            font.legend = c(12)
                            )
SUC.FAMD.4

# save SUC FAMD clustering results to data frame
SUC.cluster.res$SUC_famd_clust <- hc.SUC.famd$data.clust$clust

##### Model Based Clustering on SUC
names(data.SUC)
data.SUC.norm <- scale(data.SUC[, -c(1:6,16)])
mc.SUC <- Mclust(data.SUC.norm)
summary(mc.SUC)
mc.SUC$classification
# what about 2 clusters?
mc.SUC.2 <- Mclust(data.SUC.norm, G = 2)
fviz_mclust(mc.SUC.2, "classification", 
            ggtheme = theme_bw(), palette = "jco")
summary(mc.SUC.2)

# save SUC FAMD clustering results to data frame
SUC.cluster.res$SUC_mbc_2clust <- mc.SUC.2$classification
SUC.cluster.res$SUC_mbc_3clust <- mc.SUC$classification

# plot the BIC scores
SUC.MBC.1 <- fviz_mclust(mc.SUC, "BIC",
                     ggtheme = theme_bw(), palette = "jco")
SUC.MBC.1 <- ggpubr::ggpar(SUC.MBC.1, title = "Phase 2 Model-based Clust.: Model Selection",
                   font.main = c(18, "bold", "black"),
                   font.x = 12, font.y = 12)
SUC.MBC.1

# plot 2 clusters in clean format
SUC.MBC.2 <- fviz_mclust(mc.SUC.2, "classification",
                     show.clust.cent = FALSE,
                     geom = "point", pointsize = 3,
                     repel = TRUE, ellipse.type = "convex",
                     ggtheme = theme_bw(), palette = "jco")
SUC.MBC.2 <- ggpubr::ggpar(SUC.MBC.2, title = "Model-based Clustering of B. succulenta s.l.: K = 2",
                   font.x = 12, font.y = 12,
                   legend = "none",
                   ggtheme = theme_bw()) +
  theme(plot.subtitle = element_blank())
SUC.MBC.2

# plot 3 clusters with sample numbers included
SUC.MBC.3 <- fviz_mclust(mc.SUC, "classification",
                         show.clust.cent = FALSE,
                         geom = c("point","text"), pointsize = 2.5,
                         repel = TRUE, ellipse.type = "convex",
                         ggtheme = theme_bw(), palette = "jco")
SUC.MBC.3 <- ggpubr::ggpar(SUC.MBC.3, title = "Phase 2 Model-based Clust.: n = 3, Sample Numbers",
                           font.main = c(18, "bold", "black"),
                           font.x = 12, font.y = 12,
                           legend = c(0.08, 0.14),
                           font.legend = 12,
                           ggtheme = theme_bw()) +
  theme(plot.subtitle = element_blank())
SUC.MBC.3

# plot 3 clusters in clean format
SUC.MBC.4 <- fviz_mclust(mc.SUC, "classification",
                 show.clust.cent = FALSE,
                 geom = "point", pointsize = 3,
                 repel= TRUE, ellipse.type = "convex",
                 ggtheme = theme_bw(), palette = "jco")
SUC.MBC.4 <- ggpubr::ggpar(SUC.MBC.4, title = "Model-based Clustering of B. succulenta s.l.: K = 3", 
                   font.x = 12, font.y = 12,
                   legend = "none",
                   ggtheme = theme_bw()) + 
  theme(plot.subtitle = element_blank())
SUC.MBC.4

### save the clustering results as a csv
#write.csv(SUC.cluster.res, "SUC_cluster_results.csv")

### Figure 5
## export individual plots for final figures
#SUC - PCA with variables plotted
#ggsave(SUC.PCA.2, device = "svg", filename = "F5-A.svg",
#       height = 6, width = 6, dpi = 600)
#RCC - PCA with clusters
#ggsave(SUC.PCA.3, device = "svg", filename = "F5-B.svg",
#       height = 6, width = 6, dpi = 600)
##RCC - FAMD with variables plotted
#ggsave(SUC.FAMD.2, device = "svg", filename = "F5-C.svg",
#       height = 6, width = 6, dpi = 600)
##RCC - FAMD with clusters
#ggsave(SUC.FAMD.3, device = "svg", filename = "F5-D.svg",
#       height = 6, width = 6, dpi = 600)
##RCC - MClust with 2 clusters
#ggsave(SUC.MBC.2, device = "svg", filename = "F5-E.svg",
#       height = 6, width = 6, dpi = 600)
##RCC - MClust with 3 clusters
#ggsave(SUC.MBC.4, device = "svg", filename = "F5-F.svg",
#       height = 6, width = 6, dpi = 600)


### Figure S4
multi_plot4 <- ggarrange(SUC.PCA.1,SUC.PCA.4,
                         SUC.FAMD.1,SUC.FAMD.4,
                         SUC.MBC.1,SUC.MBC.3,
                         labels = c("A","B","C","D","E","F"), #labels given each panel
                         font.label = list(size = 30, color = "black", family = NULL),
                         ncol = 2, nrow = 3, #adjust plot space 
                         common.legend = FALSE) #does the plot have a common legend
#add titles and labels to the multi-panel graph
multi_plot4 <- annotate_figure(multi_plot4,
                               top = text_grob("Figure S4. Phase 2 Multivariate Analysis", 
                                               color = "black", face = "bold", size = 26))
# save the final plot
#ggsave(multi_plot4, device = "tiff", filename = "FigureS4.tiff", 
#       height = 19, width = 14, dpi = 350, bg = "white")
#ggsave(multi_plot4, device = "pdf", filename = "FigureS4.pdf", 
#       height = 19, width = 14, dpi = 350, bg = "white")

################################################################################
##PHASE 3: Analysis of the B. sodiroana / B. crispiloba group (SCC)
################################################################################

names(data.all)
data.all$pile_id1

# keep the SCC samples
data.SCC <- data.all[c(1:68),c(1:6,9:11,14:32)]
names(data.SCC)

# make a data frame to hold clustering results
SCC.cluster.res <- na.omit(data.SCC[,c(1:6)])

##run PCA
SCC.pca <- PCA(data.SCC[ , -c(1:7,25)], ncp = 6, graph = FALSE)
SCC.pca
SCC.pca$eig
SCC.pca$var$cos2
SCC.pca$var$contrib

SCC.PCA.1 <- fviz_pca_ind(SCC.pca,
                          geom.ind = c("point", "text"), pointsize = 2.5,
                          mean.point = FALSE, repel = TRUE,
                          col.ind = data.SCC$pile_id2)
SCC.PCA.1 <- ggpubr::ggpar(SCC.PCA.1,
                       title = "Phase 3 PCA: Sample Numbers",
                       font.main = c(18, "bold", "black"),
                       font.x = 12, font.y = 12,
                       legend = c(0.875, 0.1),
                       font.legend = c(12, "italic"),
                       ggtheme = theme_bw(),
                       palette = c("#1B9E77", "#7570B3", "#FFB000")) +
  theme(legend.title = element_blank())
SCC.PCA.1

# visualize the variables driving PCA divergence
SCC.PCA.2 <- fviz_pca_biplot(SCC.pca,
                         fill.ind = "white", label = "var", col.var = "black",
                         col.ind = data.SCC$pile_id2, alpha.ind = 1,
                         mean.point = FALSE,
                         pointsize = 3, repel = TRUE)
SCC.PCA.2 <- ggpubr::ggpar(SCC.PCA.2, title = "Phase 3 PCA: Quantitative Characters",
                   font.x = 12, font.y = 12,
                   legend = "none",
                   ggtheme = theme_bw(), 
                   palette = c("#1B9E77", "#7570B3", "#FFB000")) +
  theme(legend.title = element_blank())
SCC.PCA.2

#heirarchical clustering unsupervised
hc.SCC.pca <- HCPC(SCC.pca, nb.clust = -1)
SCC.PCA.3 <- fviz_cluster(hc.SCC.pca, repel = TRUE, pointsize = 3,
                      show.clust.cent = FALSE, geom = "point", alpha.ind = 1,
                      ggtheme = theme_bw(), palette = "jco")
SCC.PCA.3 <- ggpubr::ggpar(SCC.PCA.3, title = "Phase 3 PCA: Clusters",
                   font.x = 12, font.y = 12,
                   legend = "none"
                   )
SCC.PCA.3
hc.SCC.pca$desc.var$quanti
hc.SCC.pca$data.clust$clust

#heirarchical clustering unsupervised - with sample numbers
hc.SCC.pca <- HCPC(SCC.pca, nb.clust = -1)
SCC.PCA.4 <- fviz_cluster(hc.SCC.pca, repel = TRUE, pointsize = 2.5,
                          show.clust.cent = FALSE, geom = c("point","text"),
                          ggtheme = theme_bw(), palette = "jco")
SCC.PCA.4 <- ggpubr::ggpar(SCC.PCA.4, title = "Phase 3 PCA: Clusters with Sample Numbers",
                           font.main = c(18, "bold", "black"),
                           legend = c(0.92, 0.16),
                           font.legend = c(12),
                           font.x = 12, font.y = 12,
                           ggtheme = theme_bw(),
                           )
SCC.PCA.4

# save SUC PCA clustering results to data frame
SCC.cluster.res$SCC_pca_clust <- hc.SCC.pca$data.clust$clust

###Factor Analysis of Mixed Data on SCC
names(data.SCC)
SCC.famd <- FAMD(data.SCC[,-c(1:6)], ncp = 6, graph = FALSE)

# FAMD plot with labels
SCC.FAMD.1 <- fviz_famd_ind(SCC.famd, habillage = as.factor(data.SCC$pile_id2),
                       mean.point = FALSE,
                       geom = c("point","text"), pointsize = 2.5, repel = TRUE,
                       invisible = "quali.var")
SCC.FAMD.1 <- ggpubr::ggpar(SCC.FAMD.1,
                   title = "Phase 3 FAMD: Sample Numbers",
                   font.main = c(18, "bold", "black"),
                   font.x = 12, font.y = 12,
                   legend = c(0.12, 0.9),
                   font.legend = c(12,"italic"),
                   ggtheme = theme_bw(), 
                   palette = c("#1B9E77", "#7570B3", "#FFB000")) + 
  theme(legend.title = element_blank())
SCC.FAMD.1

# look at the contribution of variables
SCC.FAMD.2 <- fviz_famd_ind(SCC.famd, habillage = as.factor(data.SCC$pile_id2),
                   col.quali.var = "black", 
                   #goem = c("point", "text"),
                   label = "quali.var",
                   mean.point = FALSE, 
                   alpha.ind = 1, repel = TRUE)
SCC.FAMD.2 <- ggpubr::ggpar(SCC.FAMD.2, title = "Phase 3 FAMD: Qualitative Variables",
                   font.x = 12, font.y = 12,
                   legend = "none",
                   ggtheme = theme_bw(), 
                   palette = c("#1B9E77", "#7570B3", "#FFB000")) + 
  theme(legend.title = element_blank())
SCC.FAMD.2

#heirarchical clustering unsupervised
hc.SCC.famd <- HCPC(SCC.famd, nb.clust = -1)
SCC.FAMD.3 <- fviz_cluster(hc.SCC.famd, repel = TRUE, geom = "point", 
                      show.clust.cent = FALSE, pointsize = 3,
                      ggtheme = theme_bw(), palette = "jco")
SCC.FAMD.3 <- ggpubr::ggpar(SCC.FAMD.3, title = "Phase 3 FAMD: Clusters",
                   font.x = 12, font.y = 12,
                   legend = "none")
SCC.FAMD.3
hc.SCC.famd$desc.var$quanti
hc.SCC.famd$data.clust$clust
#SCC.dend.famd <- fviz_dend(hc.SCC.famd, pallette = "locuszoom", show_labels = FALSE, ggtheme = theme_grey())
#SCC.dend.famd

#heirarchical clustering unsupervised with sample numbers
SCC.FAMD.4 <- fviz_cluster(hc.SCC.famd, repel = TRUE, geom = c("point","text"), 
                           show.clust.cent = FALSE, pointsize = 2.5,
                           ggtheme = theme_bw(), palette = "jco")
SCC.FAMD.4 <- ggpubr::ggpar(SCC.FAMD.4, title = "Phase 3 FAMD: Clusters with Sample Numbers",
                            font.title = c(18, "bold", "black"), 
                            font.x = 12, font.y = 12,
                            legend = c(0.08, 0.88),
                            font.legend = 12)
SCC.FAMD.4

# save SUC FAMD clustering results to data frame
SCC.cluster.res$SCC_famd_clust <- hc.SCC.famd$data.clust$clust

##### Model Based Clustering on SCC
names(data.SCC)
data.SCC.norm <- scale(data.SCC[, -c(1:7,25)])
mc.SCC <- Mclust(data.SCC.norm)
summary(mc.SCC)
mc.SCC$classification

# save SUC FAMD clustering results to data frame
SCC.cluster.res$SCC_mbc_clust <- mc.SCC$classification

#plot the results
SCC.MBC.1 <- fviz_mclust(mc.SCC, "BIC", 
                 ggtheme = theme_bw(), palette = "jco")
SCC.MBC.1 <- ggpubr::ggpar(SCC.MBC.1, title = "Phase 3 Model-based Clust.: Model Selection",
                   font.main = c(16, "bold", "black"),
                   font.x = 12, font.y = 12)
SCC.MBC.1

# plot two clusters in clean format
mc.SCC.2 <- Mclust(data.SCC.norm, G = 2)
SCC.MBC.2 <- fviz_mclust(mc.SCC.2, "classification",
                     show.clust.cent = FALSE,
                     geom = "point", pointsize = 3,
                     repel = TRUE, ellipse.type = "convex",
                     ggtheme = theme_bw(), palette = "jco")
SCC.MBC.2 <- ggpubr::ggpar(SCC.MBC.2, title = "Phase 3 Model-based Clustering: K = 2",
                   font.x = 12, font.y = 12,
                   legend = "none",
                   ggtheme = theme_bw()) +
  theme(plot.subtitle = element_blank()) +
  scale_y_reverse()
SCC.MBC.2

SCC.cluster.res$SCC_mbc_k2 <- mc.SCC.2$classification

# plot two clusters with sample numbers included
SCC.MBC.3 <- fviz_mclust(mc.SCC.2, "classification",
                         show.clust.cent = FALSE,
                         geom = c("point","text"), pointsize = 2.5,
                         repel = TRUE, ellipse.type = "convex",
                         ggtheme = theme_bw(), palette = "jco")
SCC.MBC.3 <- ggpubr::ggpar(SCC.MBC.3, title = "Phase 3 Model-based Clust.: n = 2, Sample Numbers",
                           font.main = c(18, "bold", "black"),
                           font.x = 12, font.y = 12,
                           legend = c(0.08, 0.1),
                           font.legend = 12,
                           ggtheme = theme_bw()) +
  theme(plot.subtitle = element_blank()) +
  scale_y_reverse()
SCC.MBC.3

# plot three clusters in clean format
mc.SCC.3 <- Mclust(data.SCC.norm, G = 3)
SCC.MBC.4 <- fviz_mclust(mc.SCC.3, "classification",
                     show.clust.cent = FALSE,
                     geom = "point", pointsize = 3,
                     repel = TRUE, ellipse.type = "convex",
                     ggtheme = theme_bw(), palette = "jco")
SCC.MBC.4 <- ggpubr::ggpar(SCC.MBC.4, title = "Phase 3 Model-based Clustering: K = 3",
                       font.x = 12, font.y = 12,
                       legend = "none",
                       ggtheme = theme_bw()) +
  theme(plot.subtitle = element_blank()) +
  scale_y_reverse()
SCC.MBC.4

SCC.cluster.res$SCC_mbc_k3 <- mc.SCC.3$classification

### save the clustering results as a csv
write.csv(SCC.cluster.res, "SCC_cluster_results.csv")

### Figure 6
## export individual plots to make final Figure 6 in inkscape
#SCC - PCA with variables plotted
ggsave(SCC.B, device = "svg", filename = "F6-A.svg",
       height = 6, width = 6, dpi = 600)
#SCC - PCA with clusters
ggsave(SCC.C, device = "svg", filename = "F6-B.svg",
       height = 6, width = 6, dpi = 600)
##SCC - FAMD with variables plotted
ggsave(SCC.E, device = "svg", filename = "F6-C.svg",
       height = 6, width = 6, dpi = 600)
##SCC - FAMD with clusters
ggsave(SCC.F, device = "svg", filename = "F6-D.svg",
       height = 6, width = 6, dpi = 600)
##SCC - MClust with 2 clusters
ggsave(SCC.H, device = "svg", filename = "F6-E.svg",
       height = 6, width = 6, dpi = 600)
##SCC - MClust with 3 clusters
ggsave(SCC.I, device = "svg", filename = "F6-F.svg",
       height = 6, width = 6, dpi = 600)

### Figure S5
multi_plot5 <- ggarrange(SCC.PCA.1,SCC.PCA.4,
                         SCC.FAMD.1,SCC.FAMD.4,
                         SCC.MBC.1,SCC.MBC.3,
                         labels = c("A","B","C","D","E","F"), #labels given each panel
                         font.label = list(size = 30, color = "black", family = NULL),
                         ncol = 2, nrow = 3, #adjust plot space 
                         common.legend = FALSE) #does the plot have a common legend
#add titles and labels to the multi-panel graph
multi_plot5 <- annotate_figure(multi_plot5,
                               top = text_grob("Figure S5. Phase 3 Multivariate Analysis", 
                                               color = "black", face = "bold", size = 26))
# save the final plot
ggsave(multi_plot5, device = "tiff", filename = "FigureS5.tiff", 
       height = 19, width = 14, dpi = 350, bg = "white")
#ggsave(multi_plot5, device = "pdf", filename = "FigureS5.pdf", 
#       height = 19, width = 14, dpi = 350, bg = "white")

#####
#END#
#####