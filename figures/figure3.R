# these are all run from project parent directory
library(ArchR)
suppressMessages(library(glue))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(ggrastr))
library(ggtext)
library(ggpubr)
library(Seurat)

load("rdata/chronicSub2_deviations_2025-02-26-01-56.Rdata")


source("figures/defaultFigureSettings.R")
source("scripts/graphs.R")
source("scripts/colorPalettes.R")

projUnhashedClean$stage <- ifelse(projUnhashedClean$stage == "Chronic", "Viremic", projUnhashedClean$stage)
seuMerged$stage <- ifelse(seuMerged$stage == "Chronic", "Viremic", seuMerged$stage)

tmp <- FindMarkers(seuMerged, ident.1 = "CD4 GC-Tfh #1", ident.2 = "CD4 GC-Tfh #2", min.pct = 0.2, assay = "ADT")


data.frame(
  stage = projUnhashedClean$stage,
  donor = projUnhashedClean$donorClean2,
  hivStatus = projUnhashedClean$hivStatus,
  hivDna = projUnhashedClean$haystackOut,
  manualAnnot = projUnhashedClean$manualAnnot) %>%
  dplyr::group_by(stage, manualAnnot, hivDna) %>%
  dplyr::tally(name = "nCells") %>%
  dplyr::group_by(stage, manualAnnot) %>%
  dplyr::mutate(propOfCluster = nCells / sum(nCells)) %>%
  dplyr::filter(hivDna) %>%
  ggplot(aes(x = propOfCluster, y = manualAnnot)) +
  geom_point() +
  facet_wrap(~ stage, scales = "free_x") +
  theme_bw()


data.frame(
  stage = projUnhashedClean$stage,
  donor = projUnhashedClean$donorClean2,
  hivStatus = projUnhashedClean$hivStatus,
  hivDna = projUnhashedClean$haystackOut,
  manualAnnot = projUnhashedClean$manualAnnot) %>%
  dplyr::group_by(donor, stage, manualAnnot) %>%
  dplyr::tally(name = "nCells") %>%
  dplyr::group_by(donor) %>%
  dplyr::mutate(propOfDonor = nCells / sum(nCells)) %>%
  dplyr::filter(manualAnnot == "CD4 PRDM1+ GC-Tfh") %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(propOfDonor)) %>%
  dplyr::mutate(donor = factor(donor, levels = unique(donor))) %>%
  ggplot(aes(x = propOfDonor * 100, y = donor, color = stage)) +
  labs(x = "% PRDM1+ GC-Tfh in donor", y = "Donor") +
  geom_point() +
  theme_bw() +
  theme(axis.text = element_text(color = "#000000"),
    legend.position = "bottom")





layout <- c(
  area(1, 1, 2, 4), # umaps
  area(1, 5, 2, 6), # upset plots
  area(1, 7, 2, 9), # distance
  area(3, 1, 5, 6), # dna lol
  area(6, 1, 8, 6), # rna lol
  area(3, 7, 4, 9), # tfh comp
  area(5, 7, 6, 9), # tfh corr
  area(7, 7, 8, 9) # gonzalo data?
  # area(5, 1, 8, 9) #d
)

p <- wrap_elements(full = figA) +
  wrap_elements(full = figB) +
  wrap_elements(full = figDist) +
  wrap_elements(full = figDnaLol) +
  wrap_elements(full = figRnaLol) +
  wrap_elements(full = figTfh) + theme(plot.tag = element_text(margin = margin(t = -75))) +
  wrap_elements(full = figTfhCorr) +
  wrap_elements(panel = plot_spacer()) +
  plot_annotation(tag_levels = "A") +
  plot_layout(design = layout)

saveFinalFigure(plot = p, fn = "fig3_v1", prefixDir = "figures/outs", devices = c("png"), gwidth = 8.5, gheight = 10)
