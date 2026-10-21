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
library(magick)

load("rdata/chronicSub2_deviations_2025-02-26-01-56.Rdata")

set.seed(42)

source("figures/defaultFigureSettings.R")
source("scripts/graphs.R")
source("scripts/colorPalettes.R")

projUnhashedClean$stage <- ifelse(projUnhashedClean$stage == "Chronic", "Viremic", projUnhashedClean$stage)
seuMerged$stage <- ifelse(seuMerged$stage == "Chronic", "Viremic", seuMerged$stage)
projUnhashedClean$manualAnnot <- trimws(projUnhashedClean$manualAnnot)


uniqClusters <- unique(projUnhashedClean$manualAnnot)
uniqStages <- c("Viremic", "ART")
uniqInfection <- names(COLORS$hivStatus)


################################################################################
# panel describing blimp1 cells
################################################################################
selectGenesForAtac <- c("PDCD1", "CXCR5", "BCL6", "PRDM1", "CCR5", "TGFB1", "FOXP3", "IL2RA")
atacGeneScores <- tryCatch({
  tmp <- ArchR:::.getMatrixValues(
    ArchRProj = projUnhashedClean,
    name = selectGenesForAtac,
    matrixName = "GeneScoreMatrix",
    log2Norm = FALSE,
    threads = getArchRThreads()
  )
},
  error = function(cond) {
    message(conditionMessage(cond))
    # Choose a return value in case of error
    return(NULL)
  })

atacGeneScoresImputed <- imputeMatrix(
  mat = as.matrix(atacGeneScores),
  imputeWeights = getImputeWeights(projUnhashedClean))

atacGeneScoresImputed <- log2(atacGeneScoresImputed + 1)

atacImputedDf <- data.frame(
  cbc = projUnhashedClean$cellNames,
  manualAnnot = projUnhashedClean$manualAnnot,
  hivPos = projUnhashedClean$hivPos
)
rownames(atacImputedDf) <- atacImputedDf$cbc
atacImputedDf[colnames(atacGeneScoresImputed), selectGenesForAtac] <- t(atacGeneScoresImputed)


clustersOfInterest <- c("CD4 PRDM1+ GC-Tfh", "CD4 GC-Tfh #1", "CD4 GC-Tfh #2", "CD4 Treg #1", "CD4 Treg #2")
atacImputedDf <- atacImputedDf %>%
  pivot_longer(cols = -any_of(c("cbc", "manualAnnot", "hivPos")), names_to = "feature", values_to = "exp") %>%
  dplyr::filter(manualAnnot %in% clustersOfInterest)


atacImputedDf %>%
  dplyr::mutate(feature = factor(feature, levels = selectGenesForAtac)) %>%
  dplyr::mutate(newY = paste0(hivPos, "_", manualAnnot)) %>%
  ggplot(aes(y = newY, x = exp)) +
  ggridges::geom_density_ridges(rel_min_height = 0.01, scale = 0.75, alpha = 0.8, aes(fill = manualAnnot)) + 
  facet_wrap( ~ feature, scales = "free_x", ncol = 4)


atacFigA <- atacImputedDf %>%
  dplyr::mutate(feature = factor(feature, levels = selectGenesForAtac)) %>%
  ggplot(aes(y = manualAnnot, x = exp)) +
  ggridges::geom_density_ridges(rel_min_height = 0.01, scale = 0.75, alpha = 0.8, aes(fill = manualAnnot)) +
  geom_vline(xintercept = -Inf, linewidth = BASELINEWIDTHSIZE, color = "#000000") +
  geom_hline(yintercept = -Inf, linewidth = BASELINEWIDTHSIZE, color = "#000000") +
  facet_wrap(~ feature, scales = "free_x", ncol = 4) +
  coord_cartesian(clip = "off") +
  scale_fill_brewer(palette = "Set1") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = BASEPTFONTSIZE - 2, color = "#000000", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = BASEPTFONTSIZE - 2, color = "#000000"),
    strip.text = element_text(size = BASEPTFONTSIZE - 1, color = "#000000", face = "italic"),
    strip.background = element_rect(color = "#00000000", fill = "#00000000"),
    legend.position = "blank",
    panel.grid.major = element_line(color = "#cfcfcf50"),
    axis.line = element_blank(),
    axis.title.x = element_text(size = BASEPTFONTSIZE - 1),
    axis.title.y = element_blank()) +
  labs(
    x = "imputed log2(counts + 1)")

rnaFigB <- DotPlot(
  subset(seuMerged, subset = manualAnnot %in% clustersOfInterest),
  features = selectGenesForAtac,
  group.by = "manualAnnot")

rnaFigB_mod <- rnaFigB +
  theme(
    axis.text = element_text(size = BASEPTFONTSIZE - 2, color = "#000000"),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = BASEPTFONTSIZE - 1, face = "italic"),
    legend.title = element_text(size = BASEPTFONTSIZE - 2, color = "#000000", hjust = 0.5),
    legend.text = element_text(size = BASEPTFONTSIZE - 2, color = "#000000", hjust = 0),
    panel.grid.major.y = element_line(color = "#cfcfcf30"),
    axis.line = element_line(color = "#000000", linewidth = BASELINEWIDTHSIZE),
    legend.margin = margin(l = 0),
    legend.key.size = unit(.5, "lines")
  )

rnaFigB_mod$guides$guides$colour$params$title <- "Average\n Expression"
rnaFigB_mod$guides$guides$size$params$title <- "% Expressed"


################################################################################
# all HIV+ cells
################################################################################
propClusterDf <- data.frame(
  donor = projUnhashedClean$donorClean2,
  stage = projUnhashedClean$stage,
  infection = projUnhashedClean$hivStatus,
  manualAnnot = projUnhashedClean$manualAnnot
) %>%
  dplyr::filter(stage != "HIV-") %>%
  dplyr::mutate(
    manualAnnot = factor(manualAnnot, levels = uniqClusters),
    stage = factor(stage, levels = uniqStages),
    infection = factor(infection, levels = uniqInfection)
  ) %>%
  dplyr::group_by(stage, manualAnnot, infection, .drop = FALSE) %>%
  dplyr::tally(name = "nCells") %>%
  dplyr::group_by(stage, manualAnnot) %>%
  dplyr::mutate(
    totalCells = sum(nCells),
    percent = nCells / totalCells * 100,
    percent = ifelse(is.nan(percent), 0, percent),
    percentLabel = paste0(round(percent, digits = 2), "% (n=", nCells, ")"),
    manualAnnotLabel = paste0(manualAnnot, " (", totalCells, ")")) %>%
  dplyr::filter(infection != "HIV-") %>%
  dplyr::group_by(stage, manualAnnot) %>%
  dplyr::mutate(avgPercent = mean(percent)) %>%
  dplyr::group_by(stage) %>%
  dplyr::arrange(desc(avgPercent)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(manualAnnotLabel = factor(manualAnnotLabel, levels = unique(manualAnnotLabel)))


propClusterG <- propClusterDf %>% 
  ggplot(aes(x = percent, y = manualAnnotLabel, color = infection)) +
  geom_linerange(aes(xmin = 0, xmax = percent), position = position_dodge(0.7), linetype = "dotted", show.legend = FALSE) +
  geom_point(aes(group = infection), position = position_dodge(0.7), alpha = 1, size = 0.5) +
  geom_text(aes(group = infection, label = percentLabel), vjust = 0.5, hjust = -0.15, position = position_dodge(0.7), size = 1.4, show.legend = FALSE) +
  facet_wrap(~ stage, scales = "free", nrow = 1) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.2))) +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1), byrow = TRUE)) +
  theme_classic() +
  scale_color_manual(values = COLORS[["hivStatus"]]) +
  theme(
    plot.margin = margin(l = 15, r = 30, t = -30),
    axis.line = element_line(color = "#000000", linewidth = BASELINEWIDTHSIZE),
    axis.text = element_text(color = "#000000"),
    axis.text.x = element_text(color = "#000000", size = 6),
    axis.text.y = element_text(color = "#000000", size = 6),
    panel.grid = element_blank(),
    panel.grid.major.x = element_line(color = "#eeeeee", linetype = "dotted"),
    strip.text = element_text(color = "#000000", size = 6),
    axis.title.y = element_blank(),
    axis.title.x = element_text(color = "#000000", size = 6, margin = margin(b = 0)),
    legend.text = element_text(color = "#000000", size = 6),
    legend.box.margin = margin(t = -20),
    legend.box.background = element_rect(color = "#ffffff00", fill = "#ffffff00"),
    legend.background = element_rect(color = "#ffffff00", fill = "#ffffff00"),
    legend.title = element_text(color = "#000000", size = 6, margin = margin(b = 0, t = 0)),
    legend.position = "inside",
    legend.spacing.y = unit(1, "mm"),
    legend.position.inside = c(0.95, 0.85), # c(0,0) bottom left, c(1,1) top-right
    legend.key = element_rect(color = "#00000000", fill = "#00000000"),
    strip.background = element_rect(color = "#ffffff00", fill = NA)) +
  labs(x = "% of cluster with HIV+ cells", color = "HIV status", y = "Cluster (# of total cells)")

propClusterOnlyTfhG <- propClusterDf %>% 
  dplyr::filter(manualAnnot %in% c("CD4 GC-Tfh #1", "CD4 GC-Tfh #2", "CD4 PRDM1+ GC-Tfh")) %>%
  dplyr::mutate(manualAnnot = factor(manualAnnot, levels = unique(manualAnnot))) %>%
  ggplot(aes(x = percent, y = manualAnnot, color = infection)) +
  geom_linerange(aes(xmin = 0, xmax = percent), position = position_dodge(0.7), linetype = "dotted", show.legend = FALSE) +
  geom_point(aes(group = infection), position = position_dodge(0.7), alpha = 1, size = 0.5) +
  geom_text(aes(group = infection, label = percentLabel), vjust = 0.5, hjust = -0.15, position = position_dodge(0.7), size = 1.4, show.legend = FALSE) +
  geom_vline(xintercept = -Inf, linewidth = BASELINEWIDTHSIZE, color = "#000000") +
  geom_hline(yintercept = -Inf, linewidth = BASELINEWIDTHSIZE, color = "#000000") +
  facet_wrap(~ stage, scales = "free_x", nrow = 1) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.3))) +
  guides(colour = guide_legend(override.aes = list(size = 1.5, alpha = 1), byrow = TRUE)) +
  theme_classic() +
  scale_color_manual(values = COLORS[["hivStatus"]]) +
  theme(
    panel.spacing = unit(1.2, "lines"),
    axis.line = element_blank(),
    axis.text = element_text(color = "#000000"),
    axis.text.x = element_text(color = "#000000", size = BASEPTFONTSIZE - 2),
    axis.text.y = element_text(color = "#000000", size = BASEPTFONTSIZE - 2),
    panel.grid = element_blank(),
    panel.grid.major.x = element_line(color = "#eeeeee", linetype = "dotted"),
    strip.text = element_text(color = "#000000", size = BASEPTFONTSIZE - 1),
    axis.title.y = element_blank(),
    axis.title.x = element_text(color = "#000000", size = BASEPTFONTSIZE - 1),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(color = "#000000", size = 6),
    legend.margin = margin(0, 0, 0, 0),
    strip.background = element_rect(color = "#ffffff00", fill = NA)) +
  labs(x = "% of cluster with HIV+ cells", y = "Cluster (# of total cells)")



################################################################################
# proportion of cluster infected by donor (for just prdm1+ gc-tfh)
################################################################################
prdmByDonor <- data.frame(
  donor = projUnhashedClean$donorClean2,
  stage = projUnhashedClean$stage,
  infection = projUnhashedClean$hivStatus,
  manualAnnot = projUnhashedClean$manualAnnot
) %>%
  dplyr::filter(stage != "HIV-") %>%
  dplyr::mutate(
    manualAnnot = factor(manualAnnot, levels = uniqClusters),
    stage = factor(stage, levels = uniqStages),
    infection = factor(infection, levels = uniqInfection)
  ) %>%
  dplyr::group_by(stage, donor, manualAnnot, infection, .drop = FALSE) %>%
  dplyr::tally(name = "nCells") %>%
  dplyr::group_by(stage, donor, manualAnnot) %>%
  dplyr::mutate(
    totalCells = sum(nCells),
    percent = nCells / totalCells * 100,
    percent = ifelse(is.nan(percent), 0, percent)) %>%
  dplyr::filter(infection != "HIV-")


prdmByDonorCompiledStats <- ggpubr::compare_means(percent ~ infection,
  prdmByDonor %>% dplyr::filter(manualAnnot == "CD4 PRDM1+ GC-Tfh"),
  group.by = "stage", paired = TRUE)
prdmByDonorCompiledStats$p.sym <- pValSymnum(prdmByDonorCompiledStats$p.adj)

prdmByDonorG <- prdmByDonor %>% 
  dplyr::filter(manualAnnot == "CD4 PRDM1+ GC-Tfh") %>%
  ggplot(aes(y = percent, x = infection)) +
  ggpubr::stat_pvalue_manual(
    prdmByDonorCompiledStats, label = "{p.sym}",
    y.position = c(13.5, 0.35), size = BASEFONTSIZE * 0.6, vjust = 0, tip.length = 0.01) +
  geom_vline(xintercept = -Inf, linewidth = BASELINEWIDTHSIZE, color = "#000000") +
  geom_hline(yintercept = 0, linewidth = BASELINEWIDTHSIZE, color = "#000000") +
  geom_boxplot(outliers = FALSE, aes(color = infection), width = 0.6) +
  geom_jitter(aes(fill = infection), color = "#000000", pch = 21, height = 0, alpha = 0.75, size = 0.8, show.legend = FALSE, width = 0.2) +
  facet_wrap(~ stage, scales = "free_y") +
  coord_cartesian(clip = "off") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  theme_classic() +
  scale_fill_manual(values = COLORS[["hivStatus"]]) +
  scale_color_manual(values = COLORS[["hivStatus"]]) +
  theme(
    plot.margin = margin(r = -30, l = -10),
    legend.position = "blank",
    axis.line = element_blank(),
    axis.text = element_text(color = "#000000"),
    axis.text.x = element_text(color = "#000000", size = BASEPTFONTSIZE - 2, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(color = "#000000", size = BASEPTFONTSIZE - 2),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(color = "#eeeeee", linetype = "dotted"),
    strip.text = element_text(color = "#000000", size = BASEPTFONTSIZE - 1),
    axis.title.y = element_text(color = "#000000", size = BASEPTFONTSIZE - 1),
    axis.title.x = element_blank(),
    strip.background = element_rect(color = "#ffffff00", fill = NA)) +
  labs(y = "% HIV+ of\nPRDM1+ GC-Tfh cells")


################################################################################
# depletion graph
################################################################################
depletionDf <- data.frame(
  donor = projUnhashedClean$donorClean2,
  stage = projUnhashedClean$stage,
  manualAnnot = projUnhashedClean$manualAnnot
) %>%
  dplyr::mutate(
    manualAnnot = factor(manualAnnot, levels = uniqClusters),
    stage = factor(stage, levels = c("HIV-", "Viremic", "ART")),
  ) %>%
  dplyr::group_by(stage, manualAnnot, .drop = FALSE) %>%
  dplyr::tally(name = "nCells") %>%
  dplyr::filter(grepl("Tfh", manualAnnot)) %>%
  dplyr::group_by(stage) %>%
  dplyr::mutate(
    totalCells = sum(nCells),
    percent = nCells / totalCells * 100)

depletionDf2 <- data.frame(
  donor = projUnhashedClean$donorClean2,
  stage = projUnhashedClean$stage,
  manualAnnot = projUnhashedClean$manualAnnot
) %>%
  dplyr::mutate(
    manualAnnot = factor(manualAnnot, levels = uniqClusters),
    stage = factor(stage, levels = c("HIV-", "Viremic", "ART")),
  ) %>%
  dplyr::group_by(stage, donor, manualAnnot, .drop = FALSE) %>%
  dplyr::tally(name = "nCells") %>%
  dplyr::filter(grepl("Tfh", manualAnnot)) %>%
  dplyr::group_by(stage, donor) %>%
  dplyr::mutate(manualAnnot = gsub("CD4 ", "", manualAnnot)) %>%
  dplyr::mutate(
    totalCells = sum(nCells),
    percent = nCells / totalCells * 100,
) %>%
  dplyr::arrange(manualAnnot, stage) %>%
  dplyr::mutate(combinedFactor = paste0(manualAnnot, "_", stage)) %>%
  dplyr::mutate(combinedFactor = factor(combinedFactor, levels = unique(combinedFactor)))


depletionDf2 %>%
  dplyr::mutate(prdm1 = ifelse(grepl("PRDM1", manualAnnot), "PRDM1+", "PRDM1-")) %>%
  dplyr::group_by(stage, donor, prdm1) %>%
  dplyr::summarize(
    sumOfNCells = sum(nCells)
  ) %>%
  tidyr::pivot_wider(id_cols = c("stage", "donor"), names_from = "prdm1", values_from = "sumOfNCells") %>%
  dplyr::mutate(ratio = `PRDM1-` / `PRDM1+`) %>%
  dplyr::group_by(stage) %>%
  dplyr::summarize(meanRatio = mean(ratio))


depletionCompiledStats <- ggpubr::compare_means(percent ~ combinedFactor,
  depletionDf2, paired = FALSE)

depletionCompiledStats$p.sym <- pValSymnum(depletionCompiledStats$p.adj)
depletionCompiledStats <- depletionCompiledStats %>%
  dplyr::filter(p.sym != "ns") %>%
  dplyr::mutate(p.position = 115 - 5 * (row_number() %% 6)) %>%
  dplyr::select(-method, -p.signif, -p.format)

depletionDfG <- depletionDf2 %>%
  ggplot(aes(y = percent, x = combinedFactor)) +
  geom_line(aes(group = donor), color = "#bbbbbb75", linewidth = 0.3) +
  geom_vline(xintercept = -Inf, linewidth = BASELINEWIDTHSIZE, color = "#000000") +
  geom_hline(yintercept = 0, linewidth = BASELINEWIDTHSIZE, color = "#000000") +
  geom_boxplot(outliers = FALSE, aes(color = stage), width = 0.6) +
  geom_point(color = "#000000", pch = 21, alpha = 0.75, size = 0.7, show.legend = FALSE) +
  guides(x = legendry::guide_axis_nested(
    key = "_",
    type = "bracket",
    angle = 45,
    levels_text = list(
      NULL,
      element_text(angle = 0, hjust = 0.5, vjust = 1, color = "#000000", size = 6)))) +
  ggpubr::stat_pvalue_manual(
    depletionCompiledStats,
    label = "{p.sym}", y.position = "p.position",
    size = BASEFONTSIZE * 0.7, vjust = 0.4, tip.length = 0.007, hide.ns = TRUE) +
  scale_color_manual(values = COLORS$stage) +
  coord_cartesian(clip = "off") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  theme_classic() +
  theme(
    plot.margin = margin(r = -1, l = -5),
    legend.position = "blank",
    axis.line = element_blank(),
    axis.text = element_text(color = "#000000"),
    axis.text.x = element_text(color = "#000000", size = BASEPTFONTSIZE - 2),
    axis.text.y = element_text(color = "#000000", size = BASEPTFONTSIZE - 2),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(color = "#eeeeee", linetype = "dotted"),
    strip.text = element_text(color = "#000000", size = BASEPTFONTSIZE - 1),
    axis.title.y = element_text(color = "#000000", size = BASEPTFONTSIZE - 1),
    axis.title.x = element_blank(),
    strip.background = element_rect(color = "#ffffff00", fill = NA)) +
  labs(y = "% of donor's Tfh cells")


################################################################################
# dea analysis
################################################################################
Idents(seuMerged) = "manualAnnot"
tmp <- FindMarkers(
  seuMerged,
  ident.1 = "CD4 GC-Tfh #1",
  ident.2 = "CD4 GC-Tfh #2",
  min.pct = 0.3,
  assay = "ADT",
  features = tsaCatalog[!tsaCatalog$isCtrl, "DNA_ID"])

tmp %>%
  filter(p_val_adj < 0.05) %>% 
  mutate(DNA_ID = rownames(.)) %>%
  mutate(direction = avg_log2FC > 0) %>%
  dplyr::group_by(direction) %>%
  left_join(tsaCatalog[, c("DNA_ID", "cleanName")]) %>%
  dplyr::slice_max(abs(avg_log2FC), n = 15) %>%
  print(n = 30)

tmpRna <- FindMarkers(
  seuMerged,
  ident.1 = "CD4 GC-Tfh #1",
  ident.2 = "CD4 GC-Tfh #2",
  min.pct = 0.3,
  assay = "RNA")

tmpForCs <- tmpRna %>%
  filter(p_val_adj < 0.05) %>% 
  mutate(gene = rownames(.)) %>%
  mutate(direction = avg_log2FC > 0) %>%
  # dplyr::group_by(direction) %>%
  filter(direction) %>%
  dplyr::slice_max(abs(avg_log2FC), n = 100) 

VlnPlot(seuMerged, features = c("STAT4"), group.by = "manualAnnot", pt.size = 0) + NoLegend()


################################################################################
# DEA
################################################################################
dea_blimp1 <- FindMarkers(
  object = seuMerged,
  ident.1 = "CD4 PRDM1+ GC-Tfh",
  ident.2 = c("CD4 GC-Tfh #1", "CD4 GC-Tfh #2"),
  min.pct = 0.3,
  test.use = "LR",
  latent.vars = "donorClean2",
  features = tsaCatalog[!tsaCatalog$isCtrl, "DNA_ID"],
  assay = "ADT")

dea_blimp1 <- dea_blimp1 %>%
  dplyr::mutate(DNA_ID = rownames(.)) %>%
  left_join(dplyr::select(tsaCatalog, DNA_ID, cleanName)) %>%
  dplyr::arrange(desc(avg_log2FC))

deLollipopTheme <- theme(
  strip.clip = "off",
  strip.text = element_text(color = "#000000", size = BASEPTFONTSIZE - 1),
  strip.background = element_rect(color = "#00000000", fill = "#00000000"),
  # axis.text = element_text(color = "#000000"),
  axis.text.x = element_text(color = "#000000", size = BASEPTFONTSIZE - 2),
  axis.text.y = element_text(color = "#000000", size = BASEPTFONTSIZE - 2),
  axis.title.x = element_text(color = "#000000", size = BASEPTFONTSIZE - 1),
  axis.title.y = element_blank(),
  legend.position = "blank",
  panel.grid = element_blank(),
  axis.line = element_line(color = "#000000", linewidth = BASELINEWIDTHSIZE),
  plot.background = element_rect(fill = "#ffffff00", color = "#ffffff00"),
  panel.background = element_rect(fill = "#ffffff00", color = "#ffffff00"))


dea_blimp1G <- dea_blimp1 %>%
  dplyr::mutate(padj2 = p.adjust(p_val, method = "BH")) %>%
  dplyr::filter(padj2 < 0.05) %>%
  dplyr::mutate(direction = ifelse(avg_log2FC > 0, "Up in\nPRDM1+ GC-Tfh", "Up in\nPRDM1- GC-TFh")) %>%
  dplyr::group_by(direction) %>%
  dplyr::slice_max(abs(avg_log2FC), n = 10) %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::mutate(cleanName = gsub(" \\(Neuropilin-1\\)", "", cleanName)) %>%
  dplyr::mutate(cleanName = gsub(" \\(DNAM-1\\)", "", cleanName)) %>%
  dplyr::mutate(cleanName = gsub(" \\(Ecto-5'-nucleotidase\\)", "", cleanName)) %>%
  dplyr::mutate(cleanName = gsub(" \\(Nectin-2\\)", "", cleanName)) %>%
  dplyr::mutate(cleanName = factor(cleanName, levels = unique(cleanName))) %>%
  ggplot(aes(x = avg_log2FC, y = cleanName, color = direction)) +
  geom_linerange(aes(xmin = 0, xmax = avg_log2FC), linetype = "dotted") +
  geom_point(size = 1) +
  geom_vline(xintercept = 0) +
  facet_grid(direction ~ ., scales = "free_y") +
  scale_color_manual(values = c("#0000ff", "#ff0000")) +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.1))) +
  labs(x = "Average log2(FC)") +
  deLollipopTheme +
  theme(
    axis.text.y = element_text(color = "#000000", size = BASEPTFONTSIZE - 2),
    plot.margin = margin(r = -10, l = -25, t = 5, b = 5)
  )


################################################################################
# DEA by subset
################################################################################
dea_blimp1_subsetted <- lapply(unique(seuMerged$stage), function(x) {
  dea <- FindMarkers(
    object = subset(seuMerged, subset = stage == x),
    ident.1 = "CD4 PRDM1+ GC-Tfh",
    ident.2 = c("CD4 GC-Tfh #1", "CD4 GC-Tfh #2"),
    min.pct = 0.3,
    test.use = "LR",
    latent.vars = "donorClean2",
    features = tsaCatalog[!tsaCatalog$isCtrl, "DNA_ID"],
    assay = "ADT")
  
  dea <- dea %>%
    dplyr::mutate(DNA_ID = rownames(.)) %>%
    dplyr::mutate(padj2 = p.adjust(p_val, method = "BH"))
  
  return(dea)
})

names(dea_blimp1_subsetted) <- unique(seuMerged$stage)

dea_blimp1_df <- bind_rows(dea_blimp1_subsetted, .id = "stage") %>%
  left_join(tsaCatalog[, c("DNA_ID", "cleanName")]) %>%
  dplyr::select(-p_val) %>%
  dplyr::arrange(desc(avg_log2FC))  %>%
  dplyr::mutate(direction = ifelse(avg_log2FC > 0, "Up in\nPRDM1+ GC-Tfh", "Up in\nPRDM1- GC-TFh")) %>%
  dplyr::mutate(stage = factor(stage, levels = names(COLORS$stage))) %>%
  dplyr::group_by(stage, direction) %>%
  dplyr::slice_max(abs(avg_log2FC), n = 10) %>%
  dplyr::group_by(stage) %>%
  dplyr::arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  dplyr::mutate(cleanName = gsub(" \\(Neuropilin-1\\)", "", cleanName)) %>%
  dplyr::mutate(cleanName = gsub(" \\(DNAM-1\\)", "", cleanName)) %>%
  dplyr::mutate(cleanName = gsub(" \\(Ecto-5'-nucleotidase\\)", "", cleanName)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(rowCleanName = paste0(row_number(), "_", cleanName)) %>%
  dplyr::mutate(rowCleanName = factor(rowCleanName, levels = rowCleanName))


sup_blimp1_dea_plot <- ggplot(dea_blimp1_df, aes(x = avg_log2FC, y = rowCleanName, color = direction)) +
  geom_linerange(aes(xmin = 0, xmax = avg_log2FC), linetype = "dotted") +
  geom_point(size = 1) +
  geom_vline(xintercept = 0) +
  facet_wrap( ~ stage, scales = "free_y") +
  scale_color_manual(values = c("#0000ff", "#ff0000")) +
  scale_x_continuous(expand = expansion(mult = c(0.2, 0.1))) +
  scale_y_discrete(labels = function(x) gsub("\\d+_", "", x)) +
  labs(x = "Average log2(FC)") +
  deLollipopTheme

################################################################################
# DEG
################################################################################
deg_blimp1 <- FindMarkers(
  object = seuMerged,
  ident.1 = "CD4 PRDM1+ GC-Tfh",
  ident.2 = c("CD4 GC-Tfh #1", "CD4 GC-Tfh #2"),
  test.use = "MAST",
  latent.vars = c("donorClean2", "stage"),
  min.pct = 0.3,
  assay = "RNA")

deg_blimp1G <- deg_blimp1 %>%
  dplyr::mutate(gene = rownames(.)) %>%
  dplyr::mutate(padj2 = p.adjust(p_val, method = "BH")) %>%
  dplyr::filter(padj2 < 0.05) %>%
  dplyr::mutate(direction = ifelse(avg_log2FC > 0, "Up in\nPRDM1+ GC-Tfh", "Up in\nPRDM1- GC-TFh")) %>%
  dplyr::group_by(direction) %>%
  dplyr::slice_max(abs(avg_log2FC), n = 10) %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::mutate(gene = factor(gene, levels = unique(gene))) %>%
  ggplot(aes(x = avg_log2FC, y = gene, color = direction)) +
  geom_linerange(aes(xmin = 0, xmax = avg_log2FC), linetype = "dotted") +
  geom_point(size = 1) +
  geom_vline(xintercept = 0) +
  facet_grid(direction ~ ., scales = "free_y") +
  scale_color_manual(values = c("#0000ff", "#ff0000")) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1))) +
  labs(x = "Average log2(FC)") +
  deLollipopTheme +
  theme(
    axis.text.y = element_text(color = "#000000", size = BASEPTFONTSIZE - 2, face = "italic"),
    plot.margin = margin(r = -50, l = -10, t = 5, b = 5)
  )

################################################################################
# in vitro activated PBMC
################################################################################ 
# gate strategy
gate_flow_png <- image_read("external_data/iv_infection_sort.png")
gate_flow_gg <- image_ggplot(gate_flow_png) 

# representative flow
rep_iv_png <- image_read("external_data/iv_infection.png")
rep_iv_gg <- image_ggplot(rep_iv_png)

# iv data
iv_df <- read.csv("external_data/2025_ivinfection_compiled.csv", check.names = FALSE)

iv_df_2 <- iv_df %>%
  tidyr::pivot_longer(cols = c(-sample), names_to = "gate", values_to = "percent") %>%
  dplyr::mutate(sample = gsub(".fcs", "", sample)) %>%
  dplyr::mutate(sample = gsub("_0_(\\d+)", "_0.\\1", sample)) %>%
  tidyr::separate(sample, into = c("donor", "tp", "condition", "moi"), sep = "_") %>%
  # dplyr::filter(donor != "ND651") %>%
  dplyr::mutate(tp = gsub("hours", "hrs", tp)) %>%
  dplyr::filter(is.na(moi) | moi != 0.01) %>%
  dplyr::mutate(condition = gsub("no", "No", condition)) %>%
  dplyr::mutate(condition = factor(condition, levels = unique(condition))) %>%
  dplyr::mutate(conditionX = paste0(condition, " (", tp, ")")) %>%
  dplyr::mutate(conditionX = factor(conditionX, levels = unique(conditionX))) %>%
  dplyr::mutate(parent = ifelse(grepl("p24", gate), "p24+", "non-naive CD8-")) %>%
  dplyr::mutate(gate2 = gsub("p24\\+ ", "", gate)) %>%
  dplyr::mutate(gate2 = stringr::str_replace(gate2, "(BLIMP1)(\\+|-) (CCR5)(\\+|-)", "\\1<sup>\\2</sup> \\3<sup>\\4</sup>")) %>%
  dplyr::arrange(donor, tp, gate, conditionX) %>%
  dplyr::mutate(gate2_condition = paste0(condition, "_", gate2)) %>%
  dplyr::mutate(gate2_condition = factor(gate2_condition, levels = unique(gate2_condition)))

iv_conditions_pal <- RColorBrewer::brewer.pal(length(levels(iv_df_2$condition)), "Dark2")
names(iv_conditions_pal) <- levels(iv_df_2$condition)

# nn_stats <- ggpubr::compare_means(percent ~ condition,
#   iv_df_2 %>% dplyr::filter(parent == "non-naive CD8+"),
#   group.by = c("gate", "tp"), p.adjust.method = "BH", paired = TRUE, method = "wilcox")
# nn_stats$p.sym <- pValSymnum(nn_stats$p.adj)
# print(nn_stats, n = 40)

iv_graph_settings <- list(
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))),
  # scale_fill_manual(values = iv_conditions_pal),
  scale_color_manual(values = iv_conditions_pal),
  theme_classic(),
  coord_cartesian(clip = "off"),
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)),
  theme(
    axis.line = element_blank(),
    axis.text.y = element_markdown(size = BASEPTAXISFONTSIZE, color = "#000000"),
    axis.title.y = element_text(size = BASEPTFONTSIZE - 2, color = "#000000"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_markdown(size = BASEPTFONTSIZE - 2.75, color = "#000000"),
    legend.position = "bottom",
    legend.key.size = unit(5, "pt"),
    legend.box.margin = margin(r = 20),
    panel.grid.major.y = element_line(color = "#eeeeee", linetype = "dotted"),
    legend.text = element_text(size = BASEPTFONTSIZE - 2, color = "#000000"),
    legend.title = element_blank(),
    strip.clip = "off",
    axis.title.x = element_blank())
)

nn_iv_g <- iv_df_2 %>%
  dplyr::filter(parent == "non-naive CD8-") %>%
  ggplot(aes(x = condition, y = percent)) +
  geom_hline(yintercept = -Inf, linewidth = BASELINEWIDTHSIZE) +
  geom_boxplot(width = 0.6, linewidth = 0.5, aes(color = condition), outliers = FALSE) +
  # stat_summary(geom = "errorbar", fun.data = mean_se) +
  # stat_summary(fun = mean, geom = "bar") +
  geom_point(aes(group = condition), color = "#000000", pch = 21, alpha = 0.75, size = 0.5, show.legend = FALSE) +
  ggh4x::facet_nested(~ tp + gate2,
    scales = "free_x",
    solo_line = TRUE,
    axes = FALSE,
    nest_line = element_line(colour = "#000000"),
    strip = ggh4x::strip_nested(clip = "off")) +
  labs(color = "Condition", y = "% of non-naive\nCD8- CD3+ cells") +
  iv_graph_settings

p24_df <- iv_df_2 %>%
  dplyr::filter(parent == "p24+")

# p24_stats <- ggpubr::compare_means(percent ~ condition,
#   p24_df,
#   group.by = c("gate", "tp"), p.adjust.method = "BH", paired = TRUE)
# p24_stats$p.sym <- pValSymnum(p24_stats$p.adj)
# p24_stats <- p24_stats %>%
#   dplyr::filter(p.sym != "ns") %>%
#   dplyr::mutate(p.position = 100 - 5 * (row_number() %% 5)) %>%
#   dplyr::select(-method, -p.signif, -p.format)

p24_iv_g <- p24_df %>%
  ggplot(aes(x = condition, y = percent)) +
  geom_hline(yintercept = -Inf, linewidth = BASELINEWIDTHSIZE) +
  geom_boxplot(width = 0.6, linewidth = 0.5, aes(color = condition), outliers = FALSE) +
  # stat_summary(geom = "errorbar", fun.data = mean_se) +
  # stat_summary(fun = mean, geom = "bar", aes(fill = condition)) +
  geom_point(aes(group = condition), color = "#000000", pch = 21, alpha = 0.75, size = 0.5, show.legend = FALSE) +
  ggh4x::facet_nested(~ tp + gate2, scales = "free_x", solo_line = TRUE, 
    axes = FALSE, nest_line = element_line(colour = "#000000"), strip = ggh4x::strip_nested(clip = "off")) +
  labs(color = "Condition", y = "% of non-naive\nCD8- CD3+ p24+ cells") +
  iv_graph_settings


################################################################################
# organoid infection
################################################################################
# representative flow
rep_flow_png <- image_read("external_data/organoid_rep.png")
rep_flow_gg <- image_ggplot(rep_flow_png)

# organoid methodology
organoid_png <- image_read("external_data/tonsil_organoid_workflow.png")
organoid_gg <- image_ggplot(organoid_png)

# organoid gating
organoid_gate_png <- image_read("external_data/organoid_gate_strategy.png")
organoid_gate_gg <- image_ggplot(organoid_gate_png)

# organoid data
organoid_df <- read.csv("external_data/2025_tonsil_compiled.csv", check.names = FALSE)
organoid_df <- organoid_df[, colnames(organoid_df) != ""]

organoid_df2 <- organoid_df %>%
  tidyr::pivot_longer(cols = c(-"donorID", -"day", -"condition"), names_to = "gate", values_to = "percent") %>%
  dplyr::mutate(gate2 = stringr::str_replace(gate, "(BLIMP1) (dim|low|high)  (p24)(\\+|-)", "\\1<sup>\\2</sup> \\3<sup>\\4</sup>")) %>%
  dplyr::filter(!grepl("0.01", condition)) %>%
  dplyr::mutate(blimp1Category = str_match(gate, pattern = "(dim|low|high)")[, 1]) %>%
  dplyr::mutate(blimp1Category = factor(blimp1Category, levels = c("low", "dim", "high"))) %>%
  dplyr::mutate(p24Category = grepl("\\+", gate2)) %>%
  dplyr::arrange(p24Category, blimp1Category) %>%
  dplyr::mutate(gate2 = factor(gate2, levels = unique(gate2))) %>%
  dplyr::mutate(condition = factor(condition, levels = c("Neg", "NL4.3 MOI=0.01", "NL4.3 MOI=0.1", "SUMA MOI=0.01", "SUMA MOI=0.1")))
  

org_stats <- ggpubr::compare_means(
  percent ~ condition,
  organoid_df2,
  group.by = "gate2", p.adjust.method = "BH", paired = TRUE)
org_stats

organoid_g <- organoid_df2 %>%
  ggplot(aes(x = gate2, y = percent)) +
  geom_hline(yintercept = -Inf, linewidth = BASELINEWIDTHSIZE) +
  geom_vline(xintercept = -Inf, linewidth = BASELINEWIDTHSIZE) +
  geom_boxplot(aes(color = condition), width = 0.5, linewidth = 0.5, outliers = FALSE, position = position_dodge(width = 0.95)) +
  geom_point(aes(group = condition), color = "#000000", pch = 21, alpha = 0.75, size = 0.5,position = position_dodge(width = 0.95), show.legend = FALSE) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  guides(x = legendry::guide_axis_nested(
    key = " ",
    type = "bracket")) +
  theme_classic() +
  scale_color_manual(values = c("#3d3d3d", "#bdd7e7", "#3182bd", "#fcae91", "#de2d26")) +
  labs(color = "Condition", y = "% of non-naive\nCD8- CD3+ cells") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  coord_cartesian(clip = "off") +
  theme(
    axis.line = element_blank(),
    axis.text.y = element_markdown(size = BASEPTFONTSIZE - 2, color = "#000000"),
    axis.title.y = element_text(size = BASEPTFONTSIZE - 2, color = "#000000"),
    axis.text.x = element_markdown(size = BASEPTFONTSIZE - 2, color = "#000000"),
    legend.position = "bottom",
    legend.key.size = unit(5, "pt"),
    panel.grid.major.y = element_line(color = "#eeeeee", linetype = "dotted"),
    legend.text = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    legend.title = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    axis.title.x = element_blank())

organoid_p24_df <- organoid_df2 %>%
  dplyr::filter(p24Category) %>%
  dplyr::group_by(donorID, gate2) %>%
  dplyr::mutate(percent_bg_subtract = percent - percent[1]) %>%
  dplyr::group_by(donorID, condition) %>%
  dplyr::mutate(pp2 = percent_bg_subtract / sum(percent_bg_subtract) * 100) %>%
  dplyr::filter(condition != "Neg") %>%
  dplyr::mutate(cleanX = gsub(" p24\\S+", "", gate2)) %>%
  dplyr::mutate(cleanX = factor(cleanX, levels = unique(cleanX)))

org_p24_stats <- ggpubr::compare_means(
  pp2 ~ condition,
  organoid_p24_df,
  group.by = "cleanX", p.adjust.method = "BH", paired = TRUE)
org_p24_stats$p.sym <- pValSymnum(org_p24_stats$p.adj)

organoid_g2 <- organoid_p24_df %>%
  ggplot(aes(x = condition, y = pp2)) +
  geom_hline(yintercept = -Inf, linewidth = BASELINEWIDTHSIZE) +
  geom_boxplot(aes(color = condition), width = 0.5, linewidth = 0.5, outliers = FALSE, position = position_dodge(width = 0.95)) +
  geom_point(aes(group = condition), color = "#000000", pch = 21, alpha = 0.75, size = 0.5,position = position_dodge(width = 0.95), show.legend = FALSE) +
  ggpubr::stat_pvalue_manual(
    org_p24_stats, label = "{p.sym}",
    y.position = 75, size = BASEFONTSIZE * 0.9, vjust = 0.4, tip.length = 0.01) +
  facet_wrap(~ cleanX, scales = "free_x") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_color_manual(values = c("#bdd7e7", "#3182bd", "#fcae91", "#de2d26")) +
  labs(color = "Condition", y = "Background-subtracted % of \nnon-naive CD8- CD3+ p24+ cells") +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(
    axis.line.x = element_blank(),
    axis.line.y = element_line(linewidth = BASELINEWIDTHSIZE, color = "#000000"),
    axis.text.y = element_markdown(size = BASEPTFONTSIZE - 2, color = "#000000"),
    axis.title.y = element_text(size = BASEPTFONTSIZE - 2, color = "#000000"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "blank",
    legend.key.size = unit(5, "pt"),
    strip.background = element_rect(color = "#ffffff00", fill = "#ffffff00"),
    strip.text = element_markdown(size = BASEPTFONTSIZE - 2, color = "#000000"),
    panel.grid.major.y = element_line(color = "#eeeeee", linetype = "dotted"),
    legend.text = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    legend.title = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    axis.title.x = element_blank())
# 
# organoid_df %>%
#   dplyr::group_by(condition, p24) %>%
#   dplyr::reframe(pp = percent / sum(percent) * 100) %>%
#   print(n = 200)

  
################################################################################
# final plots
################################################################################
layout <- c(
  area(1, 1, 5, 5),
  area(1, 6, 3, 9),
  area(4, 6, 5, 9),
  area(6, 1, 7, 2), # prdm1
  area(6, 3, 7, 5), # depletion
  area(6, 6, 7, 7), # prdm1 dea
  area(6, 8, 7, 9), # prdm1 deg
  area(8, 1, 10, 3), # representative plot
  area(8, 4, 10, 7),
  area(8, 8, 10, 9)
)

p <- atacFigA + theme(plot.margin = margin(t = 0, b = -90)) +
  patchwork::free(rnaFigB_mod, type = "panel", side = "lr") +
  patchwork::free(propClusterOnlyTfhG +
      theme(plot.margin = margin(t = -15, l = 0, r = 15, b = -20),
        legend.box.margin = margin(b = 0, t = -20)), type = "panel", side = "blr") +
  patchwork::free(prdmByDonorG + theme(plot.margin = margin(t = 0)), type = "panel", side = "lr") +
  patchwork::free(depletionDfG + theme(plot.margin = margin(t = 0, r = 5)), type = "panel", side = "lr") +
  patchwork::free(dea_blimp1G, type = "panel", side = "tb") +
  patchwork::free(deg_blimp1G, type = "panel", side = "tb") +
  wrap_elements(full = rep_flow_gg + theme(plot.margin = margin(l = 0, r = 0, t = 10))) +
  free(organoid_g + 
      theme(
        plot.margin = margin(l = 0, r = -150),
        legend.box.margin = margin(l = 105),
        plot.tag = element_text(margin = margin(l = 5))), type = "panel", side = "tr") +
  free(organoid_g2 + 
      theme(
        plot.margin = margin(l = -20, r = -25),
        plot.tag = element_text(margin = margin(l = 5))), type = "panel", side = "tlr") +
  plot_annotation(tag_levels = "A") +
  plot_layout(design = layout) &
  plotTagTheme

saveFinalFigure(
  plot = p,
  fn = "fig5_v4",
  prefixDir = "figures/outs",
  units = "cm",
  devices = c("png", "pdf"),
  gwidth = JCI_MAX_WIDTH,
  gheight = 8 * 2.54,
  addTimestamp = TRUE)


layoutP <- c(
  area(1, 1, 8, 9),
  area(9, 1, 10, 9)
)

supP <- patchwork::free(propClusterG + theme(plot.margin = margin(l = 5, r = 25)), type = "panel", side = "lr") +
  patchwork::free(sup_blimp1_dea_plot, type = "panel", side = "lr") +
  plot_annotation(tag_levels = "A") +
  plot_layout(design = layoutP) &
  plotTagTheme

saveFinalFigure(plot = supP,
  fn = "supfig6_v4",
  prefixDir = "figures/outs", devices = c("png", "pdf"),
  gwidth = 8, gheight = 8,
  addTimestamp = TRUE)


layoutP2 <- c(
  area(1, 1, 3, 9),
  area(4, 1, 6, 4),
  area(4, 5, 6, 9),
  area(7, 1, 9, 3),
  area(7, 4, 9, 9)
)

supP2 <- wrap_elements(full = gate_flow_gg + theme(plot.margin = margin(l = 0, r = 0))) +
  wrap_elements(full = rep_iv_gg + theme(plot.margin = margin(l = 0, r = 0))) +
  (nn_iv_g / p24_iv_g + plot_layout(guides = "collect") & theme(legend.position = "bottom")) +
  wrap_elements(full = organoid_gg + theme(plot.margin = margin(l = 0, r = 0))) +
  wrap_elements(full = organoid_gate_gg + theme(plot.margin = margin(l = 10, r = 0))) +
  plot_annotation(tag_levels = "A") +
  plot_layout(design = layoutP2) &
  plotTagTheme

saveFinalFigure(plot = supP2,
  fn = "supfig7_v4",
  prefixDir = "figures/outs", devices = c("png", "pdf"),
  gwidth = 8, gheight = 8,
  addTimestamp = TRUE)
