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

set.seed(42)

load("rdata/chronicSub2_deviations_2025-02-26-01-56.Rdata")


source("figures/defaultFigureSettings.R")
source("scripts/graphs.R")
source("scripts/colorPalettes.R")

################################################################################
# umaps
################################################################################
projUnhashedClean$stage <- ifelse(projUnhashedClean$stage == "Chronic", "Viremic", projUnhashedClean$stage)

# fix naming for CD4 Trm #1
projUnhashedClean$manualAnnot <- ifelse(projUnhashedClean$manualAnnot == "CD4 Trm #1", "CD4 Trm", projUnhashedClean$manualAnnot)

umapDf <- getEmbedding(projUnhashedClean, embedding = "Harmony_UMAP_Combined_v2_05Dist_50Neigh")

stopifnot(sum(rownames(umapDf) == projUnhashedClean$cellNames) == length(projUnhashedClean$cellNames))

umapDf$stage <- projUnhashedClean$stage
umapDf$manualAnnot <- trimws(projUnhashedClean$manualAnnot)
umapDf$donor <- projUnhashedClean$donorClean2
umapDf$hivStatus <- projUnhashedClean$hivStatus

# umapDf <- slice(umapDf, sample(1:n()))
uniqAnnots <- sort(unique(umapDf$manualAnnot))
annotColorScale <- nColorPalette(n = length(uniqAnnots))
names(annotColorScale) <- uniqAnnots

umapDf$manualAnnot <- factor(umapDf$manualAnnot, levels = uniqAnnots)
umapDf$stage <- factor(umapDf$stage, levels = names(COLORS$stage))

umapFig2Theme <- theme(
  aspect.ratio = 1,
  legend.position = "bottom",
  legend.text = element_text(size = BASEPTFONTSIZE),
  legend.title = element_blank(),
  axis.title = element_text(size = BASEPTFONTSIZE),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  legend.spacing.x = unit(BASEPTFONTSIZE / 2, 'points'),
  plot.title = element_text(size = BASEPTFONTSIZE, hjust = 0.5),
  strip.background = element_rect(color = "#ffffff00", fill = "#ffffff00"),
  strip.text = element_text(size = BASEPTFONTSIZE),
  plot.title.position = "panel",
  axis.line = element_line(linewidth = BASELINEWIDTHSIZE),
  legend.box.margin = margin(t = -10, b = -10),
  legend.key.size = unit(BASEPTFONTSIZE * 1.1, 'points'),
  legend.background = element_rect(fill = "transparent", colour = NA))

figA <- umapDf %>% 
  mutate(hivStatus = factor(hivStatus, levels = names(COLORS$hivStatus))) %>%
  arrange(hivStatus) %>%
  ggplot(aes(x = `HarmonyLsiCleanForSeurat#UMAP_Dimension_1`, y = `HarmonyLsiCleanForSeurat#UMAP_Dimension_2`, color = hivStatus)) +
  rasterize(geom_point(alpha = 0.6, size = 0.1), dpi = 300) +
  guides(colour = guide_legend(override.aes = list(size = 4), nrow = 1)) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  facet_wrap(~ stage, nrow = 1) +
  scale_color_manual(values = COLORS$hivStatusForUMAP) +
  theme_classic() +
  umapFig2Theme


uniqAnnots <- sort(unique(umapDf$manualAnnot))
annotColorScale <- nColorPalette(n = length(uniqAnnots))
names(annotColorScale) <- uniqAnnots
umapDf$manualAnnot <- factor(umapDf$manualAnnot, levels = uniqAnnots)
umapDf$hivPos <- umapDf$hivStatus != "HIV-"


## plot for david
# p1 <- umapDf %>% 
#   mutate(tmpStage = "All (by cell type)") %>%
#   mutate(hivStatus = factor(hivStatus, levels = names(COLORS$hivStatus))) %>%
#   ggplot(aes(x = `HarmonyLsiCleanForSeurat#UMAP_Dimension_1`, y = `HarmonyLsiCleanForSeurat#UMAP_Dimension_2`, color = manualAnnot)) +
#   rasterize(geom_point(alpha = 0.6, size = 0.1, show.legend = FALSE), dpi = 300) +
#   scale_color_manual(values = annotColorScale) +
#   facet_wrap(~ tmpStage, nrow = 1) +
#   labs(x = "UMAP 1", y = "UMAP 2") +
#   theme_classic() +
#   umapFig2Theme
  
# p2df <- umapDf %>% 
#   # mutate(hivStatus = factor(hivStatus, levels = names(COLORS$hivStatus))) %>%
#   arrange(hivPos) %>%
#   dplyr::group_by(stage) %>%
#   dplyr::mutate(nHivCells = sum(hivPos)) %>%
#   dplyr::mutate(nStage = n()) %>%
#   dplyr::mutate(stage2 = glue("{stage} donors\n{nHivCells} HIV+ cells ({round(nHivCells / nStage * 100, 2)}%)"))
# 
# tmp <- p2df %>%
#   group_by(stage, stage2) %>%
#   summarize(totalCells = n())
# 
# new_labels <- tmp$stage2
# names(new_labels) <- tmp$stage
# 
# p2 <- ggplot(p2df, aes(x = `HarmonyLsiCleanForSeurat#UMAP_Dimension_1`, y = `HarmonyLsiCleanForSeurat#UMAP_Dimension_2`, color = as.character(hivPos))) +
#   rasterize(geom_point(alpha = 0.6, size = 0.1, show.legend = FALSE), dpi = 300) +
#   # guides(colour = guide_legend(override.aes = list(size = 4), nrow = 1)) +
#   labs(x = "UMAP 1", y = "UMAP 2") +
#   facet_wrap(~ stage, nrow = 1, labeller = as_labeller(new_labels)) +
#   scale_color_manual(values = c("TRUE" = "#ff0000", "FALSE" = "#cccccc80")) +
#   theme_classic() +
#   umapFig2Theme
# 
# saveFinalFigure(p2, fn = "grant_fig", prefixDir = "figures/outs", devices = "png", gheight = 4, gwidth = 4)
# saveFinalFigure(p1 + p2 + plot_layout(widths = c(1, 3)), fn = "grant_fig", prefixDir = "figures/outs", devices = "png", gheight = 4, gwidth = 6)

################################################################################
# detection upset
################################################################################
hivStats <- data.frame(
  cellName = projUnhashedClean$cellNames,
  donor = projUnhashedClean$donorClean2,
  stage = projUnhashedClean$stage,
  hivDNA = projUnhashedClean$haystackOut,
  hivRNA = projUnhashedClean$hivRNA
)

hivStats$stage <- factor(hivStats$stage, levels = c("Viremic", "ART"))

# aggregate
statGroupingDict <- hivStats$stage
names(statGroupingDict) <- hivStats$cellName

figB <- plotUpsetHivStats(
  rnaCbc = hivStats[hivStats$hivRNA, "cellName"], 
  dnaCbc = hivStats[hivStats$hivDNA, "cellName"],
  groupingVector = statGroupingDict,
  barTextSize = BASEFONTSIZE - 0.4,
  returnPlot = TRUE,
  gheight = 4) + plot_layout(heights = c(6, 1))

figB[[1]] <- figB[[1]] +
  scale_fill_manual(values = COLORS$stage) +
  theme(
    strip.text = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    axis.text = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    axis.title.y = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    axis.line = element_line(linewidth = BASELINEWIDTHSIZE),
    plot.background = element_rect(fill = "#ffffff00", color = "#ffffff00"),
    panel.background = element_rect(fill = "#ffffff00", color = "#ffffff00"),
    plot.margin = margin(t = 10, l = 5))

figB[[2]] <- figB[[2]] +
  theme(
    axis.text = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    axis.title.y = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    plot.background = element_rect(fill = "#ffffff00", color = "#ffffff00"),
    panel.background = element_rect(fill = "#ffffff00", color = "#ffffff00"),
    plot.margin = margin(t = 10, b = 0, l = 5))




statGroupingDict <- hivStats$donor
names(statGroupingDict) <- hivStats$cellName

topTheme <- list(
  scale_fill_manual(values = rep("#000000", 15)),
  theme(
    strip.text = element_text(size = BASEPTFONTSIZE, color = "#000000", hjust = 0.5),
    axis.text = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    axis.title.y = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    axis.line = element_line(linewidth = BASELINEWIDTHSIZE),
    axis.text.y = element_text(size = BASEPTFONTSIZE - 2, color = "#000000"),
    plot.background = element_rect(fill = "#ffffff00", color = "#ffffff00"),
    panel.background = element_rect(fill = "#ffffff00", color = "#ffffff00")))

bottomTheme <- theme(
  axis.text = element_text(size = BASEPTFONTSIZE, color = "#000000"),
  axis.title.y = element_text(size = BASEPTFONTSIZE, color = "#000000"),
  plot.background = element_rect(fill = "#ffffff00", color = "#ffffff00"),
  panel.background = element_rect(fill = "#ffffff00", color = "#ffffff00"))


supFigUpsetArt <- plotUpsetHivStats(
  rnaCbc = hivStats[hivStats$hivRNA & hivStats$stage == "ART", "cellName"], 
  dnaCbc = hivStats[hivStats$hivDNA & hivStats$stage == "ART", "cellName"],
  groupingVector = statGroupingDict,
  returnPlot = TRUE,
  barTextSize = BASEFONTSIZE,
  stripPlacement = "right",
  gheight = 4)

supFigUpsetArt[[1]] <- supFigUpsetArt[[1]] + topTheme
supFigUpsetArt[[2]] <- supFigUpsetArt[[2]] + bottomTheme

supFigUpsetVir <- plotUpsetHivStats(
  rnaCbc = hivStats[hivStats$hivRNA & hivStats$stage == "Viremic", "cellName"], 
  dnaCbc = hivStats[hivStats$hivDNA & hivStats$stage == "Viremic", "cellName"],
  groupingVector = statGroupingDict,
  returnPlot = TRUE,
  barTextSize = BASEFONTSIZE,
  stripPlacement = "right",
  gheight = 4)

supFigUpsetVir[[1]] <- supFigUpsetVir[[1]] + topTheme
supFigUpsetVir[[2]] <- supFigUpsetVir[[2]] + bottomTheme




################################################################################
# lollipops for hiv dna and rna
################################################################################
# fig C - hiv dna
# lollipop for HIV DNA+
figCData <- data.frame(
  cluster = projUnhashedClean$manualAnnot,
  stage = projUnhashedClean$stage,
  haystackOut = projUnhashedClean$hivStatus
) %>%
  dplyr::mutate(cluster = factor(cluster, levels = unique(cluster))) %>%
  dplyr::group_by(stage, cluster, .drop = FALSE) %>%
  dplyr::filter(haystackOut == "HIV DNA+ RNA-") %>%
  dplyr::mutate(stage = factor(stage, levels = c("Viremic", "ART"))) %>%
  dplyr::tally(name = "nCells") %>%
  dplyr::group_by(stage) %>%
  dplyr::mutate(percentStage = nCells / sum(nCells) * 100) %>%
  dplyr::mutate(percentStageLabel = glue("{round(percentStage, digits = 2)}% ({nCells})")) %>%
  dplyr::arrange(desc(percentStage), .by_group = TRUE) %>%
  dplyr::mutate(cluster = factor(cluster, levels = unique(cluster)))
  
figCSums <- figCData %>%
  dplyr::summarize(cellsInStage = sum(nCells)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(stageLabel = glue("{stage} (n={cellsInStage})"))

figCSumDict <- figCSums$stageLabel
names(figCSumDict) <- figCSums$stage

fig2LollipopTheme <- list(
  geom_vline(xintercept = 0, linewidth = BASELINEWIDTHSIZE, color = "#000000"),
  geom_hline(yintercept = -Inf, linewidth = BASELINEWIDTHSIZE, color = "#000000"),
    theme(
      legend.position = "blank",
      axis.line = element_blank(),
      axis.text = element_text(color = "#000000", size = BASEPTFONTSIZE),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(color = "#000000", size = BASEPTFONTSIZE),
      panel.grid = element_blank(),
      axis.text.x = element_text(color = "#000000", size = BASEPTFONTSIZE),
      strip.text = element_text(color = "#000000", size = BASEPTFONTSIZE),
      panel.spacing = unit(10, "pt"),
      strip.background = element_rect(color = "#ffffff00", fill = NA)))
  

figDnaLol <- ggplot(figCData, aes(x = percentStage, y = cluster, color = stage)) +
  theme_classic() +
  fig2LollipopTheme +
  geom_linerange(aes(xmin = 0, xmax = percentStage), linetype = "dotted") +
  geom_point(size = 0.75) +
  geom_text(aes(label = percentStageLabel, x = percentStage), hjust = -0.1, size = 6 / ggplot2:::.pt) +
  facet_wrap(~ stage, labeller = as_labeller(figCSumDict)) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(limits = c(0, 25), expand = expansion(mult = c(0, 0))) +
  scale_color_manual(values = COLORS$stage) +
  labs(x = "% of HIV DNA+ RNA- cells per stage")

# fig d RNA lollipop
figDData <- data.frame(
  cluster = projUnhashedClean$manualAnnot,
  stage = projUnhashedClean$stage,
  haystackOut = projUnhashedClean$hivStatus
) %>%
  dplyr::mutate(cluster = factor(cluster, levels = levels(figCData$cluster))) %>%
  dplyr::filter(haystackOut == "HIV DNA+/- RNA+") %>%
  dplyr::mutate(stage = factor(stage, levels = c("Viremic", "ART"))) %>%
  dplyr::group_by(stage, cluster, .drop = FALSE) %>%
  dplyr::tally(name = "nCells") %>%
  dplyr::group_by(stage) %>%
  dplyr::mutate(percentStage = nCells / sum(nCells) * 100) %>%
  dplyr::mutate(percentStageLabel = glue("{round(percentStage, digits = 2)}% ({nCells})")) %>%
  dplyr::arrange(desc(percentStage), .by_group = TRUE)

figDSums <- figDData %>%
  dplyr::summarize(cellsInStage = sum(nCells)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(stageLabel = glue("{stage} (n={cellsInStage})"))

figDSumDict <- figDSums$stageLabel
names(figDSumDict) <- figDSums$stage

figRnaLol <- ggplot(figDData, aes(x = percentStage, y = cluster, color = stage)) +
  theme_classic() +
  fig2LollipopTheme +
  geom_linerange(aes(xmin = 0, xmax = percentStage), linetype = "dotted") +
  geom_point(size = 0.75) +
  geom_text(aes(label = percentStageLabel, x = percentStage), hjust = -0.1, size = 6 / ggplot2:::.pt) +
  facet_wrap(~ stage, labeller = as_labeller(figDSumDict)) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(limits = c(0, 25), expand = expansion(mult = c(0, 0))) +
  scale_color_manual(values = COLORS$stage) +
  labs(x = "% of HIV DNA+/- RNA+ cells per stage")


################################################################################
# infection of tfh by stage
################################################################################
tmp <- data.frame(
  donor = projUnhashedClean$donorClean2,
  stage = projUnhashedClean$stage,
  infection = projUnhashedClean$hivStatus,
  manualAnnot = projUnhashedClean$manualAnnot
) %>%
  dplyr::filter(stage != "HIV-") %>%
  dplyr::group_by(donor, stage, infection, manualAnnot) %>%
  dplyr::tally(name = "totalCells") %>%
  dplyr::mutate(isTfh = grepl("Tfh", manualAnnot)) %>%
  dplyr::mutate(isTfh = ifelse(isTfh, "Tfh", "Non-Tfh")) %>%
  dplyr::mutate(isTfh = factor(isTfh, levels = c("Non-Tfh", "Tfh"))) %>%
  dplyr::mutate(infection = factor(infection, levels = names(COLORS$hivStatus))) %>%
  dplyr::group_by(donor, stage, infection, isTfh, .drop = FALSE) %>%
  dplyr::summarize(nCells = sum(totalCells)) %>%
  dplyr::filter(infection != "HIV-") %>%
  dplyr::ungroup() %>%
  dplyr::mutate(infection2 = ifelse(infection == "HIV DNA+ RNA-", "DNA+", "RNA+")) %>%
  dplyr::mutate(stage = factor(stage, levels = c("Viremic", "ART"))) %>%
  dplyr::group_by(donor, stage) %>%
  dplyr::mutate(totalCellsInGroup = sum(nCells)) %>%
  dplyr::filter(totalCellsInGroup > 10) %>% # filtering for donors with at least 10 HIV+ cells
  dplyr::mutate(propOfHIVPos = nCells / sum(nCells)) %>%
  dplyr::mutate(propOfHIVPos = ifelse(is.nan(propOfHIVPos), 0, propOfHIVPos)) %>%
  dplyr::mutate(xlbl = paste0(isTfh, "\n", infection2)) %>%
  dplyr::mutate(xlbl = factor(xlbl, levels = c("Non-Tfh\nDNA+", "Tfh\nDNA+", "Non-Tfh\nRNA+", "Tfh\nRNA+")))

figTfhCompStats <- ggpubr::compare_means(
  propOfHIVPos ~ xlbl,
  group.by = c("infection", "stage"),
  p.adjust.method = "holm",
  method = "wilcox.test",
  data = tmp)
figTfhCompStats$p.sym <- pValSymnum(figTfhCompStats$p.adj)


figTfh <- ggplot(tmp, aes(x = xlbl, y = propOfHIVPos * 100)) +
  geom_boxplot(aes(color = stage), outliers = FALSE, show.legend = FALSE) +
  geom_point(alpha = 0.8, pch = 21, color = "#000000", fill = "#ffffff", size = 0.8) +
  geom_line(aes(group = donor), linewidth = 0.5, color = "#cfcfcf50") +
  lemon::facet_rep_wrap(~ stage, nrow = 2) +
  theme_bw() +
  ggpubr::stat_pvalue_manual(figTfhCompStats, label = "{p.sym}",
    y.position = 75, bracket.shorten = 0.4, size = BASEFONTSIZE * 0.85, vjust = 0) +
  scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0.05, 0.05))) +
  scale_color_manual(values = COLORS$stage) +
  theme_classic() +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(t = -180),
    panel.spacing = unit(0, "lines"),
    panel.spacing.y = unit(-0.1, "lines"),
    axis.line = element_line(linewidth = BASELINEWIDTHSIZE),
    axis.text.y = element_text(color = "#000000", size = BASEPTFONTSIZE),
    axis.text.x = element_text(color = "#000000", size = BASEPTFONTSIZE - 1),
    panel.grid = element_blank(),
    panel.background = element_rect(color = "#00000000", fill = "#00000000"),
    plot.background = element_rect(color = "#00000000", fill = "#00000000"),
    strip.text = element_text(color = "#000000", size = BASEPTFONTSIZE),
    axis.title.y = element_text(color = "#000000", size = BASEPTFONTSIZE),
    axis.title.x = element_blank(),
    strip.background = element_rect(color = "#ffffff00", fill = "#ffffff00")) +
  labs(
    x = "Cell phenotype",
    y = "% of donor's HIV+ cells")


library(magick)
flow_sort_strat_png <- image_read("external_data/tfh_flow_same.png")
flow_sort_strat_gg <- image_ggplot(flow_sort_strat_png)

################################################################################
# tfh correlation
################################################################################
ipdaFlow <- read.csv("external_data/ipda_flow_combined_2024-11-22.csv", check.names = FALSE) %>%
  mutate(`3' Defective` = as.numeric(gsub("<", "", `3' Defective`))) %>%
  mutate(`5' Defective` = as.numeric(gsub("<", "", `5' Defective`)))

dogmaCondensed <- data.frame(
  donor = projUnhashedClean$donorClean2,
  hivStatus = projUnhashedClean$hivStatus,
  manualAnnot = projUnhashedClean$manualAnnot,
  stage = projUnhashedClean$stage
) %>%
  dplyr::mutate(stage = factor(stage, levels = c("HIV-", "Viremic", "ART"))) %>%
  dplyr::group_by(donor, stage) %>%
  dplyr::summarize(
    totalCells = n(),
    `HIV DNA+/- RNA+` = sum(hivStatus == "HIV DNA+ RNA+"),
    `HIV-` = sum(hivStatus == "HIV-"),
    `HIV DNA+ RNA-` = sum(hivStatus == "HIV DNA+ RNA-"),
    `HIV+ in Tfh` = sum(grepl("Tfh", manualAnnot) & hivStatus != "HIV-"),
    `Tfh` = sum(grepl("Tfh", manualAnnot))
  )

figTfhCorr <- dogmaCondensed %>%
  left_join(ipdaFlow, by = c("donor" = "Sample")) %>%
  mutate(percentHIVPos = (totalCells - `HIV-`) / totalCells * 100) %>%
  mutate(percentHIVTfhPos = `HIV+ in Tfh` / (totalCells - `HIV-`) * 100) %>%
  mutate(percentTfh = `Tfh` / totalCells * 100) %>%
  ggplot(aes(x = `GC-TFH (% of Memory CD4 T cells)`, y = `percentTfh`, shape = stage, color = stage, fill = stage)) +
  geom_point(size = 0.8, show.legend = FALSE) +
  geom_smooth(method = lm, alpha = 0.3, se = TRUE, linewidth = 0.5) +
  stat_cor(
    aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")), p.digits = 3,
    label.x = 1,
    label.y = c(50, 45, 40),
    size = 6 / ggplot2:::.pt,
    show.legend = FALSE
  ) +
  scale_color_manual(values = COLORS$stage) +
  scale_fill_manual(values = COLORS$stage) +
  coord_cartesian(clip = "off") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  theme_classic() +
  theme(
    text = element_text(size = BASEPTFONTSIZE, color = "#000000"), 
    axis.title = element_text(size = BASEPTFONTSIZE, color = "#000000"), 
    axis.line = element_line(linewidth = BASELINEWIDTHSIZE),
    axis.text = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    legend.title = element_blank(),
    legend.text = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    legend.key.height = unit(2, "pt"),
    legend.key.width = unit(10, "pt"),
    plot.margin = margin(t = -20),
    plot.title = element_text(color = "#ffffff00", size = BASEPTFONTSIZE + 2),
    legend.position = "bottom",
    legend.box.margin = margin(t = -5, b = 0),
    legend.direction = "horizontal") +
  labs(
    title = "",
    shape = "Stage",
    x = "GC-TFH (% of Mem CD4 T) by flow",
    y = "% Tfh by scDOGMAseq")


################################################################################
# fragment coverage
################################################################################
hxb2Tracker <- array(0L, dim = c(9719, 2))
stageDict <- projUnhashedClean$stage
names(stageDict) <- projUnhashedClean$cellNames

for (i in seq_along(haystackResult$filteredviralFrags$sample)) {
  startPos <- as.numeric(haystackResult$filteredviralFrags[i, "startBp"])
  endPos <- as.numeric(haystackResult$filteredviralFrags[i, "endBp"])
  cbc <- haystackResult$filteredviralFrags[i, "newCbc"]
  
  if (i %% 1000 == 0) {
    message(i)
  }
  
  if (is.na(stageDict[cbc])) {
    next
  }
  
  if (stageDict[cbc] == "Viremic") {
    hxb2Tracker[startPos:endPos, 1] <- hxb2Tracker[startPos:endPos, 1] + 1
  } else if (stageDict[cbc] == "ART") {
    hxb2Tracker[startPos:endPos, 2] <- hxb2Tracker[startPos:endPos, 2] + 1
  }
}

hxb2GeneAnnots <- read.csv("viralGenomeAnnotations/hxb2.csv") %>%
  mutate(annotation = factor(annotation, levels = rev(unique(annotation)))) %>%
  mutate(annotation2 = ifelse(annotation == "LTR", paste0(annotation, row_number()), annotation))

geneAnnot <- ggplot(hxb2GeneAnnots, aes(x = startPos, xend = endPos, y = annotation, yend = annotation)) +
  geom_line(aes(group = annotation2), color = "#CCCCCC") +
  geom_segment(linewidth = 1.5) +
  scale_x_continuous(expand = c(0, 0), limits = c(0,NA)) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  labs(y = "Annotation") +
  theme(
    axis.text.y = element_text(size = BASEPTFONTSIZE - 2),
    axis.text.x = element_text(size = BASEPTFONTSIZE),
    axis.title.x = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_line(linewidth = BASELINEWIDTHSIZE),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(size = BASEPTFONTSIZE, color = "#000000"))

hxb2Df <- data.frame(
  bp = rep(c(1:9719), times = 2),
  counts = as.numeric(c(hxb2Tracker[, 1], hxb2Tracker[, 2])),
  stage = rep(c("Viremic", "ART"), each = 9719)
) %>%
  dplyr::group_by(stage) %>%
  dplyr::mutate(
    density = counts / sum(counts),
    stage = factor(stage, levels = c("Viremic", "ART")))

covgPlot <- ggplot(hxb2Df, aes(x = bp, color = stage, y = density)) +
  geom_line() +
  theme_classic() + 
  labs(y = "Density", x = "HIV proviral position (bp)", color = "Stage") +
  scale_color_manual(values = COLORS$stage) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(limits = c(0, 9719), expand = expansion(mult = c(0, 0))) +
  theme(
    legend.title = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    legend.text = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    axis.text.y = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    axis.text.x = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    axis.title.x = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    axis.line = element_line(linewidth = BASELINEWIDTHSIZE),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(size = BASEPTFONTSIZE, color = "#000000"))

fullCovgPlot <- covgPlot / geneAnnot + plot_layout(heights = c(3.5, 1))


################################################################################
# entropy calcs
################################################################################
props <- data.frame(
  cluster = projUnhashedClean$manualAnnot,
  stage = projUnhashedClean$stage,
  donor = projUnhashedClean$donorClean2,
  hivPos = projUnhashedClean$hivPos
) %>%
  dplyr::mutate(cluster = factor(cluster, levels = unique(cluster))) %>%
  dplyr::mutate(stage = factor(stage, levels = c("HIV-", "Viremic", "ART"))) %>%
  dplyr::group_by(stage, donor, cluster, hivPos, .drop = TRUE) %>%
  dplyr::tally(name = "nCells") %>%
  dplyr::group_by(stage, donor, hivPos) %>%
  dplyr::mutate(propStage = nCells / sum(nCells))

props2 <- props %>%
  dplyr::group_by(stage, donor, hivPos) %>%
  dplyr::summarize(
    totalProp = sum(propStage),
    totalSpecies = n(),
    shannon = -1 * sum(propStage * log(propStage)),
    gs = 1 - sum(propStage ** 2),
    gs2 = sum(propStage * (1 - propStage)),
    rgs = n() * sum(propStage * (1 - propStage)),
    rgs2 = n() * (1 - sum(propStage ** 2))) %>%
  dplyr::mutate(hivPos = factor(hivPos, levels = c("HIV-", "HIV+"))) %>%
  dplyr::arrange(stage, hivPos) %>%
  dplyr::mutate(xaxislbl = paste0(hivPos, "_", stage)) %>%
  dplyr::mutate(xaxislbl = factor(xaxislbl, levels = unique(xaxislbl)))

propStats <- ggpubr::compare_means(rgs2 ~ xaxislbl,
  props2, paired = FALSE, p.adjust.method = "holm")

propStats$p.sym <- pValSymnum(propStats$p.adj)
propStats <- propStats %>%
  dplyr::filter(p.sym != "ns") %>%
  dplyr::mutate(p.position = 30 + 3 * (row_number() %% 7)) %>%
  dplyr::select(-method, -p.signif, -p.format)

pRGS <- ggplot(props2, aes(x = xaxislbl, y = rgs2)) +
  geom_boxplot(aes(color = stage), outliers = FALSE) +
  geom_point(size = 0.8, alpha = 0.8, color = "#000000", pch = 21, fill = "#ffffff") +
  guides(x = legendry::guide_axis_nested(
    key = "_",
    type = "bracket",
    angle = 0,
    levels_text = list(
      NULL,
      element_text(angle = 0, hjust = 0.5, vjust = 1, color = "#000000", size = BASEPTFONTSIZE)))) +
  scale_color_manual(values = COLORS$stage) +
  ggpubr::stat_pvalue_manual(
    propStats,
    label = "{p.sym}", y.position = "p.position",
    size = BASEFONTSIZE * 0.7, vjust = 0.4, tip.length = 0.007, hide.ns = TRUE) +
  labs(y = "Rich-Gini-Simpson\nQuadratic Index") +
  theme_classic() +
  theme(
    axis.line = element_line(linewidth = BASELINEWIDTHSIZE, color = "#000000"),
    axis.text.x = element_text(color = "#000000", size = BASEPTFONTSIZE),
    axis.text.y = element_text(color = "#000000", size = BASEPTFONTSIZE),
    axis.title.x = element_blank(),
    text = element_text(color = "#000000", size = BASEPTFONTSIZE),
    legend.position = "blank",
    strip.background = element_blank())


################################################################################
# scope data from TEACH-B
################################################################################
teachDf_base <- read.csv("external_data/teach_b_scope.csv", check.names = FALSE)
teachDf <- teachDf_base %>%
  dplyr::select(ID, timepoint, `% vDNA+ in BCF`, `% vDNA+ EF`) %>%
  tidyr::pivot_longer(cols = c(`% vDNA+ in BCF`, `% vDNA+ EF`), names_to = "param", values_to = "percent_vDNA") %>%
  dplyr::mutate(param = ifelse(grepl("BCF", param), "B cell follicle", "Extra follicular")) %>%
  dplyr::mutate(param = factor(param, levels = c("Extra follicular", "B cell follicle")))

teachDf2 <- teachDf_base %>%
  dplyr::select(ID, timepoint, `vDNA+/10^5 cells`) %>%
  dplyr::rename(dna = `vDNA+/10^5 cells`) %>%
  dplyr::mutate(timepoint = ifelse(timepoint == 0, "0", "24")) %>%
  dplyr::mutate(timepoint = factor(timepoint, levels = c("0", "24")))


figTeachDnaStats <- ggpubr::compare_means(
  dna ~ timepoint,
  p.adjust.method = "holm",
  paired = TRUE,
  method = "wilcox.test",
  data = teachDf2)
figTeachDnaStats$p.sym <- pValSymnum(figTeachDnaStats$p.adj)


figTeachDna <- ggplot(teachDf2, aes(x = timepoint, y = dna)) +
  geom_boxplot(outliers = FALSE) +
  geom_point() +
  geom_line(aes(group = ID)) +
  theme_classic() +
  ggpubr::stat_pvalue_manual(figTeachDnaStats, label = "{p.sym}",
    y.position = 6.5, bracket.shorten = 0.4, size = BASEFONTSIZE * 0.85, vjust = 0) +
  labs(
    y = "vDNA+/10<sup>5</sup> cells",
    x = "Timepoint (weeks)") +
  theme(
    axis.title.y = element_markdown(color = "#000000", size = BASEPTFONTSIZE),
    axis.text = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    text = element_text(size = BASEPTFONTSIZE, color = "#000000"), 
    axis.title = element_text(size = BASEPTFONTSIZE, color = "#000000"), 
    axis.line = element_line(linewidth = BASELINEWIDTHSIZE),
    plot.title = element_text(color = "#ffffff00", size = BASEPTFONTSIZE + 2))


figTeachCompStats <- ggpubr::compare_means(
  percent_vDNA ~ param,
  p.adjust.method = "holm",
  method = "wilcox.test",
  data = teachDf)
figTeachCompStats$p.sym <- pValSymnum(figTeachCompStats$p.adj)

figTeach <- ggplot(teachDf) +
  geom_boxplot(aes(x = param, y = percent_vDNA), outliers = FALSE, fill = NA, width = 0.5) +
  geom_jitter(aes(x = param, y = percent_vDNA, fill = as.character(timepoint)),
    height = 0, width = 0.3, shape = 21, alpha = 0.8, color = "#000000") +
  theme_classic() +
  scale_color_manual(values = c("blue", "red")) +
  scale_y_continuous(limits = c(0, 110), expand = expansion(mult = c(0.05, 0.05)), breaks = seq(0, 100, by = 20)) +
  ggpubr::stat_pvalue_manual(figTeachCompStats, label = "{p.sym}",
    y.position = 110, bracket.shorten = 0.4, size = BASEFONTSIZE * 0.85, vjust = 0) +
  labs(
    y  = "% vDNA+",
    # x = "Anatomical location",
    fill = "Timepoint (weeks)") +
  theme(
    axis.text = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    text = element_text(size = BASEPTFONTSIZE, color = "#000000"), 
    axis.title.y = element_text(size = BASEPTFONTSIZE, color = "#000000"), 
    axis.title.x = element_blank(), 
    axis.line = element_line(linewidth = BASELINEWIDTHSIZE),
    legend.title = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    legend.text = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    legend.key.height = unit(2, "pt"),
    legend.key.width = unit(2, "pt"),
    plot.title = element_text(color = "#ffffff00", size = BASEPTFONTSIZE + 2),
    legend.position = "bottom",
    legend.box.margin = margin(t = -5, b = 0),
    legend.direction = "horizontal")

range(unlist(teachDf[teachDf$param == "B cell follicle", "percent_vDNA"]))
median(unlist(teachDf[teachDf$param == "B cell follicle", "percent_vDNA"]))

range(unlist(teachDf[teachDf$param == "Extra follicular", "percent_vDNA"]))
median(unlist(teachDf[teachDf$param == "Extra follicular", "percent_vDNA"]))

tissue_tiff <- image_read("external_data/New Figure - TEACH_110_BL_CD4_vDNA_CD20_01(3).tif")
tissue_gg <- image_ggplot(tissue_tiff)

rightSidePanel <- pRGS / 
  figTfh / 
  (figTfhCorr + theme(plot.tag = element_text(margin = margin(t = -30)))) / 
  (figTeach + theme(plot.tag = element_text(margin = margin(t = -70)))) +
  plot_layout(heights = c(1, 2, 1, 1))


####
# PD1 stuff
####
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
  hivPos = projUnhashedClean$hivPos,
  stage = projUnhashedClean$stage
)
rownames(atacImputedDf) <- atacImputedDf$cbc
atacImputedDf[colnames(atacGeneScoresImputed), selectGenesForAtac] <- t(atacGeneScoresImputed)


atacImputedDf <- atacImputedDf %>%
  pivot_longer(cols = -any_of(c("cbc", "manualAnnot", "hivPos", "stage")), names_to = "feature", values_to = "exp") %>%
  dplyr::mutate(stage = factor(stage, levels = names(COLORS$stage))) %>%
  dplyr::mutate(hivPos = factor(hivPos, levels = c("HIV+", "HIV-")))

supAtacPd1 <- atacImputedDf %>%
  dplyr::filter(feature == "PDCD1") %>%
  ggplot(aes(y = hivPos, x = exp)) +
  ggridges::geom_density_ridges(rel_min_height = 0.01, scale = 0.75, alpha = 0.8, aes(fill = hivPos)) +
  geom_vline(xintercept = -Inf, linewidth = BASELINEWIDTHSIZE, color = "#000000") +
  geom_hline(yintercept = -Inf, linewidth = BASELINEWIDTHSIZE, color = "#000000") +
  facet_wrap(~ stage, ncol = 1) +
  coord_cartesian(clip = "off") +
  scale_fill_manual(values = c("HIV+" = "#ff0000", "HIV-" = "#999999")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = BASEPTFONTSIZE - 1, color = "#000000", angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = BASEPTFONTSIZE - 1, color = "#000000"),
    strip.text = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    strip.background = element_rect(color = "#00000000", fill = "#00000000"),
    legend.position = "blank",
    panel.grid.major = element_line(color = "#cfcfcf50"),
    axis.line = element_blank(),
    axis.title.x = element_text(size = BASEPTFONTSIZE - 1, color = "#000000"),
    axis.title.y = element_blank()) +
  labs(
    x = expression("imputed log2(counts + 1) for" ~ italic("PDCD1")))

pd1_rna_df <- FetchData(seuMerged, vars = c("PDCD1", "A0088", "manualAnnot", "stage", "hivPos"), layer = "data")

pd1_rna_df2 <- pd1_rna_df %>%
  dplyr::mutate(hivPos = factor(hivPos, levels = c("HIV-", "HIV+"))) %>%
  dplyr::filter(stage != "HIV-") %>%
  dplyr::mutate(stage = ifelse(stage == "Chronic", "Viremic", stage)) %>%
  dplyr::mutate(stage = factor(stage, levels = c("Viremic", "ART"))) %>%
  tidyr::complete(stage, hivPos)

supRnaPd1 <- ggplot(pd1_rna_df2, aes(x = stage, y = PDCD1, fill = hivPos)) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    colour = "black",
    linewidth = 0.3,
    width = 0.75,
    position = position_dodge(width = 0.75)) +
  rasterize(geom_jitter(
    position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.3),
    size = 0.2,
    alpha = 0.3, show.legend = FALSE), dpi = 300) +
  scale_fill_manual(values = c("HIV+" = "#ff0000", "HIV-" = "#999999")) +
  labs(
    title = "*PDCD1* RNA",
    y = "Normalized expression level",
    x = "Stage") +
  guides(fill = guide_legend(override.aes = list(size = 1))) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    axis.title.x = element_text(size = BASEPTFONTSIZE - 1),
    axis.title.y = element_text(size = BASEPTFONTSIZE - 1),
    axis.text.x = element_text(size = BASEPTFONTSIZE - 1, color = "#000000"),
    axis.text.y = element_text(size = BASEPTFONTSIZE - 1, color = "#000000"),
    legend.position = "inside",
    legend.position.inside = c(0.85, 0.85),
    legend.key.size = unit(0.3, units = "cm"),
    plot.title = element_markdown(hjust = 0.5, size = BASEPTFONTSIZE))



################################################################################
# final compilation
################################################################################
layout <- c(
  area(1, 1, 2, 4), # umaps
  area(1, 5, 2, 6), # upset plots
  area(3, 1, 5, 6), # dna lol
  area(6, 1, 8, 6), # rna lol
  area(1, 7, 8, 9) # right side
)

p <- wrap_elements(full = figA) +
  wrap_elements(full = figB) +
  wrap_elements(full = figDnaLol) +
  wrap_elements(full = figRnaLol) +
  free(rightSidePanel) +
  plot_annotation(tag_levels = list(c("A", "B", "D", "E", "C", "F", "G", "H"))) +
  plot_layout(design = layout) &
  plotTagTheme

saveFinalFigure(
  plot = p,
  fn = "fig2_v3",
  prefixDir = "figures/outs",
  devices = c("png", "pdf"),
  gwidth = JCI_MAX_WIDTH,
  gheight = JCI_MAX_HEIGHT,
  units = "cm",
  addTimestamp = TRUE)


supLayout <- c(
  area(1, 1, 6, 3),
  area(1, 4, 6, 6),
  area(7, 1, 9, 6)
  # area(10, 1, 12, 6),
  # area(1, 5, 2, 6),
  # area(3, 5, 6, 6)
)

supP <- wrap_elements(full = supFigUpsetVir) +
  wrap_elements(full = supFigUpsetArt) +
  wrap_elements(full = fullCovgPlot) +
  plot_annotation(tag_levels = "A") +
  plot_layout(design = supLayout) &
  plotTagTheme

saveFinalFigure(
  plot = supP,
  fn = "supfig3-1_v4",
  prefixDir = "figures/outs",
  devices = c("png", "pdf"),
  gwidth = 6.5,
  gheight = 8,
  addTimestamp = TRUE)

supLayout2 <- c(
  area(1, 1, 3, 4),
  area(4, 1, 5, 2),
  area(4, 3, 5, 4),
  area(1, 5, 2, 6),
  area(3, 5, 5, 6)
)


supP2 <- wrap_elements(full = flow_sort_strat_gg) +
  wrap_elements(full = supAtacPd1) +
  wrap_elements(full = supRnaPd1) +
  wrap_elements(full = figTeachDna) +
  wrap_elements(panel = tissue_gg) +
  plot_annotation(tag_levels = "A") +
  plot_layout(design = supLayout2) &
  plotTagTheme

saveFinalFigure(
  plot = supP2,
  fn = "supfig3-2_v1",
  prefixDir = "figures/outs",
  devices = c("png", "pdf"),
  gwidth = 8,
  gheight = 6,
  addTimestamp = TRUE)

