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


# A HIV infection UMAP (separated by stage)
# B HIV RNA/DNA upset plot (separated by stage)
# C Cell counts by cell identity (separated by donor, stage) for HIV DNA+/- RNA+
# D Cell counts by cell identity (separated by donor, stage) for HIV DNA+ RNA-
# E Correlation with flow data
# F? Gonzalo data (ca HIV DNA by subset) 

################################################################################
# umaps
################################################################################
projUnhashedClean$stage <- ifelse(projUnhashedClean$stage == "Chronic", "Viremic", projUnhashedClean$stage)

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
    strip.text = element_text(size = BASEPTFONTSIZE, color = "#000000", margin = margin(b = 15)),
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
  haystackOut = projUnhashedClean$haystackOut
) %>%
  dplyr::mutate(cluster = factor(cluster, levels = unique(cluster))) %>%
  dplyr::filter(haystackOut) %>%
  dplyr::mutate(stage = factor(stage, levels = c("Viremic", "ART"))) %>%
  dplyr::group_by(stage, cluster, .drop = FALSE) %>%
  dplyr::tally(name = "nCells") %>%
  dplyr::group_by(stage) %>%
  dplyr::mutate(percentStage = nCells / sum(nCells) * 100) %>%
  dplyr::mutate(percentStageLabel = paste0(round(percentStage, digits = 2), "%")) %>%
  dplyr::arrange(desc(percentStage), .by_group = TRUE) %>%
  dplyr::mutate(cluster = factor(cluster, levels = unique(cluster))) %>%
  dplyr::filter(nCells != 0)
  

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
  geom_point() +
  geom_text(aes(label = percentStageLabel, x = percentStage + 1), hjust = 0, size = BASEFONTSIZE) +
  facet_wrap(~ stage) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.2))) +
  scale_color_manual(values = COLORS$stage) +
  labs(x = "% of HIV DNA+ cells per stage")

# savePlot(figC, fn = "croi_hivDnaLollipopByStage", gheight = 6, gwidth = 8, devices = "png")


# fig d RNA lollipop
figDData <- data.frame(
  cluster = projUnhashedClean$manualAnnot,
  stage = projUnhashedClean$stage,
  haystackOut = projUnhashedClean$hivRNA
) %>%
  dplyr::filter(haystackOut) %>%
  dplyr::mutate(stage = factor(stage, levels = c("Viremic", "ART"))) %>%
  dplyr::group_by(stage, cluster) %>%
  dplyr::tally(name = "nCells") %>%
  dplyr::group_by(stage) %>%
  dplyr::mutate(percentStage = nCells / sum(nCells) * 100) %>%
  dplyr::mutate(percentStageLabel = paste0(round(percentStage, digits = 2), "%")) %>%
  dplyr::arrange(desc(percentStage), .by_group = TRUE) %>%
  dplyr::mutate(cluster = factor(cluster, levels = levels(figCData$cluster)))

figRnaLol <- ggplot(figDData, aes(x = percentStage, y = cluster, color = stage)) +
  theme_classic() +
  fig2LollipopTheme +
  geom_linerange(aes(xmin = 0, xmax = percentStage), linetype = "dotted") +
  geom_point() +
  geom_text(aes(label = percentStageLabel, x = percentStage + 1), hjust = 0, size = BASEFONTSIZE) +
  facet_wrap(~ stage) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.2))) +
  scale_color_manual(values = COLORS$stage) +
  labs(x = "% of HIV RNA+ cells per stage")

# savePlot(figD, fn = "croi_hivRnaLollipopByStage", gheight = 6, gwidth = 8, devices = c("png"))


################################################################################
# fig E pairwise distance
################################################################################
# set.seed(42)
# 
# distAverage <- function(mat) {
#   nOnes <- matrix(data = 1, nrow = nrow(mat))
# 
#   avgDist <- 1 / (nrow(mat) * (nrow(mat) - 1)) * t(nOnes) %*% mat %*% nOnes
#   return(as.numeric(avgDist))
# }
# 
# distAverageBoot <- function(sam, full) {
#   distMat <- fields::rdist(full[sam, ])
# 
#   return(distAverage(distMat))
# }
# 
# 
# 
# # for donors with more than 10 infected cells, what was the average distance
# donorCounts <- data.frame(
#   cbc = projUnhashedClean$cellNames,
#   donor = projUnhashedClean$donorClean2,
#   hivPos = projUnhashedClean$hivPos == "HIV+",
#   stage = projUnhashedClean$stage
# ) %>%
#   group_by(donor, hivPos) %>%
#   tally(name = "nCells") %>%
#   filter((hivPos & nCells >= 10) | (grepl("U", donor) & !hivPos)) %>%
#   as.data.frame(.)
# 
# averageHivPos <- round(mean(donorCounts[donorCounts$hivPos, "nCells"]))
# 
# 
# distResults <- lapply(1:nrow(donorCounts), function(i) {
#   donor <- as.character(donorCounts[i, "donor"])
#   nPosCells <- as.integer(donorCounts[i, "nCells"])
#   print(donor)
# 
#   # uninfected, use the average hiv pos for subsetting
#   if (grepl("U", donor)) {
#     donorNegCells <- projUnhashedClean$cellNames[projUnhashedClean$donorClean2 == donor & projUnhashedClean$hivPos == "HIV-"]
#     avgPosDist <- NA
#     donorPosCells <- c(1:averageHivPos)
# 
#   } else {
#     donorPosCells <- projUnhashedClean$cellNames[projUnhashedClean$donorClean2 == donor & projUnhashedClean$hivPos == "HIV+"]
#     donorNegCells <- projUnhashedClean$cellNames[projUnhashedClean$donorClean2 == donor & projUnhashedClean$hivPos == "HIV-"]
# 
#     posDist <- fields::rdist(projUnhashedClean@reducedDims$HarmonyLsiCleanForSeurat$matRD[donorPosCells, ])
#     avgPosDist <- distAverage(posDist)
#   }
# 
# 
#   negDistDf <- data.frame(num = 1:100) %>%
#     group_by(num) %>%
#     dplyr::mutate(means = distAverageBoot(sample(donorNegCells, length(donorPosCells), replace = TRUE),
#       projUnhashedClean@reducedDims$HarmonyLsiCleanForSeurat$matRD[donorNegCells, ]))
# 
#   df <- data.frame(
#     donor = donor,
#     nCells = nPosCells,
#     avgPosDist = avgPosDist,
#     avgNegDist = mean(negDistDf$means)
#   )
# 
#   return(df)
# })
# 
# distResultsCompiled <- bind_rows(distResults) %>%
#   dplyr::mutate(stage = case_when(
#     grepl("A", donor) ~ "ART",
#     grepl("C", donor) ~ "Viremic",
#     grepl("U", donor) ~ "HIV-"
#   )) %>%
#   tidyr::pivot_longer(cols = c("avgPosDist", "avgNegDist"), names_to = "distType", values_to = "dist") %>%
#   dplyr::mutate(distType = ifelse(distType == "avgPosDist", "HIV+", "HIV-")) %>%
#   dplyr::mutate(stage = factor(stage, levels = c("HIV-", "Viremic", "ART"))) %>%
#   dplyr::filter(!is.na(dist)) %>%
#   dplyr::arrange(distType, donor)
# 
# distResultsCompiledStats <- ggpubr::compare_means(dist ~ distType, distResultsCompiled, group.by = "stage", paired = TRUE)
# distResultsCompiledStats$p.sym <- pValSymnum(distResultsCompiledStats$p.adj)
# 
# figDist <- ggplot(distResultsCompiled, aes(x = distType, y = dist)) +
#   geom_boxplot(outliers = FALSE, fill = NA, aes(color = stage), show.legend = FALSE) +
#   geom_line(aes(group = donor), color = "#1f1f1f85") +
#   geom_point(alpha = 0.8) +
#   theme_classic() +
#   labs(x = "Cell infection status", y = "Mean Euclidean distance") +
#   facet_grid(~ stage, scales = "free", space = "free") +
#   ggpubr::stat_pvalue_manual(
#     distResultsCompiledStats, label = "{p.sym}", y.position = 2.2, size = BASEFONTSIZE * 0.85, vjust = 0) +
#   scale_color_manual(values = COLORS$stage) +
#   scale_y_continuous(limits = c(NA, NA), expand = expansion(mult = c(0.05, 0.2))) +
#   theme(
#     # plot.margin = margin(b = 40),
#     axis.line = element_line(linewidth = BASELINEWIDTHSIZE),
#     strip.background = element_rect(fill = "#ffffff00", color = "#ffffff00"),
#     strip.text = element_text(size = BASEPTFONTSIZE, color = "#000000"),
#     axis.title = element_text(size = BASEPTFONTSIZE, color = "#000000"),
#     axis.text = element_text(size = BASEPTFONTSIZE, color = "#000000")
#   )

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
  dplyr::mutate(isTfh = ifelse(isTfh, "Tfh", "Other")) %>%
  dplyr::mutate(isTfh = factor(isTfh, levels = c("Other", "Tfh"))) %>%
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
  dplyr::mutate(xlbl = paste(infection2, isTfh))

figTfhCompStats <- ggpubr::compare_means(
  propOfHIVPos ~ xlbl,
  group.by = c("infection", "stage"),
  p.adjust.method = "bonferroni",
  method = "wilcox.test",
  data = tmp)
figTfhCompStats$p.sym <- pValSymnum(figTfhCompStats$p.adj)


figTfh <- ggplot(tmp, aes(x = xlbl, y = propOfHIVPos * 100)) +
  geom_boxplot(outlier.colour = NA, aes(color = stage), show.legend = FALSE) +
  geom_point(alpha = 0.8, size = 0.8) +
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
    axis.line = element_line(linewidth = BASELINEWIDTHSIZE),
    axis.text.y = element_text(color = "#000000", size = BASEPTFONTSIZE),
    axis.text.x = element_text(color = "#000000", size = BASEPTFONTSIZE - 2),
    panel.grid = element_blank(),
    panel.background = element_rect(color = "#00000000", fill = "#00000000"),
    plot.background = element_rect(color = "#00000000", fill = "#00000000"),
    strip.text = element_text(color = "#000000", size = BASEPTFONTSIZE),
    axis.title = element_text(color = "#000000", size = BASEPTFONTSIZE),
    strip.background = element_rect(color = "#ffffff00", fill = "#ffffff00")) +
  labs(
    x = "Cell phenotype",
    y = "% of donor's HIV+ cells")

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
# tfh corr
################################################################################
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
    size = BASEFONTSIZE * 0.85,
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
    legend.key.width = unit(2, "pt"),
    plot.margin = margin(t = -20),
    plot.title = element_text(color = "#ffffff00", size = BASEPTFONTSIZE + 2),
    legend.position = "bottom",
    legend.box.margin = margin(t = -5, b = 0),
    legend.direction = "horizontal") +
  labs(
    title = "",
    shape = "Stage",
    x = "GC-TFH (% of Memory CD4 T cells) from flow",
    y = "% Tfh from scDOGMAseq")

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
  props2, paired = FALSE)

propStats$p.sym <- pValSymnum(propStats$p.adj)
propStats <- propStats %>%
  dplyr::filter(p.sym != "ns") %>%
  dplyr::mutate(p.position = 30 + 3 * (row_number() %% 7)) %>%
  dplyr::select(-method, -p.signif, -p.format)

pRGS <- ggplot(props2, aes(x = xaxislbl, y = rgs2)) +
  geom_boxplot(aes(color = stage), outlier.shape = NA) +
  geom_point(size = 0.8, alpha = 0.8, color = "#000000") +
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
    axis.title.x = element_blank(),
    text = element_text(color = "#000000", size = BASEPTFONTSIZE),
    legend.position = "blank",
    strip.background = element_blank())


################################################################################
# scope data from TEACH-B
################################################################################
teachDf <- read.csv("external_data/teach_b_scope.csv", check.names = FALSE)
teachDf <- teachDf %>%
  dplyr::select(ID, timepoint, `% vDNA+ in BCF`, `% vDNA+ EF`) %>%
  tidyr::pivot_longer(cols = c(`% vDNA+ in BCF`, `% vDNA+ EF`), names_to = "param", values_to = "percent_vDNA") %>%
  dplyr::mutate(param = ifelse(grepl("BCF", param), "B cell follicle", "Extra follicular")) %>%
  dplyr::mutate(param = factor(param, levels = c("Extra follicular", "B cell follicle")))


figTeachCompStats <- ggpubr::compare_means(
  percent_vDNA ~ param,
  p.adjust.method = "bonferroni",
  method = "wilcox.test",
  data = teachDf)
figTeachCompStats$p.sym <- pValSymnum(figTeachCompStats$p.adj)


figTeach <- ggplot(teachDf) +
  geom_boxplot(aes(x = param, y = percent_vDNA), outliers = FALSE, fill = NA, width = 0.5) +
  geom_jitter(aes(x = param, y = percent_vDNA, fill = as.character(timepoint)),
    height = 0, width = 0.3, shape = 21, color = "#000000") +
  theme_classic() +
  scale_color_manual(values = c("blue", "red")) +
  scale_y_continuous(limits = c(0, 110), expand = expansion(mult = c(0.05, 0.05)), breaks = seq(0, 100, by = 20)) +
  ggpubr::stat_pvalue_manual(figTeachCompStats, label = "{p.sym}",
    y.position = 110, bracket.shorten = 0.4, size = BASEFONTSIZE * 0.85, vjust = 0) +
  labs(
    y  = "% vDNA+",
    x = "Anatomical location",
    fill = "Timepoint (weeks)") +
  theme(
    axis.text = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    text = element_text(size = BASEPTFONTSIZE, color = "#000000"), 
    axis.title = element_text(size = BASEPTFONTSIZE, color = "#000000"), 
    axis.line = element_line(linewidth = BASELINEWIDTHSIZE),
    legend.title = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    legend.text = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    legend.key.height = unit(2, "pt"),
    legend.key.width = unit(2, "pt"),
    plot.title = element_text(color = "#ffffff00", size = BASEPTFONTSIZE + 2),
    legend.position = "bottom",
    legend.box.margin = margin(t = -5, b = 0),
    legend.direction = "horizontal")


rightSidePanel <- pRGS / 
  figTfh / 
  (figTfhCorr + theme(plot.tag = element_text(margin = margin(t = -30)))) / 
  (figTeach + theme(plot.tag = element_text(margin = margin(t = -70)))) +
  plot_layout(heights = c(1, 2, 1, 1))


################################################################################
# final compilation
################################################################################
layout <- c(
  area(1, 1, 2, 4), # umaps
  area(1, 5, 2, 6), # upset plots
  # area(1, 7, 2, 9), # distance
  area(3, 1, 5, 6), # dna lol
  area(6, 1, 8, 6), # rna lol
  area(1, 7, 8, 9) # right side
  # area(3, 7, 4, 9), # tfh comp
  # area(5, 7, 6, 9), # tfh corr
  # area(7, 7, 8, 9) # gonzalo data?
  # area(5, 1, 8, 9) #d
)

p <- wrap_elements(full = figA) +
  wrap_elements(full = figB) +
  wrap_elements(full = figDnaLol) +
  wrap_elements(full = figRnaLol) +
  # wrap_elements(full = figDist) +
  # wrap_elements(full = figTfh) + 
    # theme(plot.tag = element_text(margin = margin(t = -75))) +
  # wrap_elements(full = figTfhCorr) +
  # wrap_elements(full = figTeach) + 
    # theme(plot.tag = element_text(margin = margin(t = -20))) +
  free(rightSidePanel) +
  plot_annotation(tag_levels = list(c("A", "B", "D", "E", "C", "F", "G", "H"))) +
  plot_layout(design = layout) &
  plotTagTheme

saveFinalFigure(plot = p, fn = "fig2_v1", prefixDir = "figures/outs", devices = c("png"), gwidth = 8.5, gheight = 10)


supLayout <- c(
  area(1, 1, 6, 2),
  area(1, 3, 6, 4),
  area(7, 1, 9, 4)
)

supP <- wrap_elements(full = supFigUpsetVir) +
  wrap_elements(full = supFigUpsetArt) +
  wrap_elements(full = fullCovgPlot) +
  plot_annotation(tag_levels = "A") +
  plot_layout(design = supLayout) &
  plotTagTheme

saveFinalFigure(plot = supP, fn = "supfig3_v2", prefixDir = "figures/outs", devices = c("png"), gwidth = 6, gheight = 8.5)

