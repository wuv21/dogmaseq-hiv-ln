# these are all run from project parent directory
library(ggtext)
library(ggrastr)
library(cowplot)
library(ggpubr)
library(ArchR)
library(Seurat)

load("rdata/chronicSub2_deviations_2025-02-26-01-56.Rdata")

set.seed(42)

source("figures/defaultFigureSettings.R")
source("scripts/graphs.R")
source("scripts/colorPalettes.R")

# A Donor UMAP
# B Disease State UMAP
# C Cell Identity UMAP
# D distribution of cells by stage

################################################################################
# fig 1 umap
################################################################################
projUnhashedClean$stage <- ifelse(projUnhashedClean$stage == "Chronic", "Viremic", projUnhashedClean$stage)
# fix naming for CD4 Trm #1
projUnhashedClean$manualAnnot <- ifelse(projUnhashedClean$manualAnnot == "CD4 Trm #1", "CD4 Trm", projUnhashedClean$manualAnnot)
seuMerged$manualAnnot <- ifelse(seuMerged$manualAnnot == "CD4 Trm #1", "CD4 Trm", seuMerged$manualAnnot)


umapDf <- getEmbedding(projUnhashedClean, embedding = "Harmony_UMAP_Combined_v2_05Dist_50Neigh")

stopifnot(sum(rownames(umapDf) == projUnhashedClean$cellNames) == length(projUnhashedClean$cellNames))

umapDf$stage <- projUnhashedClean$stage
umapDf$manualAnnot <- trimws(projUnhashedClean$manualAnnot)
umapDf$tissue <- projUnhashedClean$tissue
umapDf$donor <- projUnhashedClean$donorClean2

umapDf <- dplyr::slice(umapDf, sample(1:n()))

uniqDonors <- unique(umapDf$donor)
donorColorScale <- nColorPalette(n = length(uniqDonors))
names(donorColorScale) <- uniqDonors

uniqAnnots <- sort(unique(umapDf$manualAnnot))
annotColorScale <- nColorPalette(n = length(uniqAnnots))
names(annotColorScale) <- uniqAnnots

umapDf$manualAnnot <- factor(umapDf$manualAnnot, levels = uniqAnnots)
umapDf$stage <- factor(umapDf$stage, levels = names(COLORS$stage))

umapTheme <- theme(
  legend.position = "bottom",
  legend.text = element_text(size = BASEPTFONTSIZE),
  legend.title = element_blank(),
  axis.title = element_text(size = BASEPTFONTSIZE),
  axis.text = element_text(size = BASEPTAXISFONTSIZE),
  legend.spacing.x = unit(BASEPTFONTSIZE / 2, 'points'),
  axis.line = element_line(color = "#000000", linewidth = BASELINEWIDTHSIZE),
  plot.title = element_text(size = BASEPTFONTSIZE, hjust = 0.5),
  plot.title.position = "panel",
  legend.key.size = unit(BASEPTFONTSIZE * 1.1, 'points'),
  legend.background = element_rect(fill = "transparent", colour = NA))

umapPlotThemeNoLeg <- list(
  theme(
    aspect.ratio = 1,
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.text = element_blank())
)

umapLegTheme <- theme(
  legend.position = "bottom",
  legend.direction = "horizontal",
  legend.title = element_text(hjust = 0.5, size = BASEPTFONTSIZE + 1, color = "#000000"),
  legend.key.spacing.y = unit(0.4, 'pt')
)

figA <- ggplot(umapDf,
  aes(x = `HarmonyLsiCleanForSeurat#UMAP_Dimension_1`, y = `HarmonyLsiCleanForSeurat#UMAP_Dimension_2`, color = donor)) +
  rasterize(geom_point(alpha = 0.6, size = 0.25), dpi = 300) +
  guides(colour = guide_legend(override.aes = list(size = 4), nrow = 7)) +
  labs(x = "UMAP 1", y = "UMAP 2", title = "Donor") +
  scale_color_manual(values = donorColorScale) +
  theme_classic() +
  umapTheme
  
figB <- ggplot(umapDf,
  aes(x = `HarmonyLsiCleanForSeurat#UMAP_Dimension_1`, y = `HarmonyLsiCleanForSeurat#UMAP_Dimension_2`, color = stage)) +
  rasterize(geom_point(alpha = 0.6, size = 0.25), dpi = 300) +
  guides(colour = guide_legend(override.aes = list(size = 4), nrow = 7)) +
  labs(x = "UMAP 1", y = "UMAP 2", title = "HIV Stage") +
  scale_color_manual(values = COLORS$stage) +
  theme_classic() +
  umapTheme

figC <- ggplot(umapDf,
  aes(x = `HarmonyLsiCleanForSeurat#UMAP_Dimension_1`, y = `HarmonyLsiCleanForSeurat#UMAP_Dimension_2`, color = manualAnnot)) +
  rasterize(geom_point(alpha = 0.6, size = 0.25), dpi = 300) +
  guides(colour = guide_legend(override.aes = list(size = 4), nrow = 7)) +
  scale_color_manual(values = annotColorScale) +
  labs(x = "UMAP 1", y = "UMAP 2", title = "Cell Identity") +
  theme_classic() +
  umapTheme
  
a_legend <- as_ggplot(get_legend(figA +
    guides(colour = guide_legend(title.position = "top", title.hjust = 0.5, override.aes = list(size = 2.5, alpha = 1), byrow =TRUE, nrow = 7)) +
    labs(color = "Legend for (a)") +
    theme(
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(0,0,t = 20,0),
      legend.text = element_text(margin = margin(l = 2, r = 12, unit = "pt"), size = 5.5)) +
    umapLegTheme +
    subplotTheme))

b_legend <- as_ggplot(get_legend(figB +
    guides(colour = guide_legend(title.position = "top", title.hjust = 0.5, override.aes = list(size = 2.5, alpha = 1), byrow =TRUE)) +
    labs(color = "Legend for (b)") +
    theme(
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(b=-20,0,0,0),
      legend.text = element_text(margin = margin(l = 2, r = 12, unit = "pt"), size = 5.5)) +
    umapLegTheme +
    subplotTheme))

c_legend <- as_ggplot(get_legend(figC +
    guides(colour = guide_legend(title.position = "top", title.hjust = 0.5, override.aes = list(size = 2.5, alpha = 1), byrow = TRUE, nrow = 7)) +
    labs(color = "Legend for (c)") +
    theme(
      legend.margin=margin(0,0,0,0),
      legend.box.margin=margin(b=-30, l=10, t=50, r = 10),
      legend.text = element_text(margin = margin(l = 2, r = 12, unit = "pt"), size = 5.5)) +
    umapLegTheme +
    subplotTheme))


figABLegend <- (a_legend / b_legend)

figA <- figA + umapPlotThemeNoLeg
figB <- figB + umapPlotThemeNoLeg
figC <- figC + umapPlotThemeNoLeg

figABC <- figA + figB + figC + plot_annotation(tag_levels = "A") & plotTagTheme

################################################################################
# fig d cluster prop by stage and donor
################################################################################
figDData <- umapDf %>% 
  dplyr::group_by(stage, manualAnnot, donor, tissue) %>%
  dplyr::tally(name = "totalCells") %>%
  dplyr::group_by(donor) %>%
  dplyr::mutate(percentOfDonor = totalCells / sum(totalCells) * 100)

figDStats <- ggpubr::compare_means(percentOfDonor ~ stage, data = figDData, group.by = c("manualAnnot"), test.method = "BH")

figDFC <- figDData %>%
  dplyr::group_by(stage, manualAnnot) %>%
  dplyr::summarize(meanVal = mean(percentOfDonor)) %>%
  tidyr::pivot_wider(names_from = stage, values_from = meanVal) %>%
  dplyr::mutate(
    "U/C" = `HIV-` / Viremic,
    "U/A" = `HIV-` / ART,
    "C/A" = Viremic / ART
  ) %>%
  tidyr::pivot_longer(cols = any_of(contains("/")), names_to = "condition", values_to = "fc")

figDStats <- figDStats %>%
  dplyr::mutate(p.signif = pValSymnum(p.adj, showNs = FALSE)) %>%
  dplyr::mutate(condition = case_when(
    group1 == "HIV-" & group2 == "Viremic" ~ "U/C",
    group1 == "HIV-" & group2 == "ART" ~ "U/A",
    group1 == "Viremic" & group2 == "ART" ~ "C/A",
  ))

figDFC <- figDFC %>%
  dplyr::left_join(figDStats, by = c("manualAnnot", "condition")) %>%
  dplyr::mutate(p.signif = ifelse(p.adj >= 0.05 & p.adj < 0.1, p.adj, p.signif)) %>%
  dplyr::mutate(condition = stringr::str_wrap(condition, width = 8)) %>%
  dplyr::mutate(condition = factor(condition, levels = c("U/C", "U/A", "C/A")))

abbrevForFig <- c(
  "HIV-" = "HIV- (U)",
  "Viremic" = "Viremic (C)",
  "ART" = "ART (A)"
)

figD <- figDData %>%
  ggplot(aes(x = percentOfDonor, y = manualAnnot)) +
  geom_vline(xintercept = -Inf, linewidth = BASELINEWIDTHSIZE, color = "#000000") +
  geom_hline(yintercept = -Inf, linewidth = BASELINEWIDTHSIZE, color = "#000000") +
  geom_boxplot(alpha = 1, outliers = FALSE, width = 0.5) +
  geom_point(aes(fill = stage, shape = tissue), color = "#000000", alpha = 0.5, stroke = 0.5, size = 1.2) +
  facet_wrap(~ stage, nrow = 1, labeller = as_labeller(abbrevForFig)) +
  scale_fill_manual(values = COLORS$stage) +
  scale_shape_manual(values = c("CLN" = 21, "ILN" = 24)) +
  theme_classic() +
  theme(
    axis.line.x = element_blank(),
    strip.background = element_rect(color = "#ffffff", fill = "#ffffff"),
    strip.text = element_text(color = "#000000", size = BASEPTFONTSIZE),
    axis.text = element_text(color = "#000000", size = BASEPTFONTSIZE),
    axis.title.y = element_blank(),
    axis.title.x = element_text(color = "#000000", size = BASEPTFONTSIZE),
    legend.position = "blank",
    panel.grid.major.x = element_blank(),
    plot.margin = margin(t = 30,0,0,0),
    panel.grid.major.y = element_line(color = "#dddddd50"),
    plot.background = element_rect(color = "#ffffff00", fill = NA),
    panel.background = element_rect(color = "#ffffff00", fill = "#ffffff00"),
    panel.grid.minor.x = element_blank()
  ) +
  scale_y_discrete(limits = rev) +
  labs(
    x = "% of each donor's cells")

figDStatPanel <- ggplot(figDFC, aes(x = condition, y = manualAnnot)) +
  geom_tile(aes(fill = log2(fc))) +
  geom_text(aes(label = p.signif), vjust = 0.5, size = 1.2) +
  scale_fill_gradient2(
    low = "#40004BFF",
    mid = "#ffffff",
    high = "#00441BFF",
    midpoint = 0,
    guide = guide_colorbar(barwidth = 0.5, title.position = "right")) +
  theme_classic() +
  scale_x_discrete(position = "bottom") +
  theme(
    axis.line.x = element_blank(),
    axis.text = element_text(color = "#000000", size = BASEPTFONTSIZE),
    axis.title.y = element_blank(),
    axis.title.x = element_text(color = "#000000", size = BASEPTFONTSIZE),
    legend.position = "right",
    legend.title = element_text(color = "#000000", size = BASEPTFONTSIZE, angle = -90, vjust = 0.5, hjust = 0.5),
    legend.box.margin = margin(b = 0, l = 0),
    legend.margin = margin(b = 0, l = 0),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(t = 30,0,0,0),
    panel.grid.major.y = element_blank(),
    plot.background = element_rect(color = "#ffffff00", fill = NA),
    panel.background = element_rect(color = "#ffffff00", fill = "#ffffff00"),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_blank(),
  ) +
  scale_y_discrete(limits = rev) +
  labs(x = "Stage comparison", fill = "log2(fold change)")


################################################################################
# sup fig 2
################################################################################
baseGenesForAnnotation_v2 <- read.table("manualAnnot/baseGenesForAnnotation_v2.txt")$V1

seuMerged$manualAnnot <- trimws(seuMerged$manualAnnot)

# base rna
rnaBasePlot <- DotPlot(
  seuMerged,
  dot.scale = 5,
  col.max = 3,
  col.min = -3,
  features = paste0("rna_", baseGenesForAnnotation_v2),
  group.by = "manualAnnot") +
  scale_y_discrete(limits = uniqAnnots) +
  scale_x_discrete(labels = ~ gsub("rna_", "", .x)) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = BASEPTFONTSIZE - 2, color = "#000000", face = "italic"),
    axis.text.y = element_text(size = BASEPTFONTSIZE - 3, color = "#000000"),
    panel.grid.major.y = element_line(color = "#efefef"),
    axis.title = element_blank(),
    legend.key.size = unit(6, "pt"),
    legend.title = element_text(size = BASEPTFONTSIZE - 2, color = "#000000", hjust = 0.5),
    legend.text = element_text(size = BASEPTFONTSIZE - 2, color = "#000000"))

rnaBasePlot$guides$guides$colour$params$title <- "Average\nExpression"
rnaBasePlot$guides$guides$size$params$title <- "% Expressed"


# base atac
baseGenesForAtac <- c("CD3D", "CD69", "CCR7", "CD27", "CD28",
  "IL2RA", "ITGAE", "PDCD1", "TGFB1", "FAS", "PRDM1", "CXCR5", "BCL6",
  "CCR5", "FOXP3", "IL2", "IL4", "STAT3", "ZNF683", "KLF2",
  "GATA3", "RORC", "NFKB1", "TOX")
atacGeneScoreForAnnot_v3 <- tryCatch({
  tmp <- ArchR:::.getMatrixValues(
    ArchRProj = projUnhashedClean,
    name = baseGenesForAtac,
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

atacGeneScoreForAnnotImputed_v3 <- imputeMatrix(
  mat = as.matrix(atacGeneScoreForAnnot_v3),
  imputeWeights = getImputeWeights(projUnhashedClean))

atacGeneScoreForAnnotImputed_v3 <- log2(atacGeneScoreForAnnotImputed_v3 + 1)

atacImputedDf_v3 <- data.frame(
  cbc = projUnhashedClean$cellNames,
  manualAnnot = projUnhashedClean$manualAnnot
)
rownames(atacImputedDf_v3) <- atacImputedDf_v3$cbc
atacImputedDf_v3[colnames(atacGeneScoreForAnnotImputed_v3), baseGenesForAtac] <- t(atacGeneScoreForAnnotImputed_v3)

atacImputedDf_v3 <- atacImputedDf_v3 %>%
  pivot_longer(cols = -any_of(c("cbc", "manualAnnot")), names_to = "feature", values_to = "exp")

atacImputePlot_v3 <- atacImputedDf_v3 %>%
  mutate(manualAnnot = trimws(manualAnnot)) %>%
  ggplot(aes(x = manualAnnot, y = exp)) +
  geom_violin(aes(color = manualAnnot), width = 1, fill = "#ffffff00") +
  facet_wrap(~ feature, scales = "free_y", nrow = 10) +
  scale_color_manual(values = annotColorScale) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = BASEPTFONTSIZE - 2.5, color = "#000000"),
    axis.text.y = element_text(size = BASEPTFONTSIZE - 2.5, color = "#000000"),
    strip.text = element_text(size = BASEPTFONTSIZE - 1, color = "#000000", face = "italic"),
    plot.margin = margin(l = 18),
    strip.background = element_rect(color = "#00000000", fill = "#00000000"),
    legend.position = "blank",
    panel.grid = element_blank(),
    axis.title.y = element_text(size = BASEPTFONTSIZE),
    axis.title.x = element_blank()) +
  labs(
    y = "imputed log2(counts + 1)")

################################################################################
# sup table of ipda
################################################################################
# ipda values are per million cells

ipdaFlow <- read.csv("external_data/ipda_flow_combined_2024-11-22.csv", check.names = FALSE) %>%
  select(Sample, `Intact`, `3' Defective`, `5' Defective`,`Total Proviruses Detected`)

supTableX <- data.frame(
  donor = projUnhashedClean$donorClean2,
  hivStatus = projUnhashedClean$hivStatus,
  manualAnnot = projUnhashedClean$manualAnnot,
  stage = projUnhashedClean$stage
) %>%
  dplyr::mutate(stage = factor(stage, levels = c("HIV-", "Viremic", "ART"))) %>%
  dplyr::group_by(donor, stage) %>%
  dplyr::summarize(
    `Total cells` = n(),
    `HIV DNA+/- RNA+` = sum(hivStatus == "HIV DNA+/- RNA+"),
    `HIV-` = sum(hivStatus == "HIV-"),
    `HIV DNA+ RNA-` = sum(hivStatus == "HIV DNA+ RNA-")
  ) %>%
  dplyr::left_join(ipdaFlow, by = c("donor" = "Sample")) %>%
  dplyr::mutate(
    `HIV+` = `HIV DNA+/- RNA+` + `HIV DNA+ RNA-`,
    `HIV+ (% of non-naive CD8- T cells)` = round(`HIV+` / `Total cells` * 100, 4),
    `Intact\n(% of total CD4+ T cells)` = round(`Intact` / 1e6 * 100, 4),
    `Total Proviruses Detected \n(% of total CD4+ T cells)` = round(`Total Proviruses Detected` / 1e6 * 100, 4),
  ) %>%
  dplyr::select(
    `donor`,
    `Total cells`,
    `stage`,
    `HIV-`,
    `HIV DNA+/- RNA+`,
    `HIV DNA+ RNA-`,`HIV+`,
    `HIV+ (% of non-naive CD8- T cells)`,
    `Intact\n(% of total CD4+ T cells)`,
    `Total Proviruses Detected \n(% of total CD4+ T cells)`) 


################################################################################
# final plot saving
################################################################################

layout <- c(
  area(1, 1, 2, 6), #abc
  area(1, 7, 2, 9), #legend
  area(2, 1, 3, 9), #legend
  area(4, 1, 8, 9) #d
)

p <- wrap_elements(full = figABC, ignore_tag = TRUE) +
  wrap_elements(panel = figABLegend, ignore_tag = TRUE, clip = FALSE) +
  wrap_elements(panel = c_legend, ignore_tag = TRUE, clip = FALSE) +
  wrap_elements(full = figD + figDStatPanel + plot_layout(widths = c(0.85, 0.15)) & 
      theme(plot.background = element_rect(fill = "#ffffff00", color = "#ffffff00"))) + 
  theme(
    plot.tag = element_text(margin = margin(t = 30))) +
  plot_annotation(tag_levels = list(c("D"))) +
  plot_layout(design = layout) &
  plotTagTheme

saveFinalFigure(
  plot = p,
  fn = "fig1_v4",
  prefixDir = "figures/outs",
  devices = c("png", "pdf"),
  gwidth = 7 * 2.54,
  gheight = 7 * 2.54,
  units = "cm",
  addTimestamp = TRUE)


supLayout <- c(
  area(1, 1, 4, 10),
  area(5, 1, 9, 10)
)

supP <- free(rnaBasePlot, type = "panel", side = "lr") +
  free(atacImputePlot_v3, type = "space", side = "lr") +
  plot_annotation(tag_levels = "A") +
  plot_layout(design = supLayout) &
  plotTagTheme

saveFinalFigure(
  plot = supP,
  fn = "supfig2_v4",
  prefixDir = "figures/outs",
  devices = c("png"),
  gwidth = 8,
  gheight = 10,
  addTimestamp = TRUE)
