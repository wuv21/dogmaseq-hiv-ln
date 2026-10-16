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
# tx analysis
################################################################################
seuMerged_art <- subset(seuMerged, subset = stage == "ART" & manualAnnotBin != "Naive/early memory" & manualAnnotBin != "Other")
seuMerged_art$manualAnnotBinWithHiv <- paste0(seuMerged_art$manualAnnotBin, "_", seuMerged_art$hivPos)
seuMerged_art$manualAnnotBinWithHivStatus <- paste0(seuMerged_art$manualAnnotBin, "_", seuMerged_art$hivStatus)

lowPop <- table(seuMerged_art$manualAnnotBinWithHivStatus) < 25
lowPop <- names(lowPop[lowPop])

Idents(seuMerged_art) <- "manualAnnotBinWithHiv"
uniq_deg_hivPosVsNeg <- FindAllMarkers(
  seuMerged_art,
  group.by = "manualAnnotBinWithHiv",
  min.pct = 0.3,
  assay = "RNA",
  only.pos = TRUE
)


Idents(seuMerged_art) <- "manualAnnotBinWithHiv"
uniq_deg_hivPosVsNeg_mast <- FindAllMarkers(
  seuMerged_art,
  group.by = "manualAnnotBinWithHiv",
  min.pct = 0.3,
  test.use = "MAST",
  latent.vars = c("donorClean2"),
  assay = "RNA",
  only.pos = TRUE
)

uniq_deg_hivPosVsNeg_filt <- uniq_deg_hivPosVsNeg_mast %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::group_by(cluster) %>%
  dplyr::slice_max(avg_log2FC, n = 25)

uniq_deg_mat <- AverageExpression(seuMerged_art, features = uniq_deg_hivPosVsNeg_filt$gene, assay = "RNA", group.by = c("manualAnnotBinWithHiv"))
uniq_deg_mat <- uniq_deg_mat$RNA
uniq_deg_mat <- uniq_deg_mat[, order(colnames(uniq_deg_mat))]

colnamesExtract <- str_match(colnames(uniq_deg_mat), "(.*)-(HIV.*)")

figA <- ComplexHeatmap::Heatmap(
  matrix = t(scale(t(uniq_deg_mat))),
  name = " ",
  row_km = 12,
  row_names_gp = gpar(fontsize = 4),
  row_title_gp = gpar(fontsize = 6),
  row_dend_width = unit(8, "points"),
  column_title_gp = gpar(fontsize = 0),
  cluster_columns = FALSE,
  column_split = rep(c(1:6), each = 2),
  show_column_names = FALSE,
  col = circlize::colorRamp2(c(-4, 0, 4), c("#0000ff", "#f7f7f7", "#ff0000")),
  heatmap_legend_param = list(
    labels_gp = gpar(fontsize = 6),
    direction = "horizontal",
    grid_height = unit(0.2, "cm")),
  top_annotation = ComplexHeatmap::HeatmapAnnotation(
    `Phenotype` = anno_block(
      gp = gpar(fill = "#ffffff", col = "#ffffff"),
      height = unit(BASEPTFONTSIZE, "pt"),
      labels = unique(colnamesExtract[, 2]),
      labels_gp = gpar(col = "#000000", fontsize = BASEPTFONTSIZE - 2),
      show_name = FALSE),
    simple_anno_size = unit(0.3, "cm")
  ),
  bottom_annotation = ComplexHeatmap::columnAnnotation(
    `Phenotype` = anno_block(
      gp = gpar(fill = "#ffffff", col = "#ffffff"),
      height = unit(BASEPTFONTSIZE, "pt"),
      labels = unique(colnamesExtract[, 2]),
      labels_gp = gpar(col = "#000000", fontsize = BASEPTFONTSIZE - 2),
      show_name = TRUE),
    `HIV detection` = colnamesExtract[, 3],
    annotation_legend_param = list(
      labels_gp = gpar(fontsize = 6),
      title_gp = gpar(fontsize = 6),
      ncol = 2
    ),
    col = list(
      `HIV detection` = c("HIV+" = "#ff0000", "HIV-" = "#807c7c")
    ),
    annotation_name_gp = gpar(fontsize = 6),
    simple_anno_size = unit(0.3, "cm")
  )
)


################################################################################
# prev reported
################################################################################
Idents(seuMerged_art) <- "manualAnnotBinWithHiv"
art_goi <- c("IKZF3", "SERPINB9", "BACH2", "GZMB", "PRF1", "FCGR2A", "JUN", "FOSB", "NFKB1", "BCL2", "SLC3A2", "ARPP19")
uniqBins <- unique(seuMerged_art$manualAnnotBin)

prevReported_tests <- lapply(uniqBins, function(x) {
  deg <- FindMarkers(
    seuMerged_art,
    features = art_goi,
    ident.1 = paste0(x, "_HIV+"),
    ident.2 = paste0(x, "_HIV-"),
    min.pct = 0
  )
  
  deg$gene <- rownames(deg)
  return(deg)
})
names(prevReported_tests) <- uniqBins
prevReported_tests_compiled <- bind_rows(prevReported_tests, .id = "bin")
rownames(prevReported_tests_compiled) <- NULL

prevReported_tests_compiled$p_val_adj2 <- p.adjust(prevReported_tests_compiled$p_val, method = "BH")


figB <- DotPlot(
  seuMerged_art,
  features = art_goi,
  group.by = "manualAnnotBinWithHiv")

figB$guides$guides$colour$params$title <- "Average\n Expression"
figB$guides$guides$size$params$title <- "% Expressed"

cleanDotPlotY <- function(x) {
  matches <- stringr::str_match(x, pattern = "(.*)_(HIV[+-])")
  colors <- ifelse(matches[, 3] == "HIV-", "#807c7c", "#ff0000")
  
  return(glue("<span style='font-size:6pt'>{matches[, 2]} <span style='color:{colors}'>{matches[, 3]}</span></span>"))
}

figB_mod <- figB +
  coord_flip() +
  theme(
    axis.text = element_text(size = BASEPTFONTSIZE - 2, color = "#000000"),
    axis.title = element_blank(),
    axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1, size = BASEPTFONTSIZE),
    legend.title = element_text(size = BASEPTFONTSIZE - 2, color = "#000000", hjust = 0.5),
    legend.text = element_text(size = BASEPTFONTSIZE - 2, color = "#000000", hjust = 0),
    panel.grid.major.y = element_line(color = "#cfcfcf30"),
    axis.line = element_line(color = "#000000", linewidth = BASELINEWIDTHSIZE),
    legend.margin = margin(l = 0),
    legend.key.size = unit(.5, "lines")
  ) +
  scale_y_discrete(labels = cleanDotPlotY)


# cell death sup fig
# https://pmc.ncbi.nlm.nih.gov/articles/PMC7169302/
proApopGOI <- c(
  "TNFSF10", # trail
  "FASLG", # fasL
  "DISC1",
  "TRADD",
  "FADD",
  "BBC3",
  "BCL2L11", # Bim
  "BMF",
  "BIK",
  "BAD",
  "BID",
  "PMAIP1", # noxa
  "BAX",
  "BAK1", # bak
  "AIFM1", # aif
  "ENDOG",
  "DIABLO"
)

antiApopGOI <- c(
  "BCL2",
  "BCL2L1", # Bcl-xl
  "CFLAR", # c-FLIP
  "BIRC5", # survivin
  "XIAP"
)

supFigA <- DotPlot(
  seuMerged_art,
  features = c(proApopGOI, antiApopGOI),
  group.by = "manualAnnotBinWithHiv")


supFigA_mod <- supFigA$data %>%
  dplyr::mutate(category = ifelse(features.plot %in% antiApopGOI, "Anti-apoptosis", "Pro-apoptosis")) %>%
  dplyr::mutate(new.features.plot = paste0(features.plot, "_", category)) %>%
  dplyr::mutate(new.features.plot = factor(new.features.plot, levels = unique(new.features.plot))) %>%
  ggplot() +
  geom_point(aes(x = id, y = new.features.plot, size = pct.exp, color = avg.exp.scaled)) +
  guides(
    y = legendry::guide_axis_nested(
      key = "_",
      type = "bracket",
      angle = 0,
      levels_text = list(
        element_text(face = "italic", color = "#000000", size = 8),
        element_text(angle = 90, hjust = 0.5, vjust = 0.5, color = "#000000", size = 8)))) +
  scale_color_distiller(palette = "PuRd", direction = 1) +
  theme_classic() +
  theme(
    axis.text = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    axis.title = element_blank(),
    axis.text.x = element_markdown(angle = 45, hjust = 1, vjust = 1, size = BASEPTFONTSIZE),
    legend.title = element_text(size = BASEPTFONTSIZE, color = "#000000", hjust = 0.5),
    legend.text = element_text(size = BASEPTFONTSIZE, color = "#000000", hjust = 0),
    panel.grid.major.y = element_line(color = "#cfcfcf30"),
    axis.line = element_line(color = "#000000", linewidth = BASELINEWIDTHSIZE),
    plot.margin = margin(10, 10, 10, 10),
    legend.key.size = unit(.5, "lines")
  ) +
  scale_x_discrete(labels = cleanDotPlotY) +
  labs(size = "Percent\nExpressed", color = "Average\nExpression")


celldeath_tests <- lapply(uniqBins, function(x) {
  Idents(seuMerged_art) <- "manualAnnotBinWithHiv"
  
  deg <- FindMarkers(
    seuMerged_art,
    features = c(proApopGOI, antiApopGOI),
    ident.1 = paste0(x, "_HIV+"),
    ident.2 = paste0(x, "_HIV-"),
    min.pct = 0
  )
  
  deg$gene <- rownames(deg)
  return(deg)
})
names(celldeath_tests) <- uniqBins
celldeath_tests_compiled <- bind_rows(celldeath_tests, .id = "bin")
rownames(celldeath_tests_compiled) <- NULL

celldeath_tests_compiled$p_val_adj2 <- p.adjust(celldeath_tests_compiled$p_val, method = "BH")


################################################################################
# deg
################################################################################
# 
Idents(seuMerged_art) <- "hivPos"
# for mast, we need to adjust based on how many genes were inputted into
# mast, since not all genes were tested/used.
# hypothesis is then testing on genes with min.pct = 0.3 and logfc.threshold = 0.2
all_deg_mast <- FindMarkers(
  seuMerged_art,
  assay = "RNA",
  ident.1 = "HIV+",
  ident.2 = "HIV-",
  logfc.threshold = 0.1,
  test.use = "MAST",
  latent.vars = c("donorClean2"),
  min.pct = 0.3)

# mast is not good for dea
all_dea_lr <- FindMarkers(
  seuMerged_art,
  features = tsaCatalog[!tsaCatalog$isCtrl, "DNA_ID"],
  assay = "ADT",
  ident.1 = "HIV+",
  ident.2 = "HIV-",
  test.use = "LR",
  latent.vars = c("donorClean2"),
  min.pct = 0.3)

tem_deg_mast <- FindMarkers(
  subset(seuMerged_art, subset = manualAnnotBin == "Tem/Th1"),
  assay = "RNA",
  ident.1 = "HIV+",
  ident.2 = "HIV-",
  test.use = "MAST",
  latent.vars = c("donorClean2"),
  logfc.threshold = 0.1,
  min.pct = 0.3)

tem_dea_lr <- FindMarkers(
  subset(seuMerged_art, subset = manualAnnotBin == "Tem/Th1"),
  features = tsaCatalog[!tsaCatalog$isCtrl, "DNA_ID"],
  assay = "ADT",
  ident.1 = "HIV+",
  ident.2 = "HIV-",
  test.use = "LR",
  latent.vars = c("donorClean2"),
  min.pct = 0.3)

tcm_deg_mast <- FindMarkers(
  subset(seuMerged_art, subset = manualAnnotBin == "Tcm/Th2"),
  assay = "RNA",
  ident.1 = "HIV+",
  ident.2 = "HIV-",
  logfc.threshold = 0.1,
  test.use = "MAST",
  latent.vars = c("donorClean2"),
  min.pct = 0.3)

tcm_dea_lr <- FindMarkers(
  subset(seuMerged_art, subset = manualAnnotBin == "Tcm/Th2"),
  features = tsaCatalog[!tsaCatalog$isCtrl, "DNA_ID"],
  assay = "ADT",
  ident.1 = "HIV+",
  ident.2 = "HIV-",
  test.use = "LR",
  latent.vars = c("donorClean2"),
  min.pct = 0.3)

deLollipop_theme <- list(
  theme_classic(),
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "#807c7c")),
  labs(x = "<span style='font-size:7pt'>Avg log2(FC)<br><span style='color:#807c7c'>&uarr; in HIV-</span> | <span style='color:#ff0000'>&uarr; in HIV+</span></span>"),
  theme(
    axis.line = element_line(linewidth = BASELINEWIDTHSIZE),
    legend.position = "blank",
    axis.text = element_text(size = BASEPTFONTSIZE - 2, color = "#000000"),
    axis.title.y = element_blank(),
    axis.title.x = element_markdown(),
    plot.title = element_text(size = BASEPTFONTSIZE, color = "#000000", hjust = 0.5)
  ))

lollipop_point_size <- 0.75

all_dea_plot <- all_dea_lr %>%
  as.data.frame() %>%
  mutate(DNA_ID = rownames(.)) %>%
  left_join(tsaCatalog[!tsaCatalog$isCtrl, c("DNA_ID", "cleanName")]) %>%
  mutate(p_val_adj2 = p.adjust(p_val, method = "BH")) %>%
  filter(p_val_adj2 < 0.05) %>%
  dplyr::mutate(direction = avg_log2FC > 0) %>%
  dplyr::group_by(direction) %>%
  dplyr::slice_max(abs(avg_log2FC), n = 10) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::mutate(cleanName = factor(cleanName, levels = cleanName)) %>%
  ggplot(aes(y = cleanName, x = avg_log2FC, color = direction)) +
  geom_vline(xintercept = 0) +
  geom_segment(aes(x = 0, xend = avg_log2FC), linetype = "dotted") +
  geom_point(size = lollipop_point_size) +
  deLollipop_theme +
  labs(y = "Surface protein", title = "All cells\nSurface protein")

all_deg_plot <- all_deg_mast %>%
  as.data.frame() %>%
  dplyr::mutate(gene = rownames(.)) %>%
  dplyr::mutate(p_val_adj2 = p.adjust(p_val, method = "BH")) %>%
  dplyr::filter(p_val_adj2 < 0.05) %>%
  dplyr::mutate(direction = avg_log2FC > 0) %>%
  dplyr::group_by(direction) %>%
  dplyr::slice_max(abs(avg_log2FC), n = 10) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::mutate(gene = factor(gene, levels = gene)) %>%
  ggplot(aes(y = gene, x = avg_log2FC, color = direction)) +
  geom_vline(xintercept = 0) +
  geom_segment(aes(x = 0, xend = avg_log2FC), linetype = "dotted") +
  geom_point(size = lollipop_point_size) +
  deLollipop_theme +
  labs(y = "Gene", title = "All cells\nRNA")

fill_sig_dict <- c(
  "TRUE_TRUE" = "#ff0000",
  "TRUE_FALSE" = "#ffffff",
  "FALSE_TRUE" = "#807c7c",
  "FALSE_FALSE" = "#ffffff")

# showing genes that have raw p value < 0.05 but with open circles
# for those have adjusted (BH) p value > 0.05
tem_deg_plot <- tem_deg_mast %>%
  as.data.frame() %>%
  dplyr::mutate(gene = rownames(.)) %>%
  dplyr::mutate(p_val_adj2 = p.adjust(p_val, method = "BH")) %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::mutate(direction = avg_log2FC > 0) %>%
  dplyr::mutate(sig = p_val_adj2 < 0.05) %>%
  dplyr::mutate(directionSigFill = paste(direction, sig, sep = "_")) %>%
  dplyr::group_by(direction) %>%
  dplyr::slice_max(abs(avg_log2FC), n = 10) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::mutate(gene = factor(gene, levels = gene)) %>%
  ggplot(aes(y = gene, x = avg_log2FC, color = direction)) +
  geom_vline(xintercept = 0) +
  geom_segment(aes(x = 0, xend = avg_log2FC), linetype = "dotted") +
  geom_point(aes(fill = directionSigFill), shape = 21, size = lollipop_point_size) +
  geom_point(fill = "#ffffff00", shape = 21, size = lollipop_point_size) +
  deLollipop_theme +
  scale_fill_manual(values = fill_sig_dict) +
  labs(y = "Gene", title = "Tem/Th1 cells\nRNA")

tem_dea_plot <- tem_dea_lr %>%
  as.data.frame() %>%
  mutate(DNA_ID = rownames(.)) %>%
  left_join(tsaCatalog[!tsaCatalog$isCtrl, c("DNA_ID", "cleanName")]) %>%
  dplyr::mutate(p_val_adj2 = p.adjust(p_val, method = "BH")) %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::mutate(direction = avg_log2FC > 0) %>%
  dplyr::mutate(sig = p_val_adj2 < 0.05) %>%
  dplyr::mutate(directionSigFill = paste(direction, sig, sep = "_")) %>%
  dplyr::group_by(direction) %>%
  dplyr::slice_max(abs(avg_log2FC), n = 10) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::mutate(cleanName = factor(cleanName, levels = cleanName)) %>%
  ggplot(aes(y = cleanName, x = avg_log2FC, color = direction)) +
  geom_vline(xintercept = 0) +
  geom_segment(aes(x = 0, xend = avg_log2FC), linetype = "dotted") +
  geom_point(aes(fill = directionSigFill), shape = 21, size = lollipop_point_size) +
  geom_point(fill = "#ffffff00", shape = 21, size = lollipop_point_size) +
  deLollipop_theme +
  scale_fill_manual(values = fill_sig_dict) +
  labs(y = "Surface protein", title = "Tem/Th1 cells\nSurface protein")

tcm_deg_plot <- tcm_deg_mast %>%
  as.data.frame() %>%
  dplyr::mutate(gene = rownames(.)) %>%
  dplyr::mutate(p_val_adj2 = p.adjust(p_val, method = "BH")) %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::mutate(direction = avg_log2FC > 0) %>%
  dplyr::mutate(sig = p_val_adj2 < 0.05) %>%
  dplyr::mutate(directionSigFill = paste(direction, sig, sep = "_")) %>%
  dplyr::group_by(direction) %>%
  dplyr::slice_max(abs(avg_log2FC), n = 10) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::mutate(gene = factor(gene, levels = gene)) %>%
  ggplot(aes(y = gene, x = avg_log2FC, color = direction)) +
  geom_vline(xintercept = 0) +
  geom_segment(aes(x = 0, xend = avg_log2FC), linetype = "dotted") +
  geom_point(aes(fill = directionSigFill), shape = 21, size = lollipop_point_size) +
  geom_point(fill = "#ffffff00", shape = 21, size = lollipop_point_size) +
  deLollipop_theme +
  scale_fill_manual(values = fill_sig_dict) +
  labs(y = "Gene", title = "Tcm/Th2 cells\nRNA")

tcm_dea_plot <- tcm_dea_lr %>%
  as.data.frame() %>%
  mutate(DNA_ID = rownames(.)) %>%
  left_join(tsaCatalog[!tsaCatalog$isCtrl, c("DNA_ID", "cleanName")]) %>%
  dplyr::mutate(p_val_adj2 = p.adjust(p_val, method = "BH")) %>%
  dplyr::filter(p_val < 0.05) %>%
  dplyr::mutate(sig = p_val_adj2 < 0.05) %>%
  dplyr::mutate(direction = avg_log2FC > 0) %>%
  dplyr::mutate(directionSigFill = paste(direction, sig, sep = "_")) %>%
  dplyr::group_by(direction) %>%
  dplyr::slice_max(abs(avg_log2FC), n = 10) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::mutate(cleanName = factor(cleanName, levels = cleanName)) %>%
  ggplot(aes(y = cleanName, x = avg_log2FC, color = direction)) +
  geom_vline(xintercept = 0) +
  geom_segment(aes(x = 0, xend = avg_log2FC), linetype = "dotted") +
  geom_point(aes(fill = directionSigFill), shape = 21, size = lollipop_point_size) +
  geom_point(fill = "#ffffff00", shape = 21, size = lollipop_point_size) +
  deLollipop_theme +
  scale_fill_manual(values = fill_sig_dict) +
  labs(y = "Surface protein", title = "Tcm/Th2 cells\nSurface protein")

################################################################################
# supplemental figures
################################################################################
# ucell
# tricycle

markerSigs <- list(
  "serpinb9" = read.table("gene_lists/serpinb9.txt")$V1,
  "bcl2" = read.table("gene_lists/bcl2.txt")$V1,
  "caspase" = read.table("gene_lists/caspase_apoptosis.txt")$V1
)

seuMerged_art <- UCell::AddModuleScore_UCell(seuMerged_art, features = markerSigs)


supFigUcell <- data.frame(
  hivPos = seuMerged_art$hivPos,
  bin = seuMerged_art$manualAnnotBin,
  `SERPINB9` = seuMerged_art$serpinb9_UCell,
  `BCL2` = seuMerged_art$bcl2_UCell,
  `Caspase` = seuMerged_art$caspase_UCell
) %>%
  tidyr::pivot_longer(c(3:5), names_to = "gene set", values_to = "score")
  
supFigUcellCompiledStats <- ggpubr::compare_means(
  score ~ hivPos,
  group.by = c("gene set", "bin"),
  paired = FALSE, data = supFigUcell)
supFigUcellCompiledStats$p.sym <- pValSymnum(supFigUcellCompiledStats$p.adj)
  
supFigUcellPlot <- supFigUcell %>%
  ggplot(aes(x = hivPos, y = score)) +
  geom_hline(yintercept = -Inf, linewidth = BASELINEWIDTHSIZE) +
  geom_vline(xintercept = -Inf, linewidth = BASELINEWIDTHSIZE) +
  geom_violin(aes(fill = hivPos)) +
  ggpubr::stat_pvalue_manual(
    supFigUcellCompiledStats, label = "{p.sym}",
    y.position = c(0.275), size = BASEFONTSIZE * 0.6, vjust = 0, tip.length = 0.01) +
  facet_grid(`gene set` ~ `bin`) +
  theme_classic() +
  labs(x = "HIV status", y = "UCell score") +
  scale_fill_manual(values = c("HIV-" = "#807c7c", "HIV+" = "#ff0000")) +
  theme(
    legend.position = "blank",
    axis.line = element_blank(),
    axis.text = element_text(size = BASEPTFONTSIZE - 2, color = "#000000"),
    axis.title.y = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    axis.title.x = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    strip.text = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    strip.background = element_rect(color = "#ffffff00", fill = "#ffffff00")
  )


################################################################################
# cell cycle
################################################################################
Idents(seuMerged_art) <- "hivPos"
library(tricycle)

seuMerged_art <- SeuratWrappers::Runtricycle(
  object = seuMerged_art,
  slot = "data",
  reduction.name = "tricycleEmbedding", 
  reduction.key = "tricycleEmbedding_",
  gname = NULL,
  gname.type = "SYMBOL",
  species = "human")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

seuMerged_art <- CellCycleScoring(
  seuMerged_art,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = FALSE)


supDataCheckpoint <- data.frame(
  cellName = seuMerged_art$cellName,
  bin = seuMerged_art$manualAnnotBin,
  serpinb9_UCell = seuMerged_art$serpinb9_UCell,
  bcl2_UCell = seuMerged_art$bcl2_UCell,
  caspase_UCell = seuMerged_art$caspase_UCell,
  tricyclePos = seuMerged_art$tricyclePosition,
  S.score = seuMerged_art$S.Score,
  G2M.score = seuMerged_art$G2M.Score
)

saveRDS(supDataCheckpoint, "rds/supFig_art_data.rds")


ccpos_plot <- tricycle::plot_ccposition_den(
  seuMerged_art$tricyclePosition, seuMerged_art$hivPos, 'HIV', bw = 30)


ccpos_plot_tmp <- tricycle::plot_ccposition_den(
  seuMerged_art$tricyclePosition, seuMerged_art$manualAnnot, 'HIV', bw = 30, palette.v = multiHueColorPalette)

ccpos_plot_tmp + 
  guides(
    color = guide_legend(override.aes = list(linewidth = 3))) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10, color = "#000000", hjust = 0),
    legend.key.width = unit(25, "pt"),
    legend.margin = margin(t = -10),
    legend.spacing.x = unit(BASEPTFONTSIZE, "pt"),
    axis.line = element_line(linewidth = BASELINEWIDTHSIZE, color = "#000000"),
    plot.title = element_blank(),
    axis.text.x = element_text(size = BASEPTFONTSIZE - 2, color = "#000000"),
    axis.text.y = element_text(size = BASEPTFONTSIZE - 2, color = "#000000"),
    axis.title.x = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    axis.title.y = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    legend.title = element_blank())

ccpos_plot_final <- ccpos_plot +
  theme_classic() +
  geom_line(aes(x = 0, y = 0, linetype = "All cells", color = "All cells"),
    show.legend = TRUE, inherit.aes = FALSE) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 2))) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = BASEPTFONTSIZE, color = "#000000", hjust = 0),
    legend.margin = margin(t = -10),
    legend.spacing.x = unit(BASEPTFONTSIZE, "pt"),
    axis.line = element_line(linewidth = BASELINEWIDTHSIZE, color = "#000000"),
    plot.title = element_blank(),
    axis.text.x = element_text(size = BASEPTFONTSIZE - 2, color = "#000000"),
    axis.text.y = element_text(size = BASEPTFONTSIZE - 2, color = "#000000"),
    axis.title.x = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    axis.title.y = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    legend.title = element_blank()) +
  scale_linetype_manual(values = c("dashed")) +
  scale_color_manual(values = c("HIV-" = "#807c7c", "HIV+" = "#ff0000", "blue"))



################################################################################
# plot compilation
################################################################################
layout <- c(
  area(1, 1, 10, 4), # art rna heatmap
  area(1, 5, 3, 8), # prev reported
  area(4, 5, 6, 6),
  area(4, 7, 6, 8),
  area(7, 5, 8, 6),
  area(7, 7, 8, 8),
  area(9, 5, 10, 6),
  area(9, 7, 10, 8)
)

set.seed(42)
p <- wrap_elements(
  full = grid.grabExpr(draw(figA,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom",
    merge_legend = TRUE,
    padding = unit(c(5, 5, 5, 5), "points")))
  ) +
  free(figB_mod + theme(
    plot.tag = element_text(margin = margin(l= -15, t = 5)),
    plot.margin = margin(l = 0, r = 20, b = -40)), type = "panel", side = "r") +
  free(all_deg_plot + theme(
    plot.tag = element_text(margin = margin(l = -15, t = 20)),
    plot.margin = margin(l = 10, r = 35)), type = "space", side = "tlr") +
  free(all_dea_plot + theme(
    plot.tag = element_text(margin = margin(l = -45, t = 20)),
    plot.margin = margin(l = 15, r = 20)), type = "space", side = "tlr") +
  free(tem_deg_plot + theme(
    plot.tag = element_text(margin = margin(l = -15, t = 20)),
    plot.margin = margin(l = 10, r = 35, t = -18)), type = "space", side = "tlr") +
  free(tem_dea_plot + theme(
    plot.tag = element_text(margin = margin(l = -45, t = 20)),
    plot.margin = margin(l = 15, r = 20, t = -18)), type = "space", side = "tlr") +
  free(tcm_deg_plot + theme(
    plot.tag = element_text(margin = margin(l = -15, t = 20)),
    plot.margin = margin(l = 10, r = 35, t = -18)), type = "space", side = "tlr") +
  free(tcm_dea_plot + theme(
    plot.tag = element_text(margin = margin(l = -45, t = 20)),
    plot.margin = margin(l = 15, r = 20, t = -18)), type = "space", side = "tlr") +
  plot_annotation(tag_levels = "A") +
  plot_layout(design = layout) &
  plotTagTheme

saveFinalFigure(
  plot = p,
  fn = "fig4_v4",
  prefixDir = "figures/outs",
  devices = c("png", "pdf"),
  addTimestamp = TRUE,
  units = "cm",
  gwidth = 7 * 2.54,
  gheight = JCI_MAX_HEIGHT)


sup_layout <- c(
  area(1, 1, 6, 8), # cell death genes of interest
  area(7, 1, 11, 8), # ucell plots
  area(12, 1, 14, 8) # tricycle
)

sup_p <-
  supFigA_mod +
  supFigUcellPlot +
  ccpos_plot_final +
  plot_annotation(tag_levels = "A") +
  plot_layout(design = sup_layout) &
  plotTagTheme

saveFinalFigure(
  plot = sup_p,
  fn = "supfig5_v3",
  prefixDir = "figures/outs",
  devices = c("png", "pdf"),
  addTimestamp = TRUE,
  gwidth = 7,
  gheight = 10)

