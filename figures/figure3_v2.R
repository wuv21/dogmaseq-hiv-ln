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
# epigenetic analysis
################################################################################
heatmapCutoff <- "FDR < 0.05 & MeanDiff > 0.5"
chronic_markerTest_deviations <- getMarkerFeatures(
  ArchRProj = projChronicSub2, 
  useMatrix = "MotifMatrix",
  groupBy = "manualAnnotBinWithHivStatus",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  verbose = TRUE,
  useSeqnames = "z",
)

chronic_markerTest_deviations_heatmapText <- bind_rows(
  lapply(getMarkers(chronic_markerTest_deviations, cutOff = heatmapCutoff), function(x) {
    return(as.data.frame(x))
  }
  ), .id = "group")

chronic_markerTest_deviations_heatmap <- plotMarkerHeatmap(
  seMarker = chronic_markerTest_deviations,
  log2Norm = FALSE,
  cutOff = heatmapCutoff,
  transpose = FALSE,
  returnMatrix = TRUE,
  scaleRows = FALSE
)

chronic_markerTest_deviations_heatmapSig <- matrix(" ",
  nrow = nrow(chronic_markerTest_deviations_heatmap),
  ncol = ncol(chronic_markerTest_deviations_heatmap))

dimnames(chronic_markerTest_deviations_heatmapSig) <- dimnames(chronic_markerTest_deviations_heatmap)

for (i in seq_along(chronic_markerTest_deviations_heatmapText$group)) {
  row <- chronic_markerTest_deviations_heatmapText[i, ]

  chronic_markerTest_deviations_heatmapSig[row$name, row$group] <- case_when(
    row$FDR < 0.01 ~ "✱✱",
    row$FDR < 0.05 ~ "✱",
    TRUE ~ " ")
}

set.seed(42)
tmp <- ComplexHeatmap::Heatmap(
  matrix = chronic_markerTest_deviations_heatmap,
  name = " ",
  border = FALSE,
  cluster_rows = TRUE,
  row_dend_side = "left",
  show_row_dend = TRUE,
  cluster_columns = FALSE,
  column_split = rep(c(1:4), each = 3),
  show_column_names = FALSE,
  show_heatmap_legend = TRUE,
  column_dend_side = "top",
  column_title = NULL,
  row_title_gp = gpar(fontsize = 6),
  layer_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    grid.text(pindex(chronic_markerTest_deviations_heatmapSig, i, j), x, y,
      gp = gpar(fontsize = 2, hjust = 0.5, vjust = 0.5))
  },
  heatmap_legend_param = list(
    legend_direction = "horizontal",
    legend_width = unit(2, "cm"),
    legend_height = unit(0.3, "cm"),
    labels_gp = gpar(fontsize = 8)),
  show_column_dend = TRUE,
  row_km = 12,
  row_km_repeats = 100,
  # row_labels = as.character(c(1:194)),
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  top_annotation = ComplexHeatmap::HeatmapAnnotation(
    `Phenotype` = anno_block(
      gp = gpar(fill = "#ffffff", col = "#ffffff"),
      height = unit(BASEPTFONTSIZE, "pt"),
      labels = unique(str_split_fixed(colnames(chronic_markerTest_deviations_heatmap), "_", n = 2)[, 1]),
      labels_gp = gpar(col = "#000000", fontsize = BASEPTFONTSIZE - 2),
      show_name = FALSE),
    simple_anno_size = unit(0.3, "cm")
  ),
  bottom_annotation = ComplexHeatmap::HeatmapAnnotation(
    `Phenotype` = anno_block(
      gp = gpar(fill = "#ffffff", col = "#ffffff"),
      height = unit(BASEPTFONTSIZE, "pt"),
      labels = unique(str_split_fixed(colnames(chronic_markerTest_deviations_heatmap), "_", n = 2)[, 1]),
      labels_gp = gpar(col = "#000000", fontsize = BASEPTFONTSIZE - 2),
      show_name = TRUE),
    # `Phenotype` = str_split_fixed(colnames(chronic_markerTest_deviations_heatmap), "_", n = 2)[, 1],
    `HIV` = str_split_fixed(colnames(chronic_markerTest_deviations_heatmap), "_", n = 2)[, 2],
    col = list(
      # `Phenotype` = c("Tcm/Th2" = "#a6611a", "Tem/Th1" = "#dfc27d", "Tfh" = "#80cdc1", "Treg" = "#018571"),
      `HIV` = COLORS$hivStatus
    ),
    annotation_legend_param = list(
      labels_gp = gpar(fontsize = 8),
      title_gp = gpar(fontsize = 8)
    ),
    annotation_name_gp = gpar(fontsize = 8),
    simple_anno_size = unit(0.3, "cm")
  )
  # right_annotation = rowAnnotation(
  #   # ind = anno_text(
  #   #   x = c(1:nrow(chronic_markerTest_deviations_heatmap)),
  #   #   rot = 0,
  #   #   gp = gpar(fontsize = 5, fontfamily = "mono")),
  #   foo = anno_block(gp = gpar(fill = "#000000"), width = unit(12/14, "pt")),
  #   bar = anno_block(
  #     graphics = function(index, levels) {
  #       lbls <- rownames(chronic_markerTest_deviations_heatmap)[index]
  #       lbls <- gsub("_\\d+", "", lbls)
  #       lbls <- str_wrap(paste(lbls, collapse = " "), width = 40)
  # 
  #       grid.rect(gp = gpar(fill = NA, col = NA))
  #       txt = paste(lbls, collapse = ",")
  #       grid.text(txt, 0.01, 0.5, rot = 0, hjust = 0, gp = gpar(fontsize = 6))
  #     },
  #     width = unit(6, "cm"))
  # )
)

set.seed(42)
heatmap_tmp <- draw(tmp)
final_row_order <- unname(unlist(row_order(heatmap_tmp)))
new_row_numbers <- rep(0, nrow(chronic_markerTest_deviations_heatmap))
new_row_numbers[final_row_order] <- c(1:length(new_row_numbers))


rightAnnot <- rowAnnotation(
  ind = anno_text(
    x = new_row_numbers,
    rot = 0,
    just = "left",
    gp = gpar(fontsize = 3, fontfamily = "mono"), width = unit(0.4, "cm")),
  foo = anno_block(gp = gpar(fill = "#000000"), width = unit(12/14, "pt")),
  bar = anno_block(
    graphics = function(index, levels) {
      lbls <- rownames(chronic_markerTest_deviations_heatmap)[index]
      lbls <- gsub("_\\d+", "", lbls)
      nums <- new_row_numbers[index]
      
      lbls <- str_wrap(paste(paste0(nums, ": ", lbls), collapse = ", "), width = 40)
      
      print("lbls")
      print(lbls)
      
      grid.rect(gp = gpar(fill = NA, col = NA))
      grid.text(lbls, 0.01, 0.5, rot = 0, hjust = 0, gp = gpar(fontsize = 5.2))
    },
    width = unit(6, "cm"))
)
# set.seed(42)
# draw(tmp + rightAnnot)


projChronicSub2$manualAnnotBinWithHivShortStatus <- paste0(projChronicSub2$manualAnnotBin, "_", projChronicSub2$hivPos)

tcmTh2_ded <- getMarkerFeatures(
  ArchRProj = projChronicSub2,
  useMatrix = "MotifMatrix",
  groupBy = "manualAnnotBinWithHivShortStatus",
  useGroups = "Tcm/Th2_HIV+",
  bgdGroups = "Tcm/Th2_HIV-",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  verbose = TRUE,
  useSeqnames = "z",
)

tcmTh2_ded_filt <- getMarkers(
  tcmTh2_ded,
  cutOff = "FDR < 0.4"
)

tem_ded <- getMarkerFeatures(
  ArchRProj = projChronicSub2,
  useMatrix = "MotifMatrix",
  groupBy = "manualAnnotBinWithHivShortStatus",
  useGroups = "Tem/Th1_HIV+",
  bgdGroups = "Tem/Th1_HIV-",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  verbose = TRUE,
  useSeqnames = "z",
)

tem_ded_filt <- getMarkers(
  tem_ded,
  cutOff = "FDR < 0.4"
)


tfh_ded <- getMarkerFeatures(
  ArchRProj = projChronicSub2,
  useMatrix = "MotifMatrix",
  groupBy = "manualAnnotBinWithHivShortStatus",
  useGroups = "Tfh_HIV+",
  bgdGroups = "Tfh_HIV-",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  verbose = TRUE,
  useSeqnames = "z",
)

tfh_ded_filt <- getMarkers(
  tfh_ded,
  cutOff = "FDR < 0.4"
)

treg_ded <- getMarkerFeatures(
  ArchRProj = projChronicSub2,
  useMatrix = "MotifMatrix",
  groupBy = "manualAnnotBinWithHivShortStatus",
  useGroups = "Treg_HIV+",
  bgdGroups = "Treg_HIV-",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  verbose = TRUE,
  useSeqnames = "z",
)

treg_ded_filt <- getMarkers(
  treg_ded,
  cutOff = "FDR < 0.4"
)


ded_list <- list(
  "Tcm/Th2" = bind_rows(lapply(tcmTh2_ded_filt@listData, as.data.frame), .id = "condition"),
  "Tem/Th1" = bind_rows(lapply(tem_ded_filt@listData, as.data.frame), .id = "condition"),
  "Tfh" = bind_rows(lapply(tfh_ded_filt@listData, as.data.frame), .id = "condition"),
  "Treg" = bind_rows(lapply(treg_ded_filt@listData, as.data.frame), .id = "condition")
)
  
ded_df <- dplyr::bind_rows(ded_list, .id = "celltype") %>%
  dplyr::mutate(direction = ifelse(MeanDiff > 0, "Enriched in HIV+", "Enriched in HIV-")) %>%
  dplyr::group_by(celltype, direction) %>%
  dplyr::slice_max(abs(MeanDiff), n = 10) %>%
  dplyr::arrange(desc(MeanDiff)) %>%
  dplyr::mutate(significantFill = case_when(
    MeanDiff > 0 & FDR < 0.05 ~ "#ff0000ff",
    MeanDiff > 0 & FDR >= 0.05 ~ "#ff000000",
    MeanDiff < 0 & FDR < 0.05 ~ "#807c7cff",
    MeanDiff < 0 & FDR >= 0.05 ~ "#807c7c00",
  )) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(name = stringr::str_match(name, "(.*)_(\\d+)")[, 2]) %>%
  dplyr::mutate(xaxisName = paste0(name, "_", row_number())) %>%
  dplyr::mutate(xaxisName = factor(xaxisName, levels = xaxisName))

ded_plots <- ded_df %>%
  ggplot(aes(y = MeanDiff, x = xaxisName, color = direction, fill = significantFill)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#cfcfcf50") +
  geom_point(size = 0.75, pch = 21) +
  facet_wrap(~ celltype, ncol = 1, scales = "free") +
  theme_classic() +
  labs(x = "chromVAR motif",
    y = "<span style='font-size:7pt'>Mean difference in chromVAR deviations<br><span style='color:#ff0000'>Enriched in HIV+</span> | <span style='color:#807c7c'>Enriched in HIV-</span></span>") +
  scale_x_discrete(breaks = ded_df$xaxisName, labels = ded_df$name) +
  scale_color_manual(values = c("Enriched in HIV+" = "#ff0000", "Enriched in HIV-" = "#807c7c")) +
  scale_fill_identity() +
  coord_cartesian(clip = "off") +
  theme(
    panel.spacing.y = unit(2.5, "pt"),
    plot.margin = margin(l = 10),
    axis.line = element_line(linewidth = BASELINEWIDTHSIZE),
    strip.text = element_text(size = BASEPTFONTSIZE, color = "#000000", margin = margin(0,0,0,0)),
    strip.background = element_rect(fill = "#00000000", color = "#00000000"),
    text = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    axis.title.y = element_markdown(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "#000000", size = BASEPTFONTSIZE - 2),
    axis.text.y = element_text(color = "#000000", size = BASEPTFONTSIZE - 2),
    legend.position = "blank"
  )


################################################################################
# dea
################################################################################
seuMerged$manualAnnotBinWithHivShortStatus <- paste0(seuMerged$manualAnnotBin, "_", seuMerged$hivPos)
Idents(seuMerged) <- "manualAnnotBinWithHivShortStatus"

seuMergedChronic <- subset(seuMerged, subset = stage == "Viremic")

Idents(seuMergedChronic) <- "hivPos"
all_dea <- FindMarkers(
  seuMergedChronic,
  features = tsaCatalog[!tsaCatalog$isCtrl, "DNA_ID"],
  ident.1 = "HIV+",
  ident.2 = "HIV-",
  min.pct = 0.3,
  assay = "ADT"
)

all_deg <- FindMarkers(
  seuMergedChronic,
  ident.1 = "HIV+",
  ident.2 = "HIV-",
  min.pct = 0.3,
  assay = "RNA"
)

deLollipop_theme <- list(
  theme_classic(),
  scale_color_manual(values = c("TRUE" = "#ff0000", "FALSE" = "#807c7c")),
  labs(x = "<span style='font-size:7pt'>Avg log2(FC)<br><span style='color:#807c7c'>&uarr; in HIV-</span> | <span style='color:#ff0000'>&uarr; in HIV+</span></span>"),
  theme(
    axis.line = element_line(linewidth = BASELINEWIDTHSIZE),
    legend.position = "blank",
    axis.text = element_text(size = BASEPTFONTSIZE - 2, color = "#000000"),
    axis.title.y = element_blank(),
    axis.title.x = element_markdown(),
    plot.title = element_text(size = BASEPTFONTSIZE, color = "#000000", hjust = 0.5)
  ))


all_dea_plot <- all_dea %>%
  as.data.frame() %>%
  dplyr::mutate(DNA_ID = rownames(.)) %>%
  left_join(tsaCatalog[, c("DNA_ID", "cleanName")]) %>%
  dplyr::mutate(p_val_adj2 = p.adjust(p_val, method = "BH")) %>%
  dplyr::filter(p_val_adj2 < 0.05) %>%
  dplyr::mutate(direction = avg_log2FC > 0) %>%
  dplyr::group_by(direction) %>%
  dplyr::slice_max(abs(avg_log2FC), n = 10) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::mutate(cleanName = factor(cleanName, levels = cleanName)) %>%
  ggplot(aes(y = cleanName, x = avg_log2FC, color = direction)) +
  geom_vline(xintercept = 0) +
  geom_segment(aes(x = 0, xend = avg_log2FC), linetype = "dotted") +
  geom_point(size = 0.75) +
  deLollipop_theme +
  labs(y = "Surface protein", title = "All cells\nSurface protein")

all_deg_plot <- all_deg %>%
  as.data.frame() %>%
  dplyr::mutate(gene = rownames(.)) %>%
  dplyr::filter(gene != "HIV") %>%
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
  geom_point(size = 0.75) +
  deLollipop_theme +
  labs(y = "Gene", title = "All cells\nRNA")


################################################################################
# cell cycle analysis
################################################################################
library(tricycle)

seuMergedChronic <- SeuratWrappers::Runtricycle(
  object = seuMergedChronic,
  slot = "data",
  reduction.name = "tricycleEmbedding", 
  reduction.key = "tricycleEmbedding_",
  gname = NULL,
  gname.type = "SYMBOL",
  species = "human")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

seuMergedChronic <- CellCycleScoring(seuMergedChronic,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = FALSE)


ccpos_plot <- tricycle::plot_ccposition_den(
  seuMergedChronic$tricyclePosition, seuMergedChronic$hivStatus, 'HIV status', bw = 30)

ccpos_plot_final <- ccpos_plot +
  theme_classic() +
  geom_line(aes(x = 0, y = 0, linetype = "All cells", color = "All cells"),
    show.legend = TRUE, inherit.aes = FALSE) +
  guides(color = guide_legend(override.aes = list(linewidth = 2))) +
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
  scale_color_manual(labels = c("HIV DNA+ RNA-", "HIV DNA+/- RNA+", "HIV-", "All cells"), values = c(COLORS$hivStatus, "#ffffff"))


ccpos_tests <- c(
  ks.test(
    x = seuMergedChronic$tricyclePosition[seuMergedChronic$hivStatus == "HIV DNA+/- RNA+"],
    y = seuMergedChronic$tricyclePosition[seuMergedChronic$hivStatus == "HIV DNA+ RNA-"])$p,
  
  
  ks.test(
    x = seuMergedChronic$tricyclePosition[seuMergedChronic$hivStatus == "HIV DNA+/- RNA+"],
    y = seuMergedChronic$tricyclePosition[seuMergedChronic$hivStatus == "HIV-"])$p,
  
  ks.test(
    x = seuMergedChronic$tricyclePosition[seuMergedChronic$hivStatus == "HIV DNA+ RNA-"],
    y = seuMergedChronic$tricyclePosition[seuMergedChronic$hivStatus == "HIV-"])$p)


ks.test(
  x = seuMergedChronic$tricyclePosition[seuMergedChronic$hivPos == "HIV+" & seuMergedChronic$manualAnnotBin == "Tfh"],
  y = seuMergedChronic$tricyclePosition[seuMergedChronic$hivPos == "HIV-" & seuMergedChronic$manualAnnotBin == "Tfh"])

ks.test(
  x = seuMergedChronic$tricyclePosition[seuMergedChronic$hivPos == "HIV+" & seuMergedChronic$manualAnnotBin == "Tfh"],
  y = seuMergedChronic$tricyclePosition[seuMergedChronic$hivPos == "HIV-" & seuMergedChronic$manualAnnotBin == "Tfh"])

################################################################################
# supplemental figure
################################################################################
Idents(seuMergedChronic) <- "manualAnnotBinWithHivShortStatus"
tfh_dea <- FindMarkers(
  seuMergedChronic,
  features = tsaCatalog[!tsaCatalog$isCtrl, "DNA_ID"],
  ident.1 = "Tfh_HIV+",
  ident.2 = "Tfh_HIV-",
  min.pct = 0.3,
  assay = "ADT")

tfh_deg <- FindMarkers(
  seuMergedChronic,
  ident.1 = "Tfh_HIV+",
  ident.2 = "Tfh_HIV-",
  assay = "RNA",
  min.pct = 0.3)

tfh_dea_plot <- tfh_dea %>%
  as.data.frame() %>%
  dplyr::mutate(DNA_ID = rownames(.)) %>%
  left_join(tsaCatalog[, c("DNA_ID", "cleanName")]) %>%
  dplyr::mutate(p_val_adj2 = p.adjust(p_val, method = "BH")) %>%
  dplyr::filter(p_val_adj2 < 0.05) %>%  dplyr::mutate(direction = avg_log2FC > 0) %>%
  dplyr::group_by(direction) %>%
  dplyr::slice_max(abs(avg_log2FC), n = 10) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::mutate(cleanName = factor(cleanName, levels = cleanName)) %>%
  ggplot(aes(y = cleanName, x = avg_log2FC, color = direction)) +
  geom_vline(xintercept = 0) +
  geom_segment(aes(x = 0, xend = avg_log2FC), linetype = "dotted") +
  geom_point(size = 0.75) +
  deLollipop_theme +
  labs(y = "Surface protein", title = "Tfh cells\nSurface protein")

tfh_deg_plot <- tfh_deg %>%
  as.data.frame() %>%
  dplyr::mutate(gene = rownames(.)) %>%
  dplyr::filter(gene != "HIV") %>%
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
  geom_point(size = 0.75) +
  deLollipop_theme +
  labs(y = "Gene", title = "Tfh cells\nRNA")



ignoredBins <- c("Other", "Naive/early memory")

hivNegColors <- c(
  "#c6dbef",
  "#9ecae1",
  "#00ff00",
  "#4292c6",
  "#2171b5",
  "#084594")

hivPosColors <- c(
  "#fcbba1",
  "#fc9272",
  "#9D00FF",
  "#ef3b2c",
  "#cb181d",
  "#99000d")

altColors <- c(rbind(hivNegColors, hivPosColors))

ccpos_plot_byBin <- tricycle::plot_ccposition_den(
  seuMergedChronic$tricyclePosition[seuMergedChronic$manualAnnotBin %ni% ignoredBins],
  seuMergedChronic$manualAnnotBinWithHivShortStatus[seuMergedChronic$manualAnnotBin %ni% ignoredBins],
  palette.v = altColors,
  'HIV status',
  bw = 10)

ccpos_plot_byBin_final <- ccpos_plot_byBin +
  geom_line(aes(x = 0, y = 0, linetype = "All cells", color = "All cells"),
    show.legend = TRUE, inherit.aes = FALSE) +
  theme_classic() +
  guides(color = guide_legend(override.aes = list(linewidth = 2), nrow = 2)) +
  scale_linetype_manual(values = c("dashed")) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = BASEPTFONTSIZE - 2, color = "#000000", hjust = 0),
    legend.key.spacing.x = unit(BASEPTFONTSIZE - 6, "pt"),
    legend.key.width = unit(6, "pt"),
    axis.line = element_line(linewidth = BASELINEWIDTHSIZE, color = "#000000"),
    plot.title = element_blank(),
    axis.text.x = element_text(size = BASEPTFONTSIZE - 2, color = "#000000"),
    axis.text.y = element_text(size = BASEPTFONTSIZE - 2, color = "#000000"),
    axis.title.x = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    axis.title.y = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    legend.title = element_blank())



################################################################################
# plot compilation
################################################################################
layout <- c(
  area(1, 1, 9, 8), # chronic heatmap
  area(1, 9, 4, 12), # ded plots
  area(5, 9, 6, 10), # deg/dea plots
  area(5, 11, 6, 12), # deg/dea plots
  area(7, 9, 8, 12) # tricycle
)

set.seed(42)
p <- wrap_elements(
  grid.grabExpr(draw(tmp + rightAnnot,
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom",
    merge_legend = TRUE))
  ) +
  free(ded_plots + theme(plot.tag = element_text(margin = margin(l = -65))), type = "space", side = "lr") +
  
  free(all_deg_plot + theme(
    plot.tag = element_text(margin = margin(l = -65, t = 10)),
    plot.margin = margin(l = 10, r = 35)), type = "space", side = "lr") +
  free(all_dea_plot + theme(
    plot.tag = element_text(margin = margin(l = -65, t = 10)),
    plot.margin = margin(l = 15, r = 20)), type = "space", side = "lr") +
  
  
  free(ccpos_plot_final + theme(
    plot.tag = element_text(margin = margin(l = -65, t = 10)),
    plot.margin = margin(l = 10, r = 35, b = -10)), type = "space", side = "blr") +
  
  plot_annotation(tag_levels = "A") +
  plot_layout(design = layout) &
  plotTagTheme

saveFinalFigure(
  plot = p,
  fn = "fig3_v2",
  prefixDir = "figures/outs",
  devices = c("png", "pdf"),
  addTimestamp = TRUE,
  gwidth = 8.5,
  gheight = 9)


sup_layout <- c(
  area(1, 1, 2, 2), # deg/dea plots
  area(1, 3, 2, 4), # deg/dea plots
  area(3, 1, 5, 4) # tricycle
)

sup_p <-
  free(tfh_deg_plot + theme(
    plot.tag = element_text(margin = margin(l = -35, t = 10)),
    plot.margin = margin(l = 30, r = 35)), type = "space", side = "l") +
  
  free(tfh_dea_plot + theme(
    plot.tag = element_text(margin = margin(l = -45, t = 10)),
    plot.margin = margin(l = 15, r = 20)), type = "space", side = "l") +
  
  free(ccpos_plot_byBin_final + theme(
    plot.tag = element_text(margin = margin(l = -35, t = 10)),
    plot.margin = margin(l = 35, r = 35)), type = "space", side = "lr") +
  
  plot_annotation(tag_levels = "A") +
  plot_layout(design = sup_layout) &
  plotTagTheme

saveFinalFigure(
  plot = sup_p,
  fn = "supfig4_v1",
  prefixDir = "figures/outs",
  devices = c("png", "pdf"),
  addTimestamp = TRUE,
  gwidth = 5,
  gheight = 6)

