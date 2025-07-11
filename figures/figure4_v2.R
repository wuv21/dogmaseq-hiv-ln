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
projUnhashedClean$manualAnnot <- trimws(projUnhashedClean$manualAnnot)


uniqClusters <- unique(projUnhashedClean$manualAnnot)
uniqStages <- c("Viremic", "ART")
uniqInfection <- names(COLORS$hivStatus)

################################################################################
# epigenetic analysis
################################################################################
heatmapCutoff <- "FDR < 0.05 & MeanDiff > 0.5"
# art temp
art_markerTest_deviations <- getMarkerFeatures(
  ArchRProj = projARTsub, 
  useMatrix = "MotifMatrix",
  groupBy = "manualAnnotBinWithHiv",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  verbose = TRUE,
  useSeqnames = "z",
)


# plots the mean z score for each motif...
art_markerTest_deviations_heatmap <- plotMarkerHeatmap(
  seMarker = art_markerTest_deviations, 
  log2Norm = FALSE,
  cutOff = "FDR < 0.05 & MeanDiff > 0.5",
  transpose = FALSE,
  returnMatrix = TRUE,
  scaleRows = FALSE
)

art_markerTest_deviations_heatmapText <- bind_rows(
  lapply(getMarkers(art_markerTest_deviations, cutOff = "FDR < 0.05 & MeanDiff > 0.5"), function(x) {
    return(as.data.frame(x))
  }
  ), .id = "group")

art_markerTest_deviations_heatmapSig <- matrix(" ",
  nrow = nrow(art_markerTest_deviations_heatmap),
  ncol = ncol(art_markerTest_deviations_heatmap))

dimnames(art_markerTest_deviations_heatmapSig) <- dimnames(art_markerTest_deviations_heatmap)

for (i in seq_along(art_markerTest_deviations_heatmapText$group)) {
  row <- art_markerTest_deviations_heatmapText[i, ]
  
  art_markerTest_deviations_heatmapSig[row$name, row$group] <- case_when(
    row$FDR < 0.01 ~ "✱✱",
    row$FDR < 0.05 ~ "✱",
    TRUE ~ " ")
}

tmp <- ComplexHeatmap::Heatmap(
  matrix = art_markerTest_deviations_heatmap,
  name = " ",
  border = FALSE,
  cluster_rows = TRUE,
  show_row_names = FALSE,
  row_dend_side = "left",
  show_row_dend = TRUE,
  cluster_columns = FALSE,
  column_split = rep(c(1:7), each = 2),
  show_column_names = FALSE,
  show_heatmap_legend = TRUE,
  column_dend_side = "top",
  column_title = NULL,
  layer_fun = function(j, i, x, y, w, h, col) { # add text to each grid
    grid.text(pindex(art_markerTest_deviations_heatmapSig, i, j), x, y,
      gp = gpar(fontsize = 2, hjust = 0.5, vjust = 0.5))
  },
  heatmap_legend_param = list(
    legend_direction = "horizontal",
    legend_width = unit(3, "cm")),
  show_column_dend = TRUE,
  # column_km = 3,
  row_km = 10,
  bottom_annotation = ComplexHeatmap::HeatmapAnnotation(
    `Phenotype` = str_split_fixed(colnames(art_markerTest_deviations_heatmap), "_", n = 2)[, 1],
    `HIV` = ifelse(str_split_fixed(colnames(art_markerTest_deviations_heatmap), "_", n = 2)[, 2] == "HIV-", "HIV-", "HIV+"),
    col = list(
      # `Phenotype` = c("Tem/effector" = "#1b9e77", "MAIT" = "#d95f02", "Tcm/Ttm" = "#7570b3"),
      `HIV` = c("HIV+" = HIVPOSCOLOR, "HIV-" = HIVNEGCOLOR)
    ),
    annotation_name_gp = gpar(fontsize = 6)
  ),
  right_annotation = rowAnnotation(
    foo = anno_block(gp = gpar(fill = "#000000"), width = unit(1, "mm")),
    bar = anno_block(
      graphics = function(index, levels) {
        lbls <- rownames(art_markerTest_deviations_heatmap)[index]
        lbls <- gsub("_\\d+", "", lbls)
        lbls <- str_wrap(paste(lbls, collapse = " "), width = 40)
        
        grid.rect(gp = gpar(fill = NA, col = NA))
        txt = paste(lbls, collapse = ",")
        grid.text(txt, 0.01, 0.5, rot = 0, hjust = 0, gp = gpar(fontsize = 8))
      },
      width = unit(7, "cm"))
  )
)

savePlot(grid.grabExpr(draw(tmp,
  heatmap_legend_side = "bottom",
  annotation_legend_side = "bottom",
  merge_legend = TRUE)), fn = "art_deviations", gheight = 11, gwidth = 8, devices = "png")

rm(tmp)


seuMerged_art <- subset(seuMerged, subset = stage == "ART")
Idents(seuMerged_art) <- "hivPos"
FindMarkers(
  seuMerged_art,
  features = c("NR4A1"),
  ident.1 = "HIV+",
  ident.2 = "HIV-",
  min.pct = 0.3,
  assay = "RNA"
)


################################################################################
# plot compilation
################################################################################
layout <- c(
  area(1, 1, 9, 8), # chronic heatmap
  area(1, 9, 4, 12), # ded plots
  area(5, 9, 6, 10), # deg/dea plots
  area(5, 11, 6, 12), # deg/dea plots
  area(7, 9, 8, 10), # deg/dea plots
  area(7, 11, 8, 12) # deg/dea plots
)

p <- wrap_elements(
  grid.grabExpr(draw(tmp,
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
  
  free(tfh_deg_plot + theme(
    plot.tag = element_text(margin = margin(l = -65, t = 10)),
    plot.margin = margin(l = 10, r = 35)), type = "space", side = "lr") +
  free(tfh_dea_plot + theme(
    plot.tag = element_text(margin = margin(l = -65, t = 10)),
    plot.margin = margin(l = 15, r = 20)), type = "space", side = "lr") +
  plot_annotation(tag_levels = "A") +
  plot_layout(design = layout) &
  plotTagTheme

saveFinalFigure(plot = p, fn = "fig4_v1", prefixDir = "figures/outs", devices = c("png"), gwidth = 8.5, gheight = 9)

