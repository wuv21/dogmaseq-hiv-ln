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
  stat_summary(geom = "bar", width = 0.6, fun = mean, fill = "#ffffff00", aes(color = infection)) +
  geom_vline(xintercept = -Inf, linewidth = BASELINEWIDTHSIZE, color = "#000000") +
  geom_hline(yintercept = 0, linewidth = BASELINEWIDTHSIZE, color = "#000000") +
  geom_jitter(aes(fill = infection), color = "#000000", pch = 21, height = 0, alpha = 0.75, size = 0.8, show.legend = FALSE, width = 0.2) +
  stat_summary(fun = mean_se_onlyPos, geom = "crossbar", width = 0.3, linewidth = 0.3) +
  stat_summary(fun.data = mean_se_onlyPosVertical, geom = "linerange", linewidth = 0.3) +
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
    axis.text.x = element_text(color = "#000000", size = 6, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(color = "#000000", size = 6),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(color = "#eeeeee", linetype = "dotted"),
    strip.text = element_text(color = "#000000", size = 6),
    axis.title = element_text(color = "#000000", size = 6),
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


depletionCompiledStats <- ggpubr::compare_means(percent ~ combinedFactor,
  depletionDf2, paired = FALSE)

depletionCompiledStats$p.sym <- pValSymnum(depletionCompiledStats$p.adj)
depletionCompiledStats <- depletionCompiledStats %>%
  dplyr::filter(p.sym != "ns") %>%
  dplyr::mutate(p.position = 105 - 4.75 * (row_number() %% 7)) %>%
  dplyr::select(-method, -p.signif, -p.format)

depletionDfG <- depletionDf2 %>%
  ggplot(aes(y = percent, x = combinedFactor)) +
  geom_line(aes(group = donor), color = "#bbbbbb75", linewidth = 0.3) +
  stat_summary(geom = "bar", width = 0.6, fun = mean, fill = "#ffffff00", aes(color = stage)) +
  geom_vline(xintercept = -Inf, linewidth = BASELINEWIDTHSIZE, color = "#000000") +
  geom_hline(yintercept = 0, linewidth = BASELINEWIDTHSIZE, color = "#000000") +
  geom_point(color = "#000000", pch = 21, alpha = 0.75, size = 0.8, show.legend = FALSE) +
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
  stat_summary(fun = mean_se_onlyPos, geom = "crossbar", width = 0.3, linewidth = 0.3) +
  stat_summary(fun.data = mean_se_onlyPosVertical, geom = "linerange", linewidth = 0.3) +
  scale_color_manual(values = COLORS$stage) +
  coord_cartesian(clip = "off") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  theme_classic() +
  theme(
    plot.margin = margin(r = -1, l = -5),
    legend.position = "blank",
    axis.line = element_blank(),
    axis.text = element_text(color = "#000000"),
    axis.text.x = element_text(color = "#000000", size = 6),
    axis.text.y = element_text(color = "#000000", size = 6),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(color = "#eeeeee", linetype = "dotted"),
    strip.text = element_text(color = "#000000", size = 6),
    axis.title = element_text(color = "#000000", size = 6),
    axis.title.x = element_blank(),
    strip.background = element_rect(color = "#ffffff00", fill = NA)) +
  labs(y = "% of donor's Tfh cells")


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
  features = tsaCatalog[!tsaCatalog$isCtrl, "DNA_ID"],
  assay = "ADT")

dea_blimp1 <- dea_blimp1 %>%
  dplyr::mutate(DNA_ID = rownames(.)) %>%
  left_join(dplyr::select(tsaCatalog, DNA_ID, cleanName)) %>%
  dplyr::select(-p_val) %>%
  dplyr::arrange(desc(avg_log2FC))

deLollipopTheme <- theme(
  strip.clip = "off",
  strip.text = element_text(color = "#000000", size = 6),
  strip.background = element_rect(color = "#00000000", fill = "#00000000"),
  axis.text = element_text(color = "#000000"),
  axis.text.x = element_text(color = "#000000", size = 6),
  axis.text.y = element_text(color = "#000000", size = 5),
  axis.title.x = element_text(color = "#000000", size = 6),
  axis.title.y = element_blank(),
  legend.position = "blank",
  panel.grid = element_blank(),
  axis.line = element_line(color = "#000000", linewidth = BASELINEWIDTHSIZE),
  plot.background = element_rect(fill = "#ffffff00", color = "#ffffff00"),
  panel.background = element_rect(fill = "#ffffff00", color = "#ffffff00"))


dea_blimp1G <- dea_blimp1 %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::mutate(direction = ifelse(avg_log2FC > 0, "Up in\nPRDM1+ GC-Tfh", "Up in\nPRDM1- GC-TFh")) %>%
  dplyr::group_by(direction) %>%
  dplyr::slice_max(abs(avg_log2FC), n = 10) %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::mutate(cleanName = gsub(" \\(Neuropilin-1\\)", "", cleanName)) %>%
  dplyr::mutate(cleanName = gsub(" \\(DNAM-1\\)", "", cleanName)) %>%
  dplyr::mutate(cleanName = gsub(" \\(Ecto-5'-nucleotidase\\)", "", cleanName)) %>%
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
    plot.margin = margin(r = -10, l = -15, t = 0, b = 25)
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
  min.pct = 0.3,
  assay = "RNA")

deg_blimp1G <- deg_blimp1 %>%
  dplyr::mutate(gene = rownames(.)) %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::mutate(direction = ifelse(avg_log2FC > 0, "Up in\nPRDM1+ GC-Tfh", "Up in\nPRDM1- GC-TFh")) %>%
  dplyr::filter(abs(avg_log2FC) > 0.3) %>%
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
    plot.margin = margin(r = -80, l = -30, t = 0, b = 25)
  )


################################################################################
# organoid infection
################################################################################
library(magick)

# rep_flow <- image_read_svg("external_data/fig5_tonsil_organoid_representative_v1.svg", width = 419*3, height = 267*3)
rep_flow_png <- image_read("external_data/fig5_tonsil_organoid_representative_v1.png")

# Produces a ggplot2 plot with the svg image:
rep_flow_gg <- image_ggplot(rep_flow_png)

organoid_df <- read.csv("external_data/20250520_data.csv", check.names = FALSE)

organoid_df <- organoid_df %>%
  pivot_longer(cols = c(2:ncol(organoid_df)), names_to = "gate", values_to = "percent") %>%
  mutate(as_tibble(str_match(gate, "(p24[+-])(blimp1(lo|mid|hi))"), .name_repair = ~ c("gate", "p24", "blimp1", "ignore"))) %>%
  mutate(p24 = factor(p24, levels = c("p24-", "p24+"))) %>%
  arrange(p24) %>%
  mutate(gate = paste0(blimp1, " ", p24)) %>%
  mutate(gate = gsub("\\+", "<sup>+</sup>", gate)) %>%
  mutate(gate = gsub("\\-", "<sup>-</sup>", gate)) %>%
  mutate(gate = gsub("blimp1(lo|mid|hi)", "Blimp1<sup>\\1</sup>", gate)) %>%
  mutate(condition = factor(condition, levels = c("Neg", "NL4.3 MOI=0.01", "NL4.3 MOI=0.1", "SUMA MOI=0.01", "SUMA MOI=0.1"))) %>%
  mutate(gate = factor(gate, levels = unique(gate)))
  
organoid_g <- organoid_df %>%
  ggplot(aes(x = gate, y = percent)) +
  geom_col(aes(color = condition), position = position_dodge(width = 0.85), width = 0.65, fill = "#ffffff") +
  geom_point(aes(group = condition), position = position_dodge(width = 0.85), size = 0.5) +
  geom_hline(yintercept = -Inf, linewidth = BASELINEWIDTHSIZE) +
  geom_vline(xintercept = -Inf, linewidth = BASELINEWIDTHSIZE) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  guides(x = legendry::guide_axis_nested(
    key = " ",
    type = "bracket")) +
  theme_classic() +
  scale_color_manual(values = c("#0c0c0c", "#bdd7e7", "#3182bd", "#fcae91", "#de2d26")) +
  labs(color = "Condition", y = "Percent of non-naive\nCD8- CD3+ cells") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(
    axis.line = element_blank(),
    axis.text.y = element_markdown(size = BASEPTAXISFONTSIZE, color = "#000000"),
    axis.title.y = element_text(size = BASEPTFONTSIZE - 2, color = "#000000"),
    axis.text.x = element_markdown(size = BASEPTAXISFONTSIZE - 2, color = "#000000"),
    legend.position = "bottom",
    legend.key.size = unit(5, "pt"),
    panel.grid.major.y = element_line(color = "#eeeeee", linetype = "dotted"),
    legend.text = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    legend.title = element_text(size = BASEPTFONTSIZE, color = "#000000"),
    axis.title.x = element_blank())

organoid_df %>%
  dplyr::group_by(condition, p24) %>%
  dplyr::reframe(pp = percent / sum(percent) * 100) %>%
  print(n = 200)

################################################################################
# final plots
################################################################################
layout <- c(
  area(1, 1, 8, 9), # prop
  area(9, 1, 10, 2), # prdm1
  area(9, 3, 10, 5), # depletion
  area(9, 6, 10, 7), # prdm1 dea
  area(9, 8, 10, 9), # prdm1 deg
  area(11, 1, 13, 6), # representative plot
  area(11, 6, 13, 9)
)

p <- patchwork::free(propClusterG + theme(plot.margin = margin(l = 5, r = 25)), type = "panel", side = "lr") +
  patchwork::free(prdmByDonorG + theme(plot.margin = margin(t = 15)), type = "panel", side = "lr") +
  patchwork::free(depletionDfG + theme(plot.margin = margin(t = 15)), type = "panel", side = "lr") +
  patchwork::free(dea_blimp1G +
      theme(
        plot.margin = margin(t = 15),
        plot.tag = element_text(margin = margin(l = 5))), type = "panel", side = "tlbr") +
  patchwork::free(deg_blimp1G, type = "panel", side = "tlbr") + theme(plot.margin = margin(t = 15)) +
  wrap_elements(full = rep_flow_gg + theme(plot.margin = margin(l = -100, r = 0))) +
  free(organoid_g + 
      theme(
        plot.margin = margin(l = 0),
        plot.tag = element_text(margin = margin(l = 5))), type = "panel", side = "tlr") +
  plot_annotation(tag_levels = "A") +
  plot_layout(design = layout)

saveFinalFigure(plot = p, fn = "fig5_v2", prefixDir = "figures/outs", devices = c("png"), gwidth = 8.5, gheight = 10)
