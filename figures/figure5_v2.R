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
  geom_text(aes(group = infection, label = percentLabel), vjust = 0.5, hjust = -0.15, position = position_dodge(0.7), size = 1.6, show.legend = FALSE) +
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

prdmByDonorG <- prdmByDonor %>% 
  dplyr::filter(manualAnnot == "CD4 PRDM1+ GC-Tfh") %>%
  ggplot(aes(y = percent, x = infection, fill = infection)) +
  stat_summary(geom = "bar", width = 0.6, fun = mean, fill = "#ffffff00", aes(color = infection)) +
  geom_vline(xintercept = -Inf, linewidth = BASELINEWIDTHSIZE, color = "#000000") +
  geom_hline(yintercept = 0, linewidth = BASELINEWIDTHSIZE, color = "#000000") +
  geom_jitter(color = "#000000", pch = 21, height = 0, alpha = 0.75, size = 0.8, show.legend = FALSE, width = 0.2) +
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


depletionDfG <- data.frame(
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
  dplyr::mutate(
    totalCells = sum(nCells),
    percent = nCells / totalCells * 100,
    manualAnnot = factor(manualAnnot, levels = c("CD4 GC-Tfh #1", "CD4 GC-Tfh #2", "CD4 PRDM1+ GC-Tfh"))) %>%
  ggplot(aes(y = percent, x = manualAnnot, fill = manualAnnot)) +
  geom_line(aes(group = donor), color = "#bbbbbb75", linewidth = 0.3) +
  stat_summary(geom = "bar", width = 0.6, fun = mean, fill = "#ffffff00", aes(color = manualAnnot)) +
  geom_vline(xintercept = -Inf, linewidth = BASELINEWIDTHSIZE, color = "#000000") +
  geom_hline(yintercept = 0, linewidth = BASELINEWIDTHSIZE, color = "#000000") +
  geom_point(color = "#000000", pch = 21, alpha = 0.75, size = 0.8, show.legend = FALSE) +
  stat_summary(fun = mean_se_onlyPos, geom = "crossbar", width = 0.3, linewidth = 0.3) +
  stat_summary(fun.data = mean_se_onlyPosVertical, geom = "linerange", linewidth = 0.3) +
  facet_wrap(~ stage, nrow = 1) +
  coord_cartesian(clip = "off") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.08))) +
  theme_classic() +
  theme(
    plot.margin = margin(r = -1, l = -5),
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
  labs(y = "% of donor's Tfh cells")

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
# DEG
################################################################################
seuMergedTfh <- subset(seuMerged, subset = manualAnnotBin == "Tfh")
seuMergedTfh
deg_blimp1 <- FindMarkers(
  object = subset(seuMerged, subset = stage == "Viremic"),
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
    plot.margin = margin(r = -80, l = -15, t = 0, b = 25)
  )


################################################################################
# testing out other graphs
################################################################################
bab_df <- read.csv("external_data/metadata_with_bab.csv")

propClusterByDonorDf <- data.frame(
  donor = projUnhashedClean$donorClean2,
  stage = projUnhashedClean$stage,
  manualAnnot = projUnhashedClean$manualAnnot
) %>%
  dplyr::filter(stage != "HIV-") %>%
  dplyr::mutate(
    manualAnnot = factor(manualAnnot, levels = uniqClusters),
    stage = factor(stage, levels = uniqStages)
  ) %>%
  dplyr::group_by(donor, manualAnnot, .drop = FALSE) %>%
  dplyr::tally(name = "nCells") %>%
  dplyr::group_by(donor) %>%
  dplyr::mutate(
    totalCells = sum(nCells),
    percent = nCells / totalCells * 100,
    percent = ifelse(is.nan(percent), 0, percent)) %>%
  tidyr::pivot_wider(id_cols = "donor", names_from = "manualAnnot", values_from = "percent")


bab_df_clean <- bab_df %>%
  dplyr::filter(!is.na(yu_gp120_bAUC)) %>%
  dplyr::mutate(sf162_ic50 = ifelse(sf162_ic50 == ">1000", NA, as.numeric(sf162_ic50))) %>%
  dplyr::mutate(mn3_ic50 = ifelse(mn3_ic50 == ">1000", NA, as.numeric(mn3_ic50))) %>%
  dplyr::mutate(mlv_ic50 = ifelse(mlv_ic50 == ">1000", NA, as.numeric(mlv_ic50))) %>%
  dplyr::filter(stage != "HIV neg") %>%
  dplyr::mutate(stage = ifelse(stage == "Chronic", "Viremic", stage)) %>%
  dplyr::mutate(viral_load = ifelse(viral_load == "UNDETECTABLE", NA, as.numeric(gsub(",","",viral_load)))) %>%
  left_join(propClusterByDonorDf, by = c("short_id" = "donor")) %>%
  dplyr::filter(!is.na(`CD4 Tcm/non-resident #1`))


bab_df_clean %>% 
  dplyr::group_by(stage) %>%
  dplyr::summarize(across(`CD4 Tcm/non-resident #1`:`CD4 Trm NR4A1/NR4A2 hi`, ~ cor(.x, yu_gp120_bAUC))) %>%
  tidyr::pivot_longer(cols = c(-1), names_to = "phenotype", values_to = "r_val") %>%
  ggplot(aes(y = phenotype, x = r_val, fill = stage)) +
  facet_wrap(~ stage, nrow = 1) +
  geom_col() +
  geom_vline(xintercept = 0) +
  theme_classic() +
  theme(
    panel.grid.major.y = element_line(color = "#efefef90"),
    legend.position = "blank") +
  labs(title = "Linear correlation of cell proportions to binding AUC")


bab_df_clean %>%
  ggplot(aes(x = yu_gp120_bAUC, y = age, shape = stage, color = stage, fill = stage)) +
  geom_point(size = 0.8, show.legend = FALSE) +
  geom_smooth(method = lm, alpha = 0.3, se = TRUE, linewidth = 0.5) +
  stat_cor(
    aes(label = paste(after_stat(rr.label), after_stat(p.label), sep = "~`,`~")), p.digits = 3,
    label.x = 0.01,
    size = BASEFONTSIZE,
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
    # plot.margin = margin(t = -20),
    plot.title = element_text(color = "#ffffff00", size = BASEPTFONTSIZE + 2),
    legend.position = "bottom",
    # legend.box.margin = margin(t = -5, b = 0),
    legend.direction = "horizontal") +
  labs()


layout <- c(
  area(1, 1, 8, 9), # prop
  area(9, 1, 10, 2), # prdm1
  area(9, 3, 10, 5), # depletion
  area(9, 6, 10, 7), # prdm1 dea
  area(9, 8, 10, 9), # prdm1 deg
  area(11, 1, 13, 9) # spacer
)

p <- patchwork::free(propClusterG, type = "panel", side = "tblr") +
  prdmByDonorG +
  depletionDfG +
  patchwork::free(dea_blimp1G, type = "panel", side = "tbr") +
  patchwork::free(deg_blimp1G, type = "panel", side = "tbr") +
  wrap_elements(full = jpeg::readJPEG("figures/flow.jpg", native = TRUE)) + labs(title = "JaYMe iS InFecTinG CelLS") +
  plot_annotation(tag_levels = "A") +
  plot_layout(design = layout)

saveFinalFigure(plot = p, fn = "fig3_v1", prefixDir = "figures/outs", devices = c("png"), gwidth = 8.5, gheight = 10)
