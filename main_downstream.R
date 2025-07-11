library(ArchR)
suppressMessages(library(Seurat))
suppressMessages(library(glue))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(ggrastr))
library(BSgenome.Hsapiens.UCSC.hg38)
library(ComplexHeatmap)
# 
# dyn.load("/usr/local/hdf5/lib/libhdf5_hl.so.310")
# library(hdf5r)

set.seed(42)

addArchRThreads(threads = 12)
addArchRGenome("hg38")


source("scripts/helperUtils.R")
source("scripts/addHaystackData.R")
source("scripts/colorPalettes.R")
source("scripts/graphs.R")
source("scripts/clusterAnnotation.R")

################################################################################
# data directories
################################################################################
cellrangerOutDir <- "../hiv_ln_dogma_data/cellranger_out/"
haystackOutDir <- "../hiv_ln_dogma_data/haystack_out"
tsaCatalog <- readRDS("rds/tsa_catalog.rds")

################################################################################
# load project
################################################################################
projUnhashed <- loadArchRProject(
  path = "dogma_lsi_Harmony_v2",
  showLogo = FALSE)

seuMerged <- readRDS("rds/adtHto_matched_final.rds")

################################################################################
# load haystack
################################################################################
samplesForHaystack <- as.character(c(1:105))
names(samplesForHaystack) <- samplesForHaystack

haystackResult <- addHaystackData(
  proj = projUnhashed,
  haystackParentDir = "../hiv_ln_dogma_data/haystack_out",
  haystackSamples = samplesForHaystack,
)

projUnhashed <- haystackResult$newProj

################################################################################
# clean donor ID
################################################################################
tmp <- gsub("(B03-16|G-LFSA|\\D+)", "", projUnhashed$donor)
projUnhashed$donorClean <- tmp
projUnhashed$stage <- str_match(projUnhashed$pool, "(ART|Chronic|HIV-)")[, 1]
projUnhashed$donorClean2 <- case_when(
  projUnhashed$stage == "ART" ~ "A",
  projUnhashed$stage == "Chronic" ~ "C",
  projUnhashed$stage == "HIV-" ~ "U",
)
projUnhashed$donorClean2 <- paste0(projUnhashed$donorClean2, projUnhashed$donorClean)
projUnhashed$tissue <- str_match(projUnhashed$pool, "(CLN|ILN| LN)")[, 1]
projUnhashed$tissue <- ifelse(projUnhashed$tissue == " LN", "ILN", projUnhashed$tissue)


################################################################################
# add gene expression back to seurat
# WAIT those values are already normalized
# need to input from 10x file instead...
################################################################################
samplesMinus98 <- samplesForHaystack[samplesForHaystack != "98"]
names(samplesMinus98) <- NULL

spMatRna <- lapply(samplesMinus98, function(x) {
  message(x)

  mat <- Read10X_h5(paste0(cellrangerOutDir, x, "/outs/filtered_feature_bc_matrix.h5"))
  mat2 <- mat$`Gene Expression`

  colnames(mat2) <- paste0(x, "#", colnames(mat2))
  return(mat2)
})

spMatRnaFinal <- BiocGenerics::Reduce(cbind, spMatRna)
spMatRnaFinal <- spMatRnaFinal[, projUnhashed$cellNames]

rnaAssay <- CreateAssayObject(counts = spMatRnaFinal)

# TODO for some reason this causes issue with R hanging but this may solve it?
detach("package:hdf5r")
seuMerged[["RNA"]] <- rnaAssay

seuMerged <- NormalizeData(seuMerged, assay = "RNA")
seuMerged <- NormalizeData(seuMerged, assay = "ADT", normalization.method = "CLR", margin = 2, scale.factor = 10000)

transferToSeu <- data.frame(
  stage = projUnhashed$stage,
  donorClean2 = projUnhashed$donorClean2,
  tissue = projUnhashed$tissue,
  sample = projUnhashed$Sample,
  gex_UMI = projUnhashed$Gex_nUMI
)
rownames(transferToSeu) <- projUnhashed$cellNames

seuMerged <- AddMetaData(seuMerged, transferToSeu)

VlnPlot(seuMerged, c("nCount_ADT", "nCount_HTO", "gex_UMI"), ncol = 1, group.by = "sample", pt.size = 0) &
#   scale_y_log10() &
#   theme(legend.position = "blank") &
#   geom_hline(yintercept = 100) &
#   annotation_logticks(sides = "l")
# 
#  VlnPlot(seuMerged, c("nFeature_ADT"), ncol = 1, group.by = "sample", pt.size = 0) +
#   theme(legend.position = "blank")
# 
# data.frame(
#   ncount = seuMerged$nCount_ADT,
#   sample = seuMerged$sample) %>%
#   group_by(sample) %>%
#   summarize(avg = mean(ncount),
#     sd = sd(ncount)) %>%
#   print(n = 80)
# 
# 
rm(spMatRna)
rm(spMatRnaFinal)
rm(transferToSeu)

gc()


currentTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
preAnnotFn <- glue("rdata/preLSI_{currentTime}.rdata")
save.image(preAnnotFn)


################################################################################
# START HERE after loading libraries.
################################################################################
load("rdata/preLSI_2024-11-16_02-16-12.rdata")

# NOTE RELOAD ALL SCRIPT FILES


################################################################################
# clean out APC/NK/B
################################################################################
DefaultAssay(seuMerged) <- "RNA"

nonTMarkers <- list(
  "b" = c("MS4A1", "IGKC", "IGLC1", "IGHM", "IGHG1", "IGHG2", "IGHD", "IGHE", "IGHA1", "IGHA2"),
  "apc_nk" = c("FCGR3A", "CD14", "LYZ", "MS4A7", "FCER1A", "CST3", "CR1", "CR2", "MFGE8", "CXCL13")
)

seuMerged <- AddModuleScore(seuMerged, features = nonTMarkers, name = "nonT_")

# will use a cutoff of 0.05 for b
RidgePlot(seuMerged, features = c("nonT_1"), group.by = "sample") +
  scale_x_continuous(limits = c(-0.1, 1)) +
  geom_vline(xintercept = 0.05)

seuMerged$isBCell <- seuMerged$nonT_1 >= 0.05

# will use a cutoff of 0.05 for apc/nk
RidgePlot(seuMerged, features = c("nonT_2"), group.by = "sample") +
  scale_x_continuous(limits = c(-0.1, 1)) +
  geom_vline(xintercept = 0.05)

seuMerged$isApcNk <- seuMerged$nonT_2 >= 0.05

# filter atac project
cellsToKeep <- names(seuMerged$orig.ident)[!seuMerged$isBCell & !seuMerged$isApcNk]
projUnhashedClean <- subsetArchRProject(
  cells = cellsToKeep,
  ArchRProj = projUnhashed,
  outputDirectory = "dogma_postFilter_v2",
  force = TRUE
)

rm(cellsToKeep)
rm(projUnhashed)
rm(rnaAssay)
gc()

currentTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
postFiltCheckptFn <- glue("rdata/preLSI_B-APC-NK_filt_{currentTime}.rdata")
save.image(postFiltCheckptFn)

# load("rdata/preLSI_B-APC-NK_filt_2024-11-16_21-30-11.rdata")

# UMAP
projUnhashedClean <- addIterativeLSI(
  ArchRProj = projUnhashedClean,
  useMatrix = "TileMatrix",
  depthCol = "nFrags",
  saveIterations = FALSE,
  name = "LSI_ATAC_clean",
  iterations = 3,
  firstSelection = "Top",
  varFeatures = 20000,
  dimsToUse = 1:30,
  sampleCellsPre = 75000,
  threads = 16,
  sampleCellsFinal = 300000, # this arg is req'd because of Matrix size limits.
  force = TRUE, # comment out when done with this
  clusterParams = list(
    resolution = 2, 
    sampleCells = 25000,
    maxClusters = 25
  )
)

projUnhashedClean <- addIterativeLSI(
  ArchRProj = projUnhashedClean, 
  clusterParams = list(
    resolution = 2, 
    sampleCells = 30000,
    maxClusters = 25
  ),
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix", 
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  sampleCellsPre = 75000,
  sampleCellsFinal = 300000, # this arg is req'd because of Matrix size limits.
  force = TRUE, # comment out when done with this
  firstSelection = "Var",
  binarize = FALSE,
  iterations = 3,
  name = "LSI_RNA_clean"
)

projUnhashedClean <- addHarmony(
  ArchRProj = projUnhashedClean,
  reducedDims = "LSI_ATAC_clean",
  name = "Harmony_LSI_ATAC_clean",
  groupBy = c("donor", "chip"),
  verbose = TRUE,
  force = TRUE,
  ncores = 10
)

projUnhashedClean <- addHarmony(
  ArchRProj = projUnhashedClean,
  reducedDims = "LSI_RNA_clean",
  name = "Harmony_LSI_RNA_clean",
  groupBy = c("donor", "chip"),
  force = TRUE,
  ncores = 10
)

# Combined Dims
projUnhashedClean <- addCombinedDims(
  projUnhashedClean,
  reducedDims = c("Harmony_LSI_ATAC_clean", "Harmony_LSI_RNA_clean"),
  name = "Harmony_LSI_Combined_clean")


# order the combined dims by variance.
projUnhashedClean@reducedDims$HarmonyLsiCleanForSeurat <- projUnhashedClean@reducedDims$Harmony_LSI_Combined_clean

# taken from ArchR MultiModal.R to calculate variance by column
# will sort these by variance explained in case FindNeighbors or FindClusters 
# are influenced by order (since FindClusters will throw an error if 
# the colnames are not labeled in sequential order)
cV <- apply(projUnhashedClean@reducedDims$HarmonyLsiCleanForSeurat$matRD, 2, var) 

projUnhashedClean@reducedDims$HarmonyLsiCleanForSeurat$matRD <- 
  projUnhashedClean@reducedDims$HarmonyLsiCleanForSeurat$matRD[, order(cV, decreasing=TRUE)]

colnames(projUnhashedClean@reducedDims$HarmonyLsiCleanForSeurat$matRD) <-
  paste0("tmp", c(1: ncol(projUnhashedClean@reducedDims$HarmonyLsiCleanForSeurat$matRD)))

projUnhashedClean <- addUMAP(projUnhashedClean,
  reducedDims = "HarmonyLsiCleanForSeurat",
  name = "Harmony_UMAP_Combined_v2",
  minDist = 0.1,
  nNeighbors = 20,
  force = TRUE)

projUnhashedClean <- addUMAP(projUnhashedClean,
  reducedDims = "Harmony_LSI_RNA_clean",
  name = "Harmony_UMAP_RNA_v2",
  minDist = 0.2,
  nNeighbors = 20,
  force = TRUE)

projUnhashedClean <- addUMAP(projUnhashedClean,
  reducedDims = "Harmony_LSI_ATAC_clean",
  name = "Harmony_UMAP_ATAC_v2",
  minDist = 0.2,
  nNeighbors = 20,
  force = TRUE)


# Initial clustering
projUnhashedClean <- addClusters(
  input = projUnhashedClean,
  reducedDims = "HarmonyLsiCleanForSeurat",
  name = "Clusters_v2",
  method = "Seurat",
  resolution = 0.6,
  force = TRUE
)

projUnhashedClean <- addClusters(
  input = projUnhashedClean,
  reducedDims = "Harmony_LSI_RNA_clean",
  name = "Clusters_RNA_v2",
  method = "Seurat",
  resolution = 0.6,
  force = TRUE
)

projUnhashedClean <- addClusters(
  input = projUnhashedClean,
  reducedDims = "Harmony_LSI_ATAC_clean",
  name = "Clusters_ATAC_v2",
  method = "Seurat",
  resolution = 0.6,
  force = TRUE
)

projUnhashedClean <- addImputeWeights(projUnhashedClean, reducedDims = "HarmonyLsiCleanForSeurat")

currentTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
postFiltCheckpt2Fn <- glue("rdata/datasetV2_preAnnot_B-APC-NK_filt2_{currentTime}.rdata")
save.image(postFiltCheckpt2Fn)


################################################################################
# begin initial graphs
################################################################################
load("rdata/datasetV2_preAnnot_B-APC-NK_filt2_2024-11-18_21-29-43.rdata")

################################################################################
# get hiv rna counts
################################################################################
hxb2Bam <- list.files("../hiv_ln_dogma_data/gex_hxb2_bam/", pattern = "*.bam", full.names = TRUE)

hivBamList <- lapply(hxb2Bam, function(x) {
  if (file.size(x) == 0L) {
    return(NULL)
  }
  
  return(Rsamtools::scanBam(x,
    param = Rsamtools::ScanBamParam(what = Rsamtools::scanBamWhat(), tag = c("NH", "HI", "AS", "nM", "CB"))))
})

hivBamDf <- lapply(hivBamList, function(x) {
  y <- x[[1]]
  
  df <- data.frame(
    qname = y$qname,
    flag = y$flag,
    rname = y$rname,
    pos = y$pos,
    qwidth = y$qwidth,
    mapq = y$mapq,
    cigar = y$cigar,
    seq = as.character(y$seq),
    qual = y$qual,
    nh = y$tag$NH,
    cb = y$tag$CB
  )
  
  return(df)
})

names(hivBamDf) <- as.character(str_match(hxb2Bam, "(\\d+).hxb2.bam")[, 2])
hivBamDf <- bind_rows(hivBamDf, .id = "Sample")

hivBamDfClean <- hivBamDf %>%
  filter(!is.na(cb)) %>%
  group_by(Sample, cb, flag, rname, pos, qwidth, cigar, nh) %>%
  summarize(
    total = n(),
    seq = ifelse(length(unique(seq)) > 1, paste0(seq, collapse = ","), c(first(seq)))) %>%
  group_by(Sample, cb, rname, pos, qwidth, cigar, nh) %>%
  group_split()


hivBamDfClean <- lapply(hivBamDfClean, function(x) {
  nFlag1024 <- 0
  if (nrow(x[x$flag == 1024, ]) > 1) {
    message("Multiple 1024s detected.")
  } else if (nrow(x[x$flag == 1024, ]) == 1) {
    nFlag1024 <- as.numeric(x[x$flag == 1024, "total"])
  }
  
  nFlag16 <- 0
  
  if (nrow(x[x$flag == 16, ]) > 1) {
    message("Multiple 16s detected.")
  } else if (nrow(x[x$flag == 16, ]) == 1) {
    nFlag16 <- as.numeric(x[x$flag == 16, "total"])
  }
  
  uniqSeq <- unique(x$seq)
  if (length(uniqSeq) > 1) {
    seqToAdd <- paste0(x$seq, collapse = ",")
  } else {
    seqToAdd <- uniqSeq
  }
  
  y <- x[1, ]
  df <- data.frame(
    sample = y$Sample,
    cb = y$cb,
    rname = y$rname,
    pos = y$pos,
    qwidth = y$qwidth,
    cigar = y$cigar,
    nh = y$nh,
    seq = seqToAdd
  )
  
  df$nFlag0 <- ifelse(nrow(x[x$flag == 0, ] != 0), 1, 0)
  df$nFlag1024 <- nFlag1024
  df$nFlag16 <- nFlag16
  
  return(df)
})

# no multiple 1024 or 16 detected per unique cbc/pos/cigar set.
hivBamDfClean <- bind_rows(hivBamDfClean)

hivBamDfClean <- hivBamDfClean %>%
  mutate(cellName = paste0(sample, "#", cb))

# ya-chi metric is "cells with a minimum of 2 HIV-1-related UMI or
#  4 reads of a single UMI were considered positive.”
hivRnaFinal <- hivBamDfClean %>%
  filter(cellName %in% projUnhashedClean$cellNames) %>%
  group_by(cellName) %>%
  mutate(
    isGenuineBy2SeparateUMI = n() >= 2,
    isGenuineBy4Reads = nFlag0 + nFlag1024 + nFlag16 >= 4)

hivRnaUmiCount <- hivRnaFinal %>%
  mutate(pseudoFlag0 = ifelse(nFlag0 == 0 & nFlag1024 > 0 & nFlag16 == 0, 1, 0)) %>%
  group_by(cellName) %>%
  summarize(totalUmi = sum(pseudoFlag0, nFlag0, nFlag16))

matchedHivRnaPos <- hivRnaFinal %>%
  filter(isGenuineBy2SeparateUMI | isGenuineBy4Reads) %>%
  distinct(cellName)

hivRnaUmiCount <- as.data.frame(hivRnaUmiCount)
rownames(hivRnaUmiCount) <- hivRnaUmiCount$cellName

projUnhashedClean$hivRNA <- FALSE
projUnhashedClean$hivRNA[projUnhashedClean$cellNames %in% matchedHivRnaPos$cellName] <- TRUE

hivStats <- data.frame(
  cellName = projUnhashedClean$cellNames,
  donor = projUnhashedClean$donorClean2,
  stage = projUnhashedClean$stage,
  hivDNA = projUnhashedClean$haystackOut,
  hivRNA = projUnhashedClean$hivRNA,
  hivRNAUmi = NA
)

hivStats[hivStats$hivRNA, "hivRNAUmi"] <- hivRnaUmiCount[hivStats[hivStats$hivRNA, "cellName"], "totalUmi"]
projUnhashedClean$hivRNAUmi <- hivStats$hivRNAUmi

# some new upset plots..
statGroupingDict <- hivStats$donor
names(statGroupingDict) <- hivStats$cellName

plotUpsetHivStats(
  rnaCbc = hivStats[hivStats$hivRNA & hivStats$stage == "ART", "cellName"], 
  dnaCbc = hivStats[hivStats$hivDNA & hivStats$stage == "ART", "cellName"],
  groupingVector = statGroupingDict,
  fn = "upsetHivStats_art",
  devices = "png",
  gheight = 4)

plotUpsetHivStats(
  rnaCbc = hivStats[hivStats$hivRNA & hivStats$stage == "Chronic", "cellName"], 
  dnaCbc = hivStats[hivStats$hivDNA & hivStats$stage == "Chronic", "cellName"],
  groupingVector = statGroupingDict,
  fn = "upsetHivStats_chronic",
  devices = "png",
  gheight = 4)

# aggregate
statGroupingDict <- hivStats$stage
names(statGroupingDict) <- hivStats$cellName

plotUpsetHivStats(
  rnaCbc = hivStats[hivStats$hivRNA, "cellName"], 
  dnaCbc = hivStats[hivStats$hivDNA, "cellName"],
  groupingVector = statGroupingDict,
  fn = "upsetHivStats_aggByStage",
  devices = "png",
  gheight = 4)


################################################################################
# annotation time!
################################################################################
DONORCOLORS <- multiHueColorPalette[c(1:length(unique(projUnhashedClean$donorClean2)))]
names(DONORCOLORS) <- unique(projUnhashedClean$donorClean2)

CLUSTERCOLORS <- multiHueColorPalette[c(1:length(unique(projUnhashedClean$Clusters_v2)))]
names(CLUSTERCOLORS) <- unique(projUnhashedClean$Clusters_v2)

# rebuild umap
projUnhashedClean <- addUMAP(
  projUnhashedClean,
  reducedDims = "HarmonyLsiCleanForSeurat",
  name = "Harmony_UMAP_Combined_v2_1Dist_20Neigh",
  minDist = 0.1,
  nNeighbors = 20,
  force = TRUE)

projUnhashedClean <- addUMAP(
  projUnhashedClean,
  reducedDims = "HarmonyLsiCleanForSeurat",
  name = "Harmony_UMAP_Combined_v2_05Dist_50Neigh",
  minDist = 0.05,
  nNeighbors = 50,
  force = TRUE)

# plot where each cluster is for help with annotation...
clusterNumericOrder <- str_sort(unique(projUnhashedClean$Clusters_v2), numeric = TRUE)

seuMerged <- subset(seuMerged, subset = cellName %in% projUnhashedClean$cellNames)

transferToSeu <- data.frame(
  Clusters_v2 = projUnhashedClean$Clusters_v2
)

rownames(transferToSeu) <- projUnhashedClean$cellNames
seuMerged <- AddMetaData(seuMerged, transferToSeu)

cellDeathOfInterest <- c(
  "NR4A1",
  "NR4A2",
  "BCL2", # apoptosis regulator
  "BCL2L1",
  "MCL1",
  "BAK1", # BCL2 antagonist killer 1
  "ATF3",
  "DUSP5",
  "DIABLO",
  "DDIT3",
  "EIF2AK3",
  "OSGIN1",
  "TNFRSF25", # DR3
  "TNFRSF10B", # DR5
  "CASP8",
  "TNFSF10", # trail
  "TRADD",
  "TP53", # p53
  "BAX",
  "BBC3",
  "PMAIP1",
  "SERPINB9"
)

baseADTMarkersForAnotation <- c(
  "A0034", #CD3
  "A0072", #CD4
  "A0046", #CD8
  "A0081", #CD14
  "A0083", #CD16
  "A0047", #CD56
  "A0050", #CD19
  "A0087", #CD45RO
  "A0063", #CD45RA
  "A0154", #CD27
  "A0386", #CD28
  "A0156", #CD95
  "A0390", #CD127
  "A0146", #CD69
  "A0085", #CD25
  "A0089", #TIGIT
  "A0088", #PD1
  "A0141", #CCR5
  "A0144", #CXCR5
  "A0149" #CD161
)

# add hiv status based on hiv dna and/or rna
hivInfo <- data.frame(
  cbc = projUnhashedClean$cellNames,
  dna = projUnhashedClean$haystackOut,
  rna = projUnhashedClean$hivRNA
) %>%
  mutate(status = case_when(
    dna & !rna ~ "HIV DNA+ RNA-",
    rna ~ "HIV DNA+/- RNA+",
    TRUE ~ "HIV-"
  )) %>%
  select(cbc, status)

projUnhashedClean <- addCellColData(
  ArchRProj = projUnhashedClean,
  data = hivInfo$status,
  cells = hivInfo$cbc,
  name = "hivStatus",
  force = TRUE)

seuMerged$hivStatus <- hivInfo$status


# discrete counts
tmp <- plotDiscreteLollipop2(
  proj = projUnhashedClean,
  colorVar = "hivStatus",
  colorScheme = COLORS[["hivStatus"]],
  groupingVar = "Clusters_v2",
  splittingVar = "stage",
  splittingVarOrder = c("HIV-", "Chronic", "ART")
)

(tmp[[1]] + theme(legend.position = "blank")) +
  (tmp[[2]] + theme(legend.position = "bottom")) +
  (tmp[[3]] + theme(legend.position = "bottom")) +
  plot_layout(guides = "collect", nrow = 1) +
  plot_annotation(theme = theme(legend.position = "bottom")) &
  labs(color = "Cell infection status") &
  theme(
    axis.text = element_text(color = "#000000", size = 10),
    plot.subtitle = element_text(color = "#000000", size = 14))

# NOTE TISSUE WAS ASSIGNED INCORRECTLY IN EARLIER CODE. FIXING IT HERE.
tissueAssigner <- data.frame(
  cbc = projUnhashedClean$cellNames,
  pool = projUnhashedClean$pool,
  donorClean2 = projUnhashedClean$donorClean2
)

tissueAssigner <- tissueAssigner %>%
  mutate(
    tissue = case_when(
      donorClean2 == "A024" ~ "ILN",
      donorClean2 == "A025" ~ "ILN",
      grepl("CLN", pool) ~ "CLN",
      grepl("ILN", pool) ~ "ILN",
      pool == "Chronic Mix 6 LN" ~ "CLN",
      pool == "Chronic Mix 3 LN" ~ "CLN",
      pool == "ART Mix 5 LN" ~ "CLN",
      pool == "ART Mix 3 LN" ~ "CLN",
      pool == "ART Mix 2 LN + Dump" ~ "CLN",
      TRUE ~ "Unknown"
    )
  )

projUnhashedClean <- addCellColData(
  ArchRProj = projUnhashedClean,
  data = tissueAssigner$tissue,
  cells = tissueAssigner$cbc,
  name = "tissue",
  force = TRUE)

seuMerged$tissue <- tissueAssigner$tissue

currentTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
preAnnotCheckptFn <- glue("rdata/preManualAnnot_finalCheckpoint_{currentTime}.rdata")
save.image(preAnnotCheckptFn)

################################################################################
# load checkpoint and reload script files
################################################################################
load("rdata/preManualAnnot_finalCheckpoint_2024-11-21_20-45-01.rdata")

source("scripts/helperUtils.R")
source("scripts/addHaystackData.R")
source("scripts/colorPalettes.R")
source("scripts/graphs.R")
source("scripts/clusterAnnotation.R")

################################################################################
# check clusterings
################################################################################
projUnhashedClean <- addClusters(
  input = projUnhashedClean,
  reducedDims = "HarmonyLsiCleanForSeurat",
  name = "Clusters_v3",
  method = "Seurat",
  resolution = 0.8,
  maxClusters = NULL,
  force = TRUE
)

# adding in mono cassette just for double checking with rna plot
baseGenesForAnnotation_v2 <- read.table("manualAnnot/baseGenesForAnnotation_v2.txt")$V1

transferToSeu <- data.frame(
  Clusters_v3 = projUnhashedClean$Clusters_v3
)

rownames(transferToSeu) <- projUnhashedClean$cellNames
seuMerged <- AddMetaData(seuMerged, transferToSeu)

clusterNumericOrder_v3 <- str_sort(unique(projUnhashedClean$Clusters_v3), numeric = TRUE)

# base
DotPlot(seuMerged,
  dot.scale = 10,
  col.max = 3,
  col.min = -3,
  features = paste0("rna_", baseGenesForAnnotation_v2),
  group.by = "Clusters_v3") +
  scale_x_discrete(labels = ~ gsub("rna_", "", .x)) +
  scale_y_discrete(limits = clusterNumericOrder_v3) +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    panel.grid.major.y = element_line(color = "#efefef"))


atacGeneScoreForAnnot_v3 <- tryCatch({
  tmp <- ArchR:::.getMatrixValues(
    ArchRProj = projUnhashedClean,
    name = baseGenesForAnnotation_v2,
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

atacGeneScoreForAnnotImputed_v3 <- imputeMatrix(mat = as.matrix(atacGeneScoreForAnnot_v3),
  imputeWeights = getImputeWeights(projUnhashedClean))

atacGeneScoreForAnnotImputed_v3 <- log2(atacGeneScoreForAnnotImputed_v3 + 1)

atacImputedDf_v3 <- data.frame(
  cbc = projUnhashedClean$cellNames,
  unannot_v3 = projUnhashedClean$Clusters_v3
)
rownames(atacImputedDf_v3) <- atacImputedDf_v3$cbc
atacImputedDf_v3[colnames(atacGeneScoreForAnnotImputed_v3), baseGenesForAnnotation_v2] <- t(atacGeneScoreForAnnotImputed_v3)

atacImputedDf_v3 <- atacImputedDf_v3 %>%
  pivot_longer(cols = -any_of(c("cbc", "unannot_v3")), names_to = "feature", values_to = "exp")

atacImputePlot_v3 <- atacImputedDf_v3 %>%
  filter(feature %in% baseGenesForAnnotation_v2[c(41:58)]) %>%
  ggplot(aes(x = unannot_v3, y = exp)) +
  geom_violin(aes(color = unannot_v3), width = 1, fill = "#ffffff00") +
  facet_wrap(~ feature, scales = "free_y", nrow = 10) +
  scale_x_discrete(limits = clusterNumericOrder_v3) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 8),
    legend.position = "blank",
    panel.grid = element_blank(),
    axis.title.x = element_blank()) +
  labs(
    y = "imputed log2(counts + 1)")

# adt plot
baseRidgePlot <- RidgePlot(
  seuMerged,
  features = baseADTMarkersForAnotation,
  group.by = "Clusters_v3",
  combine = FALSE)

baseRidgePlotGrid <- lapply(baseRidgePlot, function(x) {
  tsaID <- x$labels$title
  tsaID <- gsub("adt_", "", tsaID)
  
  tsaCDName <- tsaCatalog[tsaCatalog$DNA_ID == tsaID, "cleanName"]
  
  # plotting up to the 98th percentile for ease of visualization only because
  # of outliers.
  newXMax <- quantile(x$data[, 1], probs = 0.98)
  newXMax <- ceiling(newXMax)
  
  tmp <- x +
    guides(color = "none", fill = "none") +
    labs(title = glue("{tsaCDName} ({tsaID})")) +
    scale_y_discrete(limits = clusterNumericOrder_v3) +
    scale_x_continuous(limits = c(-0.25, as.numeric(newXMax)), expand = expansion(mult = c(0.05, 0.3))) +
    theme(axis.title = element_blank(),
      axis.text = element_text(size = 6),
      panel.grid.major.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title.position = "plot",
      plot.title = element_text(size = 8))
  
  return(tmp)
})

baseRidgePlot <- do.call(cowplot::plot_grid, c(list(ncol = 5), baseRidgePlotGrid))
savePlot(plot = baseRidgePlot, fn = "base_adt_v2", devices = c("png"), gheight = 10, gwidth = 13)


################################################################################
# add in manual annotation again
################################################################################
annotDf <- read.csv("manualAnnot/dataset2_clusters_v3_manualAnnot.csv", header = FALSE)

manualAnnots <- annotDf[, 2]
names(manualAnnots) <- annotDf[, 1]

projUnhashedClean <- addCellColData(
  ArchRProj = projUnhashedClean,
  data = manualAnnots[getCellColData(projUnhashedClean, select = "Clusters_v3")[, 1]],
  cells = projUnhashedClean$cellNames,
  name = "manualAnnot",
  force = TRUE)

seuMerged$manualAnnot <- unname(manualAnnots[seuMerged$Clusters_v3])

ANNOTCLUSTERCOLORS <- CLUSTERCOLORS
names(ANNOTCLUSTERCOLORS) <- manualAnnots[names(CLUSTERCOLORS)]

currentTime <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
postAnnotCheckptFn <- glue("rdata/postManualAnnot_{currentTime}.rdata")
save.image(postAnnotCheckptFn)

################################################################################
# downstream analysis
################################################################################
seuMerged$hivPos <- ifelse(seuMerged$hivStatus != "HIV-", "HIV+", "HIV-")
projUnhashedClean$hivPos <- ifelse(seuMerged$hivStatus != "HIV-", "HIV+", "HIV-")


CLUSTERCOLORS <- multiHueColorPalette[c(1:length(unique(projUnhashedClean$Clusters_v3)))]
names(CLUSTERCOLORS) <- unique(projUnhashedClean$Clusters_v3)

ANNOTCLUSTERCOLORS <- CLUSTERCOLORS
names(ANNOTCLUSTERCOLORS) <- manualAnnots[names(CLUSTERCOLORS)]

projUnhashedClean$isCD4 <- grepl("^CD4", projUnhashedClean$manualAnnot)
seuMerged$isCD4 <- grepl("^CD4", seuMerged$manualAnnot)
seuMerged$haystackOut <- projUnhashedClean$haystackOut

manualAnnotBins <- data.frame(manualAnnot = seuMerged$manualAnnot) %>%
  mutate(manualAnnotBin = case_when(
    grepl("(Tcm|Th2)", manualAnnot) ~ "Tcm/Th2",
    grepl("Tem", manualAnnot) ~ "Tem/Th1",
    grepl("Treg", manualAnnot) ~ "Treg",
    grepl("Trm", manualAnnot) ~ "Trm",
    grepl("Tfh", manualAnnot) ~ "Tfh",
    grepl("Th17", manualAnnot) ~ "Th17",
    grepl("(early memory|Tn|Tscm)", manualAnnot) ~ "Naive/early memory",
    TRUE ~ "Other"
  ))

seuMerged <- AddMetaData(seuMerged, metadata = manualAnnotBins)
projUnhashedClean$manualAnnotBin <- manualAnnotBins$manualAnnotBin

################################################################################
# art epigenetic analysis
################################################################################
projARTsub <- projUnhashedClean[!grepl("Other", projUnhashedClean$manualAnnotBin), ]
projARTsub <- projARTsub[projARTsub$stage == "ART", ]
projARTsub$manualAnnotBinWithHiv <- paste0(projARTsub$manualAnnotBin, "_", projARTsub$hivPos)

projARTsub <- addGroupCoverages(
  ArchRProj = projARTsub,
  groupBy = "manualAnnotBinWithHiv",
  minReplicates = 2,
  minCells = 20,
  force = TRUE,
  sampleRatio = 0.9)

projARTsub <- addReproduciblePeakSet(
  ArchRProj = projARTsub,
  groupBy = "manualAnnotBinWithHiv", 
  pathToMacs2 = "/home/wuv/anaconda3/envs/asapseq/bin/macs2")


projARTsub <- addPeakMatrix(projARTsub, force = TRUE)
projARTsub <- addMotifAnnotations(ArchRProj = projARTsub, motifSet = "cisbp", name = "Motif", force = TRUE)
projARTsub <- addBgdPeaks(projARTsub, force = TRUE)


gc()
addArchRThreads(threads = 4)
projARTsub <- addDeviationsMatrix(
  ArchRProj = projARTsub, 
  peakAnnotation = "Motif",
  force = TRUE
)

gc()
addArchRThreads(threads = 16)
projARTsub <- saveArchRProject(
  projARTsub,
  outputDirectory = "dogma_filt_art_gc",
  load = TRUE,
  dropCells = TRUE)


################################################################################
# chronic epigenetic analysis
################################################################################
projChronicSub <- projUnhashedClean[!grepl("Other", projUnhashedClean$manualAnnotBin), ]
projChronicSub <- projChronicSub[projChronicSub$stage == "Chronic", ]
projChronicSub$manualAnnotBinWithHiv <- paste0(projChronicSub$manualAnnotBin, "_", projChronicSub$hivPos)

addArchRThreads(threads = 12)
projChronicSub <- addGroupCoverages(
  ArchRProj = projChronicSub,
  groupBy = "manualAnnotBinWithHiv",
  minReplicates = 2,
  minCells = 20,
  force = TRUE,
  sampleRatio = 0.9)

projChronicSub <- addReproduciblePeakSet(
  ArchRProj = projChronicSub,
  groupBy = "manualAnnotBinWithHiv", 
  pathToMacs2 = "/home/wuv/anaconda3/envs/asapseq/bin/macs2")

projChronicSub <- addPeakMatrix(projChronicSub, force = TRUE)
projChronicSub <- addMotifAnnotations(ArchRProj = projChronicSub, motifSet = "cisbp", name = "Motif", force = TRUE)
projChronicSub <- addBgdPeaks(projChronicSub, force = TRUE)


gc()
addArchRThreads(threads = 10)
projChronicSub <- addDeviationsMatrix(
  ArchRProj = projChronicSub, 
  peakAnnotation = "Motif",
  force = TRUE
)

gc()
addArchRThreads(threads = 24)
projChronicSub <- saveArchRProject(
  projChronicSub,
  outputDirectory = "dogma_filt_chronic_gc",
  load = TRUE,
  dropCells = TRUE)

addArchRThreads(threads = 8)
chronic_markerTest_deviations <- getMarkerFeatures(
  ArchRProj = projChronicSub, 
  useMatrix = "MotifMatrix",
  groupBy = "manualAnnotBinWithHiv",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  verbose = TRUE,
  useSeqnames = "z",
)

save.image(format(Sys.time(), "rdata/chronicSub_deviations_%Y-%m-%d-%H-%M.Rdata"))

# chronic but separate from tx active and inactive
projUnhashedClean$manualAnnotBinWithHivStatus <- paste0(projUnhashedClean$manualAnnotBin, "_", projUnhashedClean$hivStatus)
projChronicSub2 <- projUnhashedClean[projUnhashedClean$stage == "Chronic" & projUnhashedClean$manualAnnotBin %in% c("Tcm/Th2", "Tem/Th1", "Tfh", "Treg"), ]


addArchRThreads(threads = 12)
projChronicSub2 <- addGroupCoverages(
  ArchRProj = projChronicSub2,
  groupBy = "manualAnnotBinWithHivStatus",
  minReplicates = 2,
  minCells = 30,
  force = TRUE,
  sampleRatio = 0.95)

projChronicSub2 <- addReproduciblePeakSet(
  ArchRProj = projChronicSub2,
  groupBy = "manualAnnotBinWithHivStatus", 
  pathToMacs2 = "/home/wuv/anaconda3/envs/asapseq/bin/macs2")

projChronicSub2 <- addPeakMatrix(projChronicSub2, force = TRUE)
projChronicSub2 <- addMotifAnnotations(ArchRProj = projChronicSub2, motifSet = "cisbp", name = "Motif", force = TRUE)
projChronicSub2 <- addBgdPeaks(projChronicSub2, force = TRUE)


gc()
addArchRThreads(threads = 8)
projChronicSub2 <- addDeviationsMatrix(
  ArchRProj = projChronicSub2, 
  peakAnnotation = "Motif",
  force = TRUE
)

gc()
addArchRThreads(threads = 8)
projChronicSub2 <- saveArchRProject(
  projChronicSub2,
  outputDirectory = "dogma_filt_chronic_gc_2",
  load = TRUE,
  dropCells = TRUE)


save.image(format(Sys.time(), "rdata/chronicSub2_deviations_%Y-%m-%d-%H-%M.Rdata"))
