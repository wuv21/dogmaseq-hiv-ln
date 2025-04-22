# renv::load("/data/hiv_ln_dogma/hiv_ln_dogma_archr")

suppressMessages(library(ArchR))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(rtracklayer))
suppressMessages(library(Seurat))
suppressMessages(library(glue))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))

set.seed(42)

addArchRThreads(threads = 36)
addArchRGenome("hg38")

source("scripts/initialSeuratListMaker.R")

################################################################################
# data directories
################################################################################
cellrangerOutDir <- "../hiv_ln_dogma_data/cellranger_out/"
adtHtoOutDir <- "../hiv_ln_dogma_data/adthto_out/"
amuletOutDir <- "../hiv_ln_dogma_data/amulet_out/"


################################################################################
# generate arrow files
################################################################################
samples <- list.files(cellrangerOutDir)

# omit sample 98 because of clog
samples <- samples[samples != "98"]

fragmentsFns <- paste0(cellrangerOutDir, samples, "/outs/atac_fragments.tsv.gz")

chrNames <- read.table("chrnames/all.chrnames.txt", header = FALSE)
chrNamesDiscard <- chrNames[grepl("(random|EBV|chrUn|chrM|chrHXB2)", chrNames[, 1]), 1]
chrNamesDiscard <- c(chrNamesDiscard, chrNames[grepl("^(KI|GL)", chrNames[, 1]), 1])

arrowfiles <- createArrowFiles(
  inputFiles = fragmentsFns,
  sampleNames = samples,
  outputNames = paste0("origArrow_v2/", samples),
  minTSS = 4,
  minFrags = 1000,
  addTileMat = TRUE,
  excludeChr = chrNamesDiscard,
  subThreading = FALSE,
  force = FALSE,
  addGeneScoreMat = TRUE)

# setting it lower to reduce RAM usage
addArchRThreads(threads = 24)

initProjDir <- "dogma_init_v2"
qcFiltProjDir <- "dogma_lsi_Harmony_v2"
if (!dir.exists(qcFiltProjDir)) {
  proj <- ArchRProject(
    ArrowFiles = arrowfiles,
    outputDirectory = initProjDir,
    copyArrows = TRUE,
    showLogo = FALSE)
  
  # filter using AMULET
  multCells <- lapply(samples, function(x) {
    amuletOut <- unlist(read.table(paste0(amuletOutDir, x, "/MultipletBarcodes_01.txt"), header = FALSE))
    amuletOut <- paste0(x, "#", amuletOut)
    
    return(amuletOut)
  })
  multCells <- unlist(multCells)
  
  singlets <- setdiff(proj$cellNames, multCells)
  projQcFilter <- proj[singlets, ]
  
  ##############################################################################
  # add RNA into ArchRProj
  ##############################################################################
  seRNA <- import10xFeatureMatrix(
    input = paste0(cellrangerOutDir, samples, "/outs/filtered_feature_bc_matrix.h5"),
    names = samples
  )
  
  # intersection of barcodes in atac and rna
  multiBcs <- intersect(projQcFilter$cellNames, colnames(seRNA))
  projQcFilter <- subsetArchRProject(
    cells = multiBcs,
    ArchRProj = projQcFilter,
    outputDirectory = "dogma_subsetMatched_v2",
  )
  
  # extract HIV RNA counts
  hivCounts <- seRNA@assays@data$data[grep("^HIV$", rownames(seRNA)), ]
  
  # this step already adds in Mito and Ribo automatically :)
  projQcFilter <- addGeneExpressionMatrix(
    input = projQcFilter,
    seRNA = seRNA,
    force = TRUE,
    strictMatch = TRUE,
    excludeChr = c(chrNamesDiscard, "chrNA"))
  
  
  tmp <- getCellColData(ArchRProj = projQcFilter, select = c("Sample", "Gex_nGenes", "Gex_nUMI", "Gex_MitoRatio"))
  
  ggplot(as.data.frame(tmp), aes(x = Sample, y = Gex_MitoRatio)) +
    geom_violin() +
    # geom_hline(yintercept = 500) +
    # geom_hline(yintercept = 6000) +
    theme_classic()
  
  
  # qc filter
  idxPass <- BiocGenerics::which(
    projQcFilter$TSSEnrichment >= 8 &
    log10(projQcFilter$nFrags) >= 3.6 &
    projQcFilter$Gex_nGenes > 500 &
    projQcFilter$Gex_nGenes < 6000 &
    projQcFilter$Gex_MitoRatio < 0.1)
  
  cellsPass <- projQcFilter$cellNames[idxPass]
  projQcFilter <- projQcFilter[cellsPass, ]
  
  
  ##############################################################################
  # generate Seurat objects with ADT and HTO
  ##############################################################################
  cellNameDict <- getCellColData(ArchRProj = projQcFilter, select = c("Sample")) %>%
    as.data.frame(.) %>%
    mutate(cellNames = rownames(.))
  
  cellNameDict <- split(cellNameDict, cellNameDict$Sample)
  
  seuList <- lapply(samples, function(x) {
    message(paste0("working on: ", x))
    adt <- createMatFromKallisto(
      sample = x,
      dir = adtHtoOutDir,
      dataMode = "ADT",
      runEmptyDrops = FALSE,
      outRDSPath = paste0("rds/ADT/", x, "_ADT.rds"),
      sampleSuffix = "_adt_output_count")
    
    colnames(adt) <- paste0(x, "#", colnames(adt), "-1")

    hto <- createMatFromKallisto(
      sample = x,
      dir = adtHtoOutDir,
      dataMode = "HTO",
      sampleSuffix = "_hto_output_count")
    
    colnames(hto) <- paste0(x, "#", colnames(hto), "-1")
    

    allowableCbcs <- Reduce(intersect, list(cellNameDict[[x]]$cellNames, colnames(adt), colnames(hto)))
    if (length(allowableCbcs) != length(cellNameDict[[x]]$cellNames)) {
      message(paste0(x, " has unequal cbc between multiome and adt/hto."))
    }
    
    seuObj <- CreateSeuratObject(
      counts = hto[, allowableCbcs],
      assay = "HTO"
    )
    seuObj[["ADT"]] <- CreateAssayObject(counts = adt[, allowableCbcs])
    
    return(seuObj)
  })
  
  seuList <- lapply(seuList, function(x) {
    message(paste0("working on ", Cells(x)[1]))

    x <- NormalizeData(x, assay = "HTO", normalization.method = "CLR", margin = 2)
    x <- MULTIseqDemux(x, assay = "HTO", autoThresh = TRUE, verbose = FALSE)

    return(x)
  })

  seuListSinglet <- lapply(seuList, function(x) {
    y <- subset(x, subset = MULTI_ID != "Doublet" & MULTI_ID != "Negative")
    return(y)
  })

  seuMerged <- merge(
    x = seuListSinglet[[1]],
    y = do.call(c, seuListSinglet[2:length(seuListSinglet)]),
    add.cell.ids = NULL)
  
  seuMerged <- JoinLayers(seuMerged)
  
  # clear memory
  gc()

  # doublet filter archrproj
  projQcFilter <- projQcFilter[names(seuMerged$MULTI_ID), ]
  
  projQcFilter <- addCellColData(projQcFilter, data = seuMerged$MULTI_ID, name = "hashedID", cells = names(seuMerged$MULTI_ID))
  
  # filter based on known antibody mixture per pool (corrects for any
  # erroneous alignment since any cells with these will not map to a unhashed
  # donor)
  tmp <- getCellColData(projQcFilter, select = c("hashedID", "Sample"))
  
  tmp <- as.data.frame(tmp) %>%
    mutate(cbc = rownames(.)) %>%
    mutate(Sample = as.numeric(Sample)) %>%
    mutate(check = case_when(
      Sample >= 1 & Sample <= 8 ~ hashedID %in% c("A0251", "A0252", "A0253"),
      Sample >= 9 & Sample <= 24 ~ hashedID %in% c("A0251", "A0252", "A0253", "A0254", "A0256", "A0257"),
      Sample >= 25 & Sample <= 32 ~ hashedID %in% c("A0251", "A0252"),
      Sample >= 33 & Sample <= 62 ~ hashedID %in% c("A0253", "A0254", "A0256", "A0257", "A0258"),
      Sample >= 63 & Sample <= 80 ~ hashedID %in% c("A0251", "A0252", "A0253"),
      Sample >= 81 & Sample <= 88 ~ hashedID %in% c("A0251", "A0252"),
      Sample >= 89 & Sample <= 96 ~ hashedID %in% c("A0253", "A0254"),
      Sample >= 97 & Sample <= 105 ~ hashedID %in% c("A0256", "A0257", "A0258"),
    )) %>%
    filter(check)
  
  seuMerged$cellName <- names(seuMerged$orig.ident)
  seuMergedFinal <- subset(seuMerged, subset = cellName %in% tmp$cbc)
  
  # save adt/hto seurat object
  saveRDS(seuMergedFinal, "rds/adtHto_matched_final.rds")
  
  projUnhashed <- subsetArchRProject(
    cells = tmp$cbc,
    ArchRProj = projQcFilter,
    outputDirectory = "dogma_unhashed_v2",
  )
  
  rm(tmp)
  gc()
  
  ##############################################################################
  # perform iterative LSI
  ##############################################################################
  # addArchRThreads(threads = 16)
  # 
  # projUnhashed <- addIterativeLSI(
  #   ArchRProj = projUnhashed,
  #   useMatrix = "TileMatrix",
  #   depthCol = "nFrags",
  #   saveIterations = FALSE,
  #   name = "LSI_ATAC",
  #   iterations = 2,
  #   clusterParams = list(
  #     resolution = 0.2, 
  #     sampleCells = 125000,
  #     n.start = 10,
  #     maxClusters = 30
  #   )
  # )
  # 
  # projUnhashed <- addIterativeLSI(
  #   ArchRProj = projUnhashed, 
  #   clusterParams = list(
  #     resolution = 0.2, 
  #     sampleCells = 125000,
  #     n.start = 10,
  #     maxClusters = 30
  #   ),
  #   saveIterations = FALSE,
  #   useMatrix = "GeneExpressionMatrix", 
  #   depthCol = "Gex_nUMI",
  #   varFeatures = 2500,
  #   firstSelection = "variable",
  #   binarize = FALSE,
  #   name = "LSI_RNA"
  # )
  # 
  # ##############################################################################
  # # add information to ArchRProj and metadata
  # ##############################################################################
  hashedMeta <- read.csv("meta/hiv_ln_dogma_completeMeta_hashed.csv")
  unhashedMeta <- read.csv("meta/hiv_ln_dogma_completeMeta_unhashed.csv")

  hashedMeta <- hashedMeta %>%
    mutate(hashedSample = as.character(hashedSample)) %>%
    rename(PoolName = Sample)

  cellDf <- getCellColData(projUnhashed,
    select = c("Sample", "hashedID")) %>%
    as.data.frame(.) %>%
    mutate(cbc = rownames(.))

  rownames(cellDf) <- NULL

  cellDf <- cellDf %>%
    left_join(hashedMeta, by = c("Sample" = "hashedSample")) %>%
    left_join(unhashedMeta, by = join_by("PoolName" == "Pool", "hashedID" == "HTO_antibody"))

  cellDf <- cellDf %>%
    select(-Sort_Day, -Old_number, -Day_of_sample_ID, -Well_number)

  projUnhashed$pool <- cellDf$PoolName
  projUnhashed$sortDate <- cellDf$sort_date
  projUnhashed$donor <- cellDf$participant_ID
  projUnhashed$chip <- cellDf$chip

  # add back HIV RNA count info
  projUnhashed$HIV_nUMI <- 0
  
  # save archr project
  projUnhashed <- saveArchRProject(
    ArchRProj = projUnhashed,
    outputDirectory = qcFiltProjDir,
    load = TRUE
  )
  
} else {
  projUnhashed <- loadArchRProject(
    path = qcFiltProjDir,
    showLogo = FALSE
  )
}