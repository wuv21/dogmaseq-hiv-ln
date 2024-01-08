suppressMessages(library(ArchR))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(rtracklayer))
suppressMessages(library(Seurat))
suppressMessages(library(glue))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))

set.seed(42)

addArchRThreads(threads = 48)
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
fragmentsFns <- paste0(cellrangerOutDir, samples, "/outs/atac_fragments.tsv.gz")

chrNames <- read.table("chrnames/all.chrnames.txt", header = FALSE)
chrNamesDiscard <- chrNames[grepl("(random|EBV|chrUn|chrM|chrHXB2)", chrNames[, 1]), 1]
chrNamesDiscard <- c(chrNamesDiscard, chrNames[grepl("^(KI|GL)", chrNames[, 1]), 1])

arrowfiles <- createArrowFiles(
  inputFiles = fragmentsFns,
  sampleNames = samples,
  outputNames = paste0("origArrow/", samples),
  minTSS = 4,
  minFrags = 1000,
  addTileMat = TRUE,
  excludeChr = chrNamesDiscard,
  subThreading = TRUE,
  force = FALSE,
  addGeneScoreMat = TRUE)

# setting it lower to reduce RAM usage
addArchRThreads(threads = 24)

initProjDir <- "dogma_init"
qcFiltProjDir <- "dogma_lsi"
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
  
  # intersection barcodes in atac and rna
  multiBcs <- intersect(projQcFilter$cellNames, colnames(seRNA))
  projQcFilter <- subsetArchRProject(
    cells = multiBcs,
    ArchRProj = projQcFilter,
    outputDirectory = "dogma_subsetMatched",
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
    )) %>%
    filter(check)
  
  projUnhashed <- subsetArchRProject(
    cells = tmp$cbc,
    ArchRProj = projQcFilter,
    outputDirectory = "dogma_unhashed",
  )
  
  rm(tmp)
  gc()
  
  ##############################################################################
  # perform iterative LSI
  ##############################################################################
  projUnhashed <- addIterativeLSI(
    ArchRProj = projUnhashed,
    useMatrix = "TileMatrix",
    depthCol = "nFrags",
    saveIterations = FALSE,
    name = "LSI_ATAC",
    iterations = 2,
    clusterParams = list(
      resolution = 0.2, 
      sampleCells = 125000,
      n.start = 10,
      maxClusters = 30
    )
  )
  
  projUnhashed <- addIterativeLSI(
    ArchRProj = projUnhashed, 
    clusterParams = list(
      resolution = 0.2, 
      sampleCells = 125000,
      n.start = 10,
      maxClusters = 30
    ),
    saveIterations = FALSE,
    useMatrix = "GeneExpressionMatrix", 
    depthCol = "Gex_nUMI",
    varFeatures = 2500,
    firstSelection = "variable",
    binarize = FALSE,
    name = "LSI_RNA"
  )
  
  ##############################################################################
  # add information to ArchRProj and metadata
  ##############################################################################
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
  projUnhashed$sortDate <- cellDf$date_of_sort
  projUnhashed$donor <- cellDf$participant_ID
  projUnhashed$chip <- cellDf$chip

  # harmonyParams <- list()
  # harmonyParams$data_mat <- getReducedDims(
  #   ArchRProj = projUnhashed, 
  #   reducedDims = "LSI_RNA", 
  #   dimsToUse = NULL, 
  #   scaleDims = NULL, 
  #   corCutOff = 0.75
  # )
  # harmonyParams$verbose <- TRUE
  # harmonyParams$meta_data <- data.frame(getCellColData(
  #   ArchRProj = projUnhashed, 
  #   select = "chip")[rownames(harmonyParams$data_mat), , drop = FALSE])
  # harmonyParams$do_pca <- FALSE
  # harmonyParams$vars_use <- "donor"
  # harmonyParams$plot_convergence <- FALSE
  # 
  # harmonyRes <- harmony::RunHarmony(
  #   data_mat = harmonyParams$data_mat,
  #   meta_data = harmonyParams$meta_data,
  #   vars_use = "donor",
  #   verbose = TRUE
  # )
  # 
  # harmonyParams$meta_data <- data.frame(getCellColData(
  #   ArchRProj = projUnhashed, 
  #   select = "Sample")[rownames(harmonyParams$data_mat), , drop = FALSE])
  # 
  # tmp <- Reduce(rbind, lapply("Sample", function(var_use) {
  #   res <- Matrix::sparse.model.matrix(~0 + as.factor(harmonyParams$meta_data[[var_use]]))
  #   Matrix::t(res)
  # }))
  
  projUnhashed <- addHarmony(
    ArchRProj = projUnhashed,
    reducedDims = "LSI_ATAC",
    name = "Harmony_LSI_ATAC",
    groupBy = c("donor", "chip")
  )
  
  projUnhashed <- addHarmony(
    ArchRProj = projUnhashed,
    reducedDims = "LSI_RNA",
    name = "Harmony_LSI_RNA",
    groupBy = c("donor", "chip")
  )
  
  # Combined Dims
  projUnhashed <- addCombinedDims(projUnhashed, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")
  projUnhashed <- addCombinedDims(projUnhashed, reducedDims = c("Harmony_LSI_ATAC", "Harmony_LSI_RNA"), name =  "Harmony_LSI_Combined")
  
  # UMAP
  projUnhashed <- addUMAP(projUnhashed, reducedDims = "Harmony_LSI_ATAC", name = "Harmony_UMAP_ATAC", minDist = 0.8, force = TRUE)
  projUnhashed <- addUMAP(projUnhashed, reducedDims = "Harmony_LSI_RNA", name = "Harmony_UMAP_RNA", minDist = 0.8, force = TRUE)
  projUnhashed <- addUMAP(projUnhashed, reducedDims = "Harmony_LSI_Combined", name = "Harmony_UMAP_Combined", minDist = 0.8, force = TRUE)
  
  # TODO Initial clustering
  projUnhashed <- addClusters(
    input = projUnhashed,
    reducedDims = "Harmony_LSI_Combined",
    name = "Clusters",
    resolution = 0.3,
    force = TRUE
  )
  
  # TODO add back HIV RNA count info
  hivCounts
  
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