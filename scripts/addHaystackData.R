addHaystackData <- function(
  proj,
  haystackParentDir,
  haystackSamples,
  sample_suffix = "_hxb2"
) {
  
  haystackList <- lapply(haystackSamples, function(i) {
    fn <- glue("{haystackParentDir}/{i}{sample_suffix}/viralFrags.tsv")
    message(paste0("working on ", fn))
    
    df <- read.csv(fn, header = TRUE, sep = "\t")
    df <- df %>%
      filter(alreadyRecordedInIntegration != "True")
    
    if (nrow(df) == 0) {
      message("No viral frags found in this sample. Be careful that haystackSamples is named!!")
      return(NULL)
    }
    
    intSiteFragsFn <- glue("{haystackParentDir}/{i}{sample_suffix}/integrationSites_viralFrags.tsv")
    if (file.exists(intSiteFragsFn)) {
      df2 <- read.csv(intSiteFragsFn, header = TRUE, sep = "\t")
      
      if (nrow(df2) > 0) {
        df <- bind_rows(df, df2)
      }
      
    }
    return(df)
  })
  
  names(haystackList) <- haystackSamples
  
  haystackDf <- bind_rows(haystackList, .id = "sample")
  haystackDf <- haystackDf %>%
    mutate(newCbc = paste0(sample, "#", cbc))
  
  # count UMI too (0.5 = one mate detecte)
  haystackPosCells <- haystackDf %>%
    group_by(sample, newCbc) %>%
    dplyr::summarize(totalUniqFrags = n() / 2)
  
  matchedHaystackArchRCells <- intersect(haystackPosCells$newCbc, proj$cellNames)
  
  haystackPosCellsViralFrags <- haystackDf %>%
    filter(newCbc %in% matchedHaystackArchRCells)
  
  umiCountFilt <- haystackPosCells %>%
    filter(newCbc %in% matchedHaystackArchRCells)
  
  haystackUmiCount <- rep(0, length(proj$cellNames))
  names(haystackUmiCount) <- proj$cellNames
  haystackUmiCount[umiCountFilt$newCbc] <- umiCountFilt$totalUniqFrags
  
  proj$haystackUmiCount <- haystackUmiCount
  proj$haystackOut <- haystackUmiCount > 0
  
  return(list(newProj = proj, allViralFrags = haystackDf, filteredviralFrags = haystackPosCellsViralFrags))
}