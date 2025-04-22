createMatFromKallisto <- function(
  sample,
  dir,
  dataMode,
  emptyDropsLower = 100,
  outRDSPath = paste0("rds/", dataMode, "/", sample, ".rds"),
  sampleSuffix = "_count_out",
  runEmptyDrops = TRUE,
  kallistoPrefix = "output")
{
  
  stopifnot(dataMode %in% c("ADT", "HTO"))
  
  if (file.exists(outRDSPath)) {
    mat <- readRDS(outRDSPath)
    
    return(mat)
  }
  
  bus_output <- read.table(file = glue("{dir}/{sample}{sampleSuffix}/{kallistoPrefix}.mtx"), header = FALSE, skip = 4)
  feats <- read.table(file = glue("{dir}/{sample}{sampleSuffix}/{kallistoPrefix}.genes.txt"), header = FALSE)
  cbcs <- read.table(file = glue("{dir}/{sample}{sampleSuffix}/{kallistoPrefix}.barcodes.txt"), header = FALSE)
  
  count_matrix <- Matrix::sparseMatrix(
    i = bus_output$V2,
    j = bus_output$V1,
    x = bus_output$V3)
  
  rownames(count_matrix) <- feats$V1
  colnames(count_matrix) <- cbcs$V1
  
  if (dataMode == "ADT") {
    
    if (runEmptyDrops) {
      empty_calc <- emptyDrops(count_matrix, lower = emptyDropsLower)
      adt_filt <- count_matrix[, which(empty_calc$FDR < 0.05)]
      
      saveRDS(adt_filt, file = outRDSPath)
      return(adt_filt)
      
    } else {
      saveRDS(count_matrix, file = outRDSPath)
      return(count_matrix)
    }
    
  } else {
    return(count_matrix)
  }
}

