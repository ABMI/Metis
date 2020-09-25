#' @export
andromedaToSparseM <- function(andromedaCovariate, fileName = NULL, batchSize = 10000, append = F){
  
  if(is.null(fileName)){
    fileName = "data.mtx.gz"
  }
  if(!grepl("mtx\\.gz$",fileName)) stop("the name of file should end with mtx.gz")
  
  maxY <- as.data.frame(andromedaCovariate %>% dplyr::summarise(max=max(covariateId, na.rm = TRUE)))$max
  ParallelLogger::logDebug(paste0('Max newCovariateId in mapping: ',maxY))
  
  maxX <- as.data.frame(andromedaCovariate %>% dplyr::summarise(max=max(rowId, na.rm = TRUE)))$max
  ParallelLogger::logDebug(paste0('Max rowId in population: ',maxX))
  
  rowNum = as.data.frame(andromedaCovariate %>% dplyr::count())[[1]]
  ParallelLogger::logDebug(paste0('Total count in covariates: ',rowNum))
  
  if(!append){
    startWritingMMgz(fileName = fileName, maxX=maxX, maxY=maxY, rowNum, mType = "real", removeIfFileExist = TRUE)
  }
  
  Andromeda::batchApply(andromedaCovariate, writeMMgz, fileName = fileName, batchSize = batchSize)
  
  covMat <- readBigMM(fileName)
  
  ParallelLogger::logDebug(paste0('The  object size of sparse matrix is: ', round(object.size(covMat)/(2^(30)),3) , "GB"))
  
  covMat <- as(covMat, "dgCMatrix") # To decrease the size of the object
  
  ParallelLogger::logDebug(paste0('The  object size of sparse matrix in dgCMatrix format is: ', round(object.size(covMat)/(2^(30)),3) , "GB"))
}

#' @export
mapCov <- function(covariateData, mapping = NULL, deletedCovariateIds = c()){
  ParallelLogger::logTrace('mapping Covariate ID to new covariate IDs')
  newCovariateData <- Andromeda::andromeda(covariateRef = covariateData$covariateRef,
                                           analysisRef = covariateData$analysisRef)
  
  # includedCovariateIds <- covariateData$covariates %>% dplyr::pull(covariateId)
  # covariateData$covariateRefTidy <- covariateData$covariateRef %>% 
  #   dplyr::filter(.data$covariateId %in% !!includedCovariateIds)
  # attr(newCovariateData, "metaData")$deletedInfrequentCovariateIds <- attr(covariateData, "metaData")$deletedInfrequentCovariateIds
  # attr(newCovariateData, "metaData")$deletedRedundantCovariateIds <- attr(covariateData, "metaData")$deletedRedundantCovariateIds
  
  # restrict to population for speed
  if(is.null(mapping)){
    ParallelLogger::logTrace('generate new mapping')
    oldCovariateId = as.data.frame(covariateData$covariateRef %>% 
                                     dplyr::select(.data$covariateId) %>%
                                     dplyr::filter(!(.data$covariateId %in% !!deletedCovariateIds)) %>%
                                     dplyr::arrange(covariateId)
    ) 
    newCovariateId = 1:nrow(oldCovariateId)
    
    
    mapping <- data.frame(oldCovariateId = oldCovariateId,
                          newCovariateId = newCovariateId)
    colnames(mapping) <- c('oldCovariateId','newCovariateId')
  }
  if(sum(colnames(mapping)%in%c('oldCovariateId','newCovariateId'))!=2){
    colnames(mapping) <- c('oldCovariateId','newCovariateId')
  }
  
  #rowIds <- covariateData$covariates %>% dplyr::pull(rowId) %>% unique()
  
  covariateData$mapping <- mapping
  # assign new ids :
  newCovariateData$covariates <- covariateData$covariates %>%
    #dplyr::filter(.data$rowId %in% !!rowIds) %>%
    dplyr::rename(oldCovariateId = covariateId) %>% 
    dplyr::inner_join(covariateData$mapping) %>% 
    dplyr::select(-oldCovariateId)  %>%
    dplyr::rename(covariateId = newCovariateId)
  
  newCovariateData$mapping <- mapping
  
  return(newCovariateData)
}

#' @export
reduceCovariateData <- function(covariateData, rowIdLimit){
  covariates <- covariateData$covariates %>%
    filter(.data$rowId <= rowIdLimit)
  
  filteredCovariateData <- Andromeda::andromeda(covariates = covariates,
                                                covariateRef = covariateData$covariateRef,
                                                analysisRef = covariateData$analysisRef)
  metaData <- attr(covariateData, "metaData")
  metaData$rowIdLimit <- rowIdLimit
  attr(filteredCovariateData, "metaData") <- metaData
  class(filteredCovariateData) <- "CovariateData"
  
  return(filteredCovariateData)
}