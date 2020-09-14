# restricts to pop and saves/creates mapping
MapCovariates <- function(cohortMethodData, population = NULL, mapping = NULL){
  
  newCovariateData <- Andromeda::andromeda(covariateRef = cohortMethodData$covariateRef,
                                           analysisRef = cohortMethodData$analysisRef)
  
  # restrict to population for speed
  ParallelLogger::logTrace('restricting to population for speed and mapping')
  if(is.null(mapping)){
    mapping <- data.frame(oldCovariateId = as.data.frame(cohortMethodData$covariateRef %>% dplyr::select(covariateId)),
                          newCovariateId = 1:as.data.frame(dplyr::count(cohortMethodData$covariateRef))[[1]])
  }
  if(sum(colnames(mapping)%in%c('oldCovariateId','newCovariateId'))!=2){
    colnames(mapping) <- c('oldCovariateId','newCovariateId')
  }
  cohortMethodData$mapping <- mapping
  if(is.null(population)){
    cohortMethodData$population <- data.frame(rowId = cohortMethodData$cohorts %>% select(rowId)) 
  } else{
    cohortMethodData$population <- data.frame(rowId = population%>%select(rowId)) 
  }
  # assign new ids :
  newCovariateData$covariates <- cohortMethodData$covariates %>%
    dplyr::inner_join(cohortMethodData$population) %>% 
    dplyr::rename(oldCovariateId = covariateId) %>% 
    dplyr::inner_join(cohortMethodData$mapping) %>% 
    dplyr::select(-oldCovariateId)  %>%
    dplyr::rename(covariateId = newCovariateId)
  cohortMethodData$population <- NULL
  cohortMethodData$mapping <- NULL
  
  newCovariateData$mapping <- mapping
  
  return(newCovariateData)
}
 

# restricts to pop and saves/creates mapping
#' @export
MapCovariates2 <- function(covariateData, population = NULL, mapping = NULL){
  
  newCovariateData <- Andromeda::andromeda(covariateRef = covariateData$covariateRef,
                                           analysisRef = covariateData$analysisRef)
  includedCovariateIds <- covariateData$covariates %>% dplyr::pull(covariateId)
  covariateData$covariateRefTidy <- covariateData$covariateRef %>% 
    dplyr::filter(.data$covariateId %in% !!includedCovariateIds)
  attr(newCovariateData, "metaData")$deletedInfrequentCovariateIds <- attr(covariateData, "metaData")$deletedInfrequentCovariateIds
  attr(newCovariateData, "metaData")$deletedRedundantCovariateIds <- attr(covariateData, "metaData")$deletedRedundantCovariateIds
  
  # restrict to population for speed
  ParallelLogger::logTrace('restricting to population for speed and mapping')
  if(is.null(mapping)){
    mapping <- data.frame(oldCovariateId = as.data.frame(covariateData$covariateRefTidy %>% dplyr::select(covariateId)),
                          newCovariateId = 1:as.data.frame(dplyr::count(covariateData$covariateRefTidy))[[1]])
  }
  if(sum(colnames(mapping)%in%c('oldCovariateId','newCovariateId'))!=2){
    colnames(mapping) <- c('oldCovariateId','newCovariateId')
  }
  
  if(is.null(population)){
    rowIds <- covariateData$covariates %>% dplyr::pull(rowId) %>% unique()
  } else{
    rowIds <- population$rowId
  }
  covariateData$mapping <- mapping
  # assign new ids :
  newCovariateData$covariates <- covariateData$covariates %>%
    dplyr::filter(.data$rowId %in% !!rowIds) %>%
    dplyr::rename(oldCovariateId = covariateId) %>% 
    dplyr::inner_join(covariateData$mapping) %>% 
    dplyr::select(-oldCovariateId)  %>%
    dplyr::rename(covariateId = newCovariateId)
  
  newCovariateData$mapping <- mapping
  attr(newCovariateData, "metaData")$covariateRefTidy <- as.data.frame(covariateData$covariateRefTidy)
  
  covariateData$mapping <- NULL
  covariateData$covariateRefTidy <- NULL
  return(newCovariateData)
}
