# Copyright 2020 Observational Health Data Sciences and Informatics
#
# This file is part of CohortMethod
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# @author Observational Health Data Sciences and Informatics
# @author Seng Chan You

#' Create propensity scores
#'
#' @description
#' Creates propensity scores using a regularized logistic regression.
#'
#' @param cohortMethodData         An object of type [CohortMethodData] as generated using
#'                                 [getDbCohortMethodData()].
#' @param population               A data frame describing the population. This should at least have a
#'                                 `rowId` column corresponding to the `rowId` column in the
#'                                 [CohortMethodData] covariates object and a `treatment` column.
#'                                 If population is not specified, the full population in the
#'                                 [CohortMethodData] will be used.
#' @param excludeCovariateIds      Exclude these covariates from the propensity model.
#' @param includeCovariateIds      Include only these covariates in the propensity model.
#' @param maxCohortSizeForFitting  If the target or comparator cohort are larger than this number, they
#'                                 will be downsampled before fitting the propensity model. The model
#'                                 will be used to compute propensity scores for all subjects. The
#'                                 purpose of the sampling is to gain speed. Setting this number to 0
#'                                 means no downsampling will be applied.
#' @param removeRedundancy         If true, the function will remove the redundant covariates by using
#'                                 `FeatureExtraction::tidyCovariateData` function. Default setting is true.
#'
#' @examples
#' data(cohortMethodDataSimulationProfile)
#' cohortMethodData <- simulateCohortMethodData(cohortMethodDataSimulationProfile, n = 1000)
#' ps <- createPs(cohortMethodData)
#'
#' @export
createDataForEncoder <- function(cohortMethodData,
                                 population = NULL,
                                 mapping = NULL,
                                 excludeCovariateIds = c(),
                                 includeCovariateIds = c(),
                                 maxCohortSizeForFitting = 250000,
                                 fileName = NULL,
                                 weight = "default",
                                 tidyCovariate = T){
  # if (exportFolder == getwd()){
  #   warning("We recommend to specify the exportFolder for weights of autoencoder")
  # }
  # if (!file.exists(exportFolder)) {
  #   dir.create(exportFolder, recursive = TRUE)
  # }
  if (is.null(population)) {
    population <- cohortMethodData$cohorts %>%
      collect()
  }
  if (!("rowId" %in% colnames(population)))
    stop("Missing column rowId in population")
  if (!("treatment" %in% colnames(population)))
    stop("Missing column treatment in population")
  
  start <- Sys.time()
  population <- population[order(population$rowId), ]
  if (cohortMethodData$cohorts %>% count() %>% pull() == 0) {
    error <- "No subjects in population, so cannot fit model"
    sampled <- FALSE
    ref <- NULL
  } else if (cohortMethodData$covariates %>% count() %>% pull() == 0) {
    error <- "No covariate data, so cannot fit model"
    sampled <- FALSE
    ref <- NULL
  } else {
    covariates <- cohortMethodData$covariates %>%
      filter(.data$rowId %in% local(population$rowId))
    
    if (length(includeCovariateIds) != 0) {
      covariates <- covariates %>%
        filter(.data$covariateId %in% includeCovariateIds)
    }
    if (length(excludeCovariateIds) != 0) {
      covariates <- covariates %>%
        filter(!.data$covariateId %in% includeCovariateIds)
    }
    filteredCovariateData <- Andromeda::andromeda(covariates = covariates,
                                                  covariateRef = cohortMethodData$covariateRef,
                                                  analysisRef = cohortMethodData$analysisRef)
    metaData <- attr(cohortMethodData, "metaData")
    metaData$populationSize <- nrow(population)
    attr(filteredCovariateData, "metaData") <- metaData
    class(filteredCovariateData) <- "CovariateData"
    if(tidyCovariate){
      covariateData <- FeatureExtraction::tidyCovariateData(filteredCovariateData)
    }else{
      covariateData <- Andromeda::copyAndromeda(filteredCovariateData)
    }
    Andromeda::close(filteredCovariateData)
    on.exit(Andromeda::close(covariateData))
    
    attr(population, "metaData")$deletedInfrequentCovariateIds <- attr(covariateData, "metaData")$deletedInfrequentCovariateIds
    attr(population, "metaData")$deletedRedundantCovariateIds <- attr(covariateData, "metaData")$deletedRedundantCovariateIds
    # covariates <- covariateData$covariates
    
    # sampled <- FALSE
    # if (maxCohortSizeForFitting != 0) {
    #   set.seed(0)
    #   targetRowIds <- population$rowId[population$treatment == 1]
    #   if (length(targetRowIds) > maxCohortSizeForFitting) {
    #     ParallelLogger::logInfo(paste0("Downsampling target cohort from ", length(targetRowIds), " to ", maxCohortSizeForFitting, " before fitting"))
    #     targetRowIds <- sample(targetRowIds, size = maxCohortSizeForFitting, replace = FALSE)
    #     sampled <- TRUE
    #   }
    #   comparatorRowIds <- population$rowId[population$treatment == 0]
    #   if (length(comparatorRowIds) > maxCohortSizeForFitting) {
    #     ParallelLogger::logInfo(paste0("Downsampling comparator cohort from ", length(comparatorRowIds), " to ", maxCohortSizeForFitting, " before fitting"))
    #     comparatorRowIds <- sample(comparatorRowIds, size = maxCohortSizeForFitting, replace = FALSE)
    #     sampled <- TRUE
    #   }
    #   if (sampled) {
    #     fullPopulation <- population
    #     fullCovariates <- covariates
    #     population <- population[population$rowId %in% c(targetRowIds, comparatorRowIds), ]
    #     covariates <- covariates %>%
    #       filter(.data$rowId %in% local(population$rowId))
    #   }
    # }
  }
  
  newcovariateData <- MapCovariates2(covariateData,
                                     population,
                                     mapping=mapping)
  
  ParallelLogger::logDebug(paste0('Max covariateId in covariates: ',as.data.frame(newcovariateData$covariates %>% dplyr::summarise(max = max(covariateId, na.rm=T)))))
  ParallelLogger::logDebug(paste0('# covariates in covariateRef: ', nrow(newcovariateData$covariateRef)))
  ParallelLogger::logDebug(paste0('Max rowId in covariates: ', as.data.frame(newcovariateData$covariates %>% dplyr::summarise(max = max(rowId, na.rm=T)))))
  
  maxY <- as.data.frame(newcovariateData$mapping %>% dplyr::summarise(max=max(newCovariateId, na.rm = TRUE)))$max
  ParallelLogger::logDebug(paste0('Max newCovariateId in mapping: ',maxY))
  maxX <- max(studyPop$rowId)
  ParallelLogger::logDebug(paste0('Max rowId in population: ',maxX))
  
  data <- Matrix::sparseMatrix(i = 1,
                               j = 1,
                               x = 0,
                               dims = c(maxX,maxY))
  convertData <- function(batch) {
    data <<- data + Matrix::sparseMatrix(i = as.data.frame(batch %>% select(rowId))$rowId,
                                         j = as.data.frame(batch %>% select(covariateId))$covariateId,
                                         x = as.data.frame(batch %>% select(covariateValue))$covariateValue,
                                         dims=c(maxX,maxY))
    return(NULL)
  }
  
  Andromeda::batchApply(newcovariateData$covariates, convertData, batchSize = 100000)
  if(!is.null(fileName)){
    Andromeda::saveAndromeda(newcovariateData, fileName = fileName, maintainConnection = TRUE)
  }
  
  return(list(data=data,
              newcovariateData = newcovariateData,
              mapping = newcovariateData$mapping))
}

#' Create propensity scores
#'
#' @description
#' Creates propensity scores using a regularized logistic regression.
#'
#' @param cohortMethodData         An object of type [CohortMethodData] as generated using
#'                                 [getDbCohortMethodData()].
#' @param population               A data frame describing the population. This should at least have a
#'                                 `rowId` column corresponding to the `rowId` column in the
#'                                 [CohortMethodData] covariates object and a `treatment` column.
#'                                 If population is not specified, the full population in the
#'                                 [CohortMethodData] will be used.
#' @param excludeCovariateIds      Exclude these covariates from the propensity model.
#' @param includeCovariateIds      Include only these covariates in the propensity model.
#' @param maxCohortSizeForFitting  If the target or comparator cohort are larger than this number, they
#'                                 will be downsampled before fitting the propensity model. The model
#'                                 will be used to compute propensity scores for all subjects. The
#'                                 purpose of the sampling is to gain speed. Setting this number to 0
#'                                 means no downsampling will be applied.
#' @param removeRedundancy         If true, the function will remove the redundant covariates by using
#'                                 `FeatureExtraction::tidyCovariateData` function. Default setting is true.
#'
#' @examples
#' data(cohortMethodDataSimulationProfile)
#' cohortMethodData <- simulateCohortMethodData(cohortMethodDataSimulationProfile, n = 1000)
#' ps <- createPs(cohortMethodData)
#'
#' @export
# createDataForEncoder <- function(cohortMethodData,
#                                  population = NULL,
#                                  excludeCovariateIds = c(),
#                                  includeCovariateIds = c(),
#                                  maxCohortSizeForFitting = 250000,
#                                  #epochNum = 50,
#                                  #targetDim = 50,
#                                  #weight = "default",
#                                  removeRedundancy = T
# ){
#   if (is.null(population)) {
#     population <- cohortMethodData$cohorts %>%
#       collect()
#   }
#   if (!("rowId" %in% colnames(population)))
#     stop("Missing column rowId in population")
#   if (!("treatment" %in% colnames(population)))
#     stop("Missing column treatment in population")
#   
#   start <- Sys.time()
#   population <- population[order(population$rowId), ]
#   if (cohortMethodData$cohorts %>% count() %>% pull() == 0) {
#     error <- "No subjects in population, so cannot fit model"
#     sampled <- FALSE
#     ref <- NULL
#   } else if (cohortMethodData$covariates %>% count() %>% pull() == 0) {
#     error <- "No covariate data, so cannot fit model"
#     sampled <- FALSE
#     ref <- NULL
#   } else {
#     covariates <- cohortMethodData$covariates %>%
#       filter(.data$rowId %in% local(population$rowId))
#     
#     if (length(includeCovariateIds) != 0) {
#       covariates <- covariates %>%
#         filter(.data$covariateId %in% includeCovariateIds)
#     }
#     if (length(excludeCovariateIds) != 0) {
#       covariates <- covariates %>%
#         filter(!.data$covariateId %in% includeCovariateIds)
#     }
#     filteredCovariateData <- Andromeda::andromeda(covariates = covariates,
#                                                   covariateRef = cohortMethodData$covariateRef,
#                                                   analysisRef = cohortMethodData$analysisRef)
#     metaData <- attr(cohortMethodData, "metaData")
#     metaData$populationSize <- nrow(population)
#     attr(filteredCovariateData, "metaData") <- metaData
#     class(filteredCovariateData) <- "CovariateData"
#     
#     covariateData <- FeatureExtraction::tidyCovariateData(filteredCovariateData, removeRedundancy = removeRedundancy)
#     Andromeda::close(filteredCovariateData)
#     on.exit(Andromeda::close(covariateData))
#     
#     attr(population, "metaData")$deletedInfrequentCovariateIds <- attr(covariateData, "metaData")$deletedInfrequentCovariateIds
#     attr(population, "metaData")$deletedRedundantCovariateIds <- attr(covariateData, "metaData")$deletedRedundantCovariateIds
#     #covariates <- covariateData$covariates
#     
#     #sampled <- FALSE
#     # if (maxCohortSizeForFitting != 0) {
#     #   set.seed(0)
#     #   targetRowIds <- population$rowId[population$treatment == 1]
#     #   if (length(targetRowIds) > maxCohortSizeForFitting) {
#     #     ParallelLogger::logInfo(paste0("Downsampling target cohort from ", length(targetRowIds), " to ", maxCohortSizeForFitting, " before fitting"))
#     #     targetRowIds <- sample(targetRowIds, size = maxCohortSizeForFitting, replace = FALSE)
#     #     sampled <- TRUE
#     #   }
#     #   comparatorRowIds <- population$rowId[population$treatment == 0]
#     #   if (length(comparatorRowIds) > maxCohortSizeForFitting) {
#     #     ParallelLogger::logInfo(paste0("Downsampling comparator cohort from ", length(comparatorRowIds), " to ", maxCohortSizeForFitting, " before fitting"))
#     #     comparatorRowIds <- sample(comparatorRowIds, size = maxCohortSizeForFitting, replace = FALSE)
#     #     sampled <- TRUE
#     #   }
#     #   if (sampled) {
#     #     fullPopulation <- population
#     #     fullCovariates <- covariates
#     #     population <- population[population$rowId %in% c(targetRowIds, comparatorRowIds), ]
#     #     covariates <- covariates %>%
#     #       filter(.data$rowId %in% local(population$rowId))
#     #   }
#     # }
#   }
#   
#   newcovariateData <- MapCovariates2(covariateData,
#                                      population,
#                                      mapping=NULL)
#   
#   ParallelLogger::logDebug(paste0('Max covariateId in covariates: ',as.data.frame(newcovariateData$covariates %>% dplyr::summarise(max = max(covariateId, na.rm=T)))))
#   ParallelLogger::logDebug(paste0('# covariates in covariateRef: ', nrow(newcovariateData$covariateRef)))
#   ParallelLogger::logDebug(paste0('Max rowId in covariates: ', as.data.frame(newcovariateData$covariates %>% dplyr::summarise(max = max(rowId, na.rm=T)))))
#   
#   maxY <- as.data.frame(newcovariateData$mapping %>% dplyr::summarise(max=max(newCovariateId, na.rm = TRUE)))$max
#   ParallelLogger::logDebug(paste0('Max newCovariateId in mapping: ',maxY))
#   maxX <- max(studyPop$rowId)
#   ParallelLogger::logDebug(paste0('Max rowId in population: ',maxX))
#   
#   data <- Matrix::sparseMatrix(i=1,
#                                j=1,
#                                x=0,
#                                dims=c(maxX,maxY))
#   convertData <- function(batch) {
#     data <<- data + Matrix::sparseMatrix(i=as.data.frame(batch %>% select(rowId))$rowId,
#                                          j=as.data.frame(batch %>% select(covariateId))$covariateId,
#                                          x=as.data.frame(batch %>% select(covariateValue))$covariateValue,
#                                          dims=c(maxX,maxY))
#     return(NULL)
#   }
#   Andromeda::batchApply(newcovariateData$covariates, convertData, batchSize = 100000)
#   return(list(data=data,
#               newcovariateData = newcovariateData))
# }
