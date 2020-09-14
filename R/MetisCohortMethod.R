#' Create a parameter object for the function reduceDim
#'
#' @details
#' Create an object defining the parameter values.
#'
#' @param pathToEncoder                Should only the first exposure per subject be included?
#' @param pathToMap           Restrict the analysis to the period when both exposures are
#'                                         observed?
#'
#' @export
createReduceDimArgs <- function(pathToAutoencoderWeight = "",
                                pathToMap = ""
) {
  analysis <- list()
  for (name in names(formals(createReduceDimArgs))) {
    analysis[[name]] <- get(name)
  }
  class(analysis) <- "args"
  return(analysis)
}

#' Create a CohortMethod analysis specification
#'
#' @details
#' Create a set of analysis choices, to be used with the \code{\link{runCmAnalyses}} function.
#'
#' @param analysisId                      An integer that will be used later to refer to this specific
#'                                        set of analysis choices.
#' @param description                     A short description of the analysis.
#' @param targetType                      If more than one target is provided for each
#'                                        drugComparatorOutcome, this field should be used to select
#'                                        the specific target to use in this analysis.
#' @param comparatorType                  If more than one comparator is provided for each
#'                                        drugComparatorOutcome, this field should be used to select
#'                                        the specific comparator to use in this analysis.
#' @param getDbCohortMethodDataArgs       An object representing the arguments to be used when calling
#'                                        the \code{\link{getDbCohortMethodData}} function.
#' @param createStudyPopArgs              An object representing the arguments to be used when calling
#'                                        the \code{\link{createStudyPopulation}} function.
#' @param createPs                        Should the \code{\link{createPs}} function be used in this
#'                                        analysis?
#' @param createPsArgs                    An object representing the arguments to be used when calling
#'                                        the \code{\link{createPs}} function.
#' @param trimByPs                        Should the \code{\link{trimByPs}} function be used in this
#'                                        analysis?
#' @param trimByPsArgs                    An object representing the arguments to be used when calling
#'                                        the \code{\link{trimByPs}} function.
#' @param trimByPsToEquipoise             Should the \code{\link{trimByPsToEquipoise}} function be used
#'                                        in this analysis?
#' @param trimByPsToEquipoiseArgs         An object representing the arguments to be used when calling
#'                                        the \code{\link{trimByPsToEquipoise}} function.
#' @param matchOnPs                       Should the \code{\link{matchOnPs}} function be used in this
#'                                        analysis?
#' @param matchOnPsArgs                   An object representing the arguments to be used when calling
#'                                        the \code{\link{matchOnPs}} function.
#' @param matchOnPsAndCovariates          Should the \code{\link{matchOnPsAndCovariates}} function be
#'                                        used in this analysis?
#' @param matchOnPsAndCovariatesArgs      An object representing the arguments to be used when calling
#'                                        the \code{\link{matchOnPsAndCovariates}} function.
#' @param stratifyByPs                    Should the \code{\link{stratifyByPs}} function be used in
#'                                        this analysis?
#' @param stratifyByPsArgs                An object representing the arguments to be used when calling
#'                                        the \code{\link{stratifyByPs}} function.
#' @param stratifyByPsAndCovariates       Should the \code{\link{stratifyByPsAndCovariates}} function
#'                                        be used in this analysis?
#' @param stratifyByPsAndCovariatesArgs   An object representing the arguments to be used when calling
#'                                        the \code{\link{stratifyByPsAndCovariates}} function.
#' @param fitOutcomeModel                 Should the \code{\link{fitOutcomeModel}} function be used in
#'                                        this analysis?
#' @param fitOutcomeModelArgs             An object representing the arguments to be used when calling
#'                                        the \code{\link{fitOutcomeModel}} function.
#'
#' @export
createCmAnalysis <- function(analysisId = 1,
                             description = "",
                             targetType = NULL,
                             comparatorType = NULL,
                             getDbCohortMethodDataArgs,
                             createStudyPopArgs,
                             createPs = FALSE,
                             createPsArgs = NULL,
                             trimByPs = FALSE,
                             trimByPsArgs = NULL,
                             trimByPsToEquipoise = FALSE,
                             trimByPsToEquipoiseArgs = NULL,
                             matchOnPs = FALSE,
                             matchOnPsArgs = NULL,
                             matchOnPsAndCovariates = FALSE,
                             matchOnPsAndCovariatesArgs = NULL,
                             stratifyByPs = FALSE,
                             stratifyByPsArgs = NULL,
                             stratifyByPsAndCovariates = FALSE,
                             stratifyByPsAndCovariatesArgs = NULL,
                             fitOutcomeModel = FALSE,
                             fitOutcomeModelArgs = NULL,
                             reduceDim = FALSE,
                             reduceDimArgs = NULL) {
  if (matchOnPs + matchOnPsAndCovariates + stratifyByPs + stratifyByPsAndCovariates > 1) {
    stop("Need to pick one matching or stratification function")
  }
  if (trimByPs && trimByPsToEquipoise) {
    stop("Cannot trim to fraction and equipoise at the same time")
  }
  if (!createPs && (trimByPs | matchOnPs | matchOnPsAndCovariates | stratifyByPs | stratifyByPsAndCovariates)) {
    stop("Must create propensity score model to use it for trimming, matching, or stratification")
  }
  if (!(matchOnPs | matchOnPsAndCovariates | stratifyByPs | stratifyByPsAndCovariates) && !is.null(fitOutcomeModelArgs) &&
      fitOutcomeModelArgs$stratified) {
    stop("Must create strata by using matching or stratification to fit a stratified outcome model")
  }
  if (!createPs) {
    createPsArgs <- NULL
  }
  if (!trimByPs) {
    trimByPsArgs <- NULL
  }
  if (!trimByPsToEquipoise) {
    trimByPsToEquipoiseArgs <- NULL
  }
  if (!matchOnPs) {
    matchOnPsArgs <- NULL
  }
  if (!matchOnPsAndCovariates) {
    matchOnPsAndCovariatesArgs <- NULL
  }
  if (!stratifyByPs) {
    stratifyByPsArgs <- NULL
  }
  if (!stratifyByPsAndCovariates) {
    stratifyByPsAndCovariatesArgs <- NULL
  }
  if (!fitOutcomeModel) {
    fitOutcomeModelArgs <- NULL
  }
  if (!reduceDim) {
    reduceDimArgs <- NULL
  }
  
  # First: get the default values:
  analysis <- list()
  for (name in names(formals(createCmAnalysis))) {
    analysis[[name]] <- get(name)
  }
  
  # Next: overwrite defaults with actual values if specified:
  values <- lapply(as.list(match.call())[-1], function(x) eval(x, envir = sys.frame(-3)))
  for (name in names(values)) {
    if (name %in% names(analysis)) {
      analysis[[name]] <- values[[name]]
    }
  }
  
  class(analysis) <- "cmAnalysis"
  return(analysis)
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
#' @param errorOnHighCorrelation   If true, the function will test each covariate for correlation with
#'                                 the treatment assignment. If any covariate has an unusually high
#'                                 correlation (either positive or negative), this will throw and
#'                                 error.
#' @param stopOnError              If an error occurr, should the function stop? Else, the two cohorts
#'                                 will be assumed to be perfectly separable.
#' @param prior                    The prior used to fit the model. See
#'                                 [Cyclops::createPrior()] for details.
#' @param control                  The control object used to control the cross-validation used to
#'                                 determine the hyperparameters of the prior (if applicable). See
#'                                 [Cyclops::createControl()] for details.
#' @param removeRedundancy         If true, the function will remove the redundant covariates by using
#'                                 `FeatureExtraction::tidyCovariateData` function. Default setting is true.
#'
#' @examples
#' data(cohortMethodDataSimulationProfile)
#' cohortMethodData <- simulateCohortMethodData(cohortMethodDataSimulationProfile, n = 1000)
#' ps <- createPs(cohortMethodData)
#'
#' @export
createEncodedPs <- function(cohortMethodData,
                            encoder = NULL,
                            population = NULL,
                            excludeCovariateIds = c(),
                            includeCovariateIds = c(),
                            maxCohortSizeForFitting = 250000,
                            errorOnHighCorrelation = TRUE,
                            stopOnError = TRUE,
                            prior = createPrior("laplace", exclude = c(0), useCrossValidation = TRUE),
                            control = createControl(noiseLevel = "silent",
                                                    cvType = "auto",
                                                    seed = 1,
                                                    tolerance = 2e-07,
                                                    cvRepetitions = 10,
                                                    startingVariance = 0.01),
                            removeRedundancy = T,
                            cmEncoder = NULL) {
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
    
    if(is.null(cmEncoder)){
      covariateData <- FeatureExtraction::tidyCovariateData(filteredCovariateData, removeRedundancy = removeRedundancy)
      close(filteredCovariateData)
      on.exit(close(covariateData))
      covariates <- covariateData$covariates
      attr(population, "metaData")$deletedInfrequentCovariateIds <- attr(covariateData, "metaData")$deletedInfrequentCovariateIds
      attr(population, "metaData")$deletedRedundantCovariateIds <- attr(covariateData, "metaData")$deletedRedundantCovariateIds
    }else{
      map = cmEncoder$map
      
    }
    
    sampled <- FALSE
    if (maxCohortSizeForFitting != 0) {
      set.seed(0)
      targetRowIds <- population$rowId[population$treatment == 1]
      if (length(targetRowIds) > maxCohortSizeForFitting) {
        ParallelLogger::logInfo(paste0("Downsampling target cohort from ", length(targetRowIds), " to ", maxCohortSizeForFitting, " before fitting"))
        targetRowIds <- sample(targetRowIds, size = maxCohortSizeForFitting, replace = FALSE)
        sampled <- TRUE
      }
      comparatorRowIds <- population$rowId[population$treatment == 0]
      if (length(comparatorRowIds) > maxCohortSizeForFitting) {
        ParallelLogger::logInfo(paste0("Downsampling comparator cohort from ", length(comparatorRowIds), " to ", maxCohortSizeForFitting, " before fitting"))
        comparatorRowIds <- sample(comparatorRowIds, size = maxCohortSizeForFitting, replace = FALSE)
        sampled <- TRUE
      }
      if (sampled) {
        fullPopulation <- population
        fullCovariates <- covariates
        population <- population[population$rowId %in% c(targetRowIds, comparatorRowIds), ]
        covariates <- covariates %>%
          filter(.data$rowId %in% local(population$rowId))
      }
    }
    population <- population[order(population$rowId), ]
    outcomes <- population
    colnames(outcomes)[colnames(outcomes) == "treatment"] <- "y"
    covariateData$outcomes <- outcomes
    floatingPoint <- getOption("floatingPoint")
    if (is.null(floatingPoint)) {
      floatingPoint <- 64
    } else {
      ParallelLogger::logInfo("Cyclops using precision of ", floatingPoint)
    }
    cyclopsData <- Cyclops::convertToCyclopsData(covariateData$outcomes, covariates, modelType = "lr", quiet = TRUE, floatingPoint = floatingPoint)
    error <- NULL
    ref <- NULL
    if (errorOnHighCorrelation) {
      suspect <- Cyclops::getUnivariableCorrelation(cyclopsData, threshold = 0.5)
      suspect <- suspect[!is.na(suspect)]
      if (length(suspect) != 0) {
        covariateIds <- as.numeric(names(suspect))
        ref <- cohortMethodData$covariateRef %>%
          filter(.data$covariateId %in% covariateIds) %>%
          collect()
        ParallelLogger::logInfo("High correlation between covariate(s) and treatment detected:")
        ParallelLogger::logInfo(paste(colnames(ref), collapse = "\t"))
        for (i in 1:nrow(ref))
          ParallelLogger::logInfo(paste(ref[i, ], collapse = "\t"))
        message <- "High correlation between covariate(s) and treatment detected. Perhaps you forgot to exclude part of the exposure definition from the covariates?"
        if (stopOnError) {
          stop(message)
        } else {
          error <- message
        }
      }
    }
  }
  if (is.null(error)) {
    cyclopsFit <- tryCatch({
      Cyclops::fitCyclopsModel(cyclopsData, prior = prior, control = control)
    }, error = function(e) {
      e$message
    })
    if (is.character(cyclopsFit)) {
      if (stopOnError) {
        stop(cyclopsFit)
      } else {
        error <- cyclopsFit
      }
    } else if (cyclopsFit$return_flag != "SUCCESS") {
      if (stopOnError) {
        stop(cyclopsFit$return_flag)
      } else {
        error <- cyclopsFit$return_flag
      }
    }
  }
  if (is.null(error)) {
    error <- "OK"
    cfs <- coef(cyclopsFit)
    if (all(cfs[2:length(cfs)] == 0)) {
      warning("All coefficients (except maybe the intercept) are zero. Either the covariates are completely uninformative or completely predictive of the treatment. Did you remember to exclude the treatment variables from the covariates?")
    }
    if (sampled) {
      # Adjust intercept to non-sampled population:
      y.bar <- mean(population$treatment)
      y.odds <- y.bar/(1 - y.bar)
      y.bar.new <- mean(fullPopulation$treatment)
      y.odds.new <- y.bar.new/(1 - y.bar.new)
      delta <- log(y.odds) - log(y.odds.new)
      cfs[1] <- cfs[1] - delta  # Equation (7) in King and Zeng (2001)
      cyclopsFit$estimation$estimate[1] <- cfs[1]
      covariateData$fullOutcomes <- fullPopulation
      population <- fullPopulation
      population$propensityScore <- predict(cyclopsFit, newOutcomes = covariateData$fullOutcomes, newCovariates = fullCovariates)
    } else {
      population$propensityScore <- predict(cyclopsFit)
    }
    attr(population, "metaData")$psModelCoef <- coef(cyclopsFit)
    attr(population, "metaData")$psModelPriorVariance <- cyclopsFit$variance[1]
  } else {
    if (sampled) {
      population <- fullPopulation
    }
    population$propensityScore <- population$treatment
    attr(population, "metaData")$psError <- error
    if (!is.null(ref)) {
      attr(population, "metaData")$psHighCorrelation <- ref
    }
  }
  population <- computePreferenceScore(population)
  delta <- Sys.time() - start
  ParallelLogger::logDebug("Propensity model fitting finished with status ", error)
  ParallelLogger::logInfo("Creating propensity scores took ", signif(delta, 3), " ", attr(delta, "units"))
  return(population)
}
