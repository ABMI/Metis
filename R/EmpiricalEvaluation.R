#' @export
EmpiricalEvaluation <- function(cohortMethodData,
                                ps,
                                outcomeIds,
                                hoiIds,
                                ncsIds,
                                analysisId,
                                plot = T,
                                cmOutputFolder= NULL){
  
  if (!file.exists(cmOutputFolder)) {
    dir.create(cmOutputFolder, recursive = TRUE)
  }
  analysisSum <- data.frame()
  for(outcomeId in outcomeIds){
    ####Study Population#####
    studyPop <- CohortMethod::createStudyPopulation(cohortMethodData = cohortMethodData,
                                                    outcomeId = outcomeId,
                                                    firstExposureOnly = FALSE,
                                                    restrictToCommonPeriod = FALSE,
                                                    washoutPeriod = 0,
                                                    removeDuplicateSubjects = "keep first",
                                                    removeSubjectsWithPriorOutcome = FALSE,
                                                    minDaysAtRisk = 1,
                                                    riskWindowStart = 0,
                                                    startAnchor = "cohort start",
                                                    riskWindowEnd = 30,
                                                    endAnchor = "cohort end")
    #CohortMethod::getAttritionTable(studyPop)
    matchedPop <- CohortMethod::matchOnPs(ps, caliper = 0.2, caliperScale = "standardized logit", maxRatio = 1)
    matchedPop<- matchedPop %>% 
      dplyr::select(rowId, propensityScore, preferenceScore,stratumId) %>%
      dplyr::inner_join(studyPop, by = "rowId")
    if(outcomeId %in% hoiIds){
      if(plot){
        balance <- CohortMethod::computeCovariateBalance(matchedPop, cohortMethodData)
        CohortMethod::plotCovariateBalanceScatterPlot(balance,
                                                      showCovariateCountLabel = TRUE, showMaxLabel = TRUE,
                                                      fileName = file.path(cmOutputFolder, sprintf("balance_a%d_o_%d.png", analysisId, outcomeId)))
        CohortMethod::plotCovariateBalanceOfTopVariables(balance,
                                                         fileName = file.path(cmOutputFolder, sprintf("top_balance_a%d_o_%d.png", analysisId, outcomeId)))
      }
    }
    outcomeModel <- CohortMethod::fitOutcomeModel(population = matchedPop,
                                                  modelType = "cox")
    analysisResult <-as.data.frame(outcomeModel$outcomeModelTreatmentEstimate)
    if(nrow(analysisResult)==1) {
      analysisResult$analysisId <- analysisId
      analysisResult$outcomeId <- outcomeId
      analysisSum <- rbind(analysisSum, analysisResult)
    }
  }
  
  negCons <- analysisSum [analysisSum$analysisId==analysisId & (analysisSum$outcomeId %in% ncsIds),]
  hoi <- analysisSum [analysisSum$analysisId==analysisId & (analysisSum$outcomeId %in% hoiIds),]
  null <- EmpiricalCalibration::fitNull(negCons$logRr, negCons$seLogRr)
  if(plot){
    EmpiricalCalibration::plotCalibrationEffect(negCons$logRr, negCons$seLogRr, hoi$logRr, hoi$seLogRr, null, 
                                                fileName = file.path(cmOutputFolder, sprintf("empirical_cal_a%d.png", analysisId)))
  }
  
  return(analysisSum)
}
