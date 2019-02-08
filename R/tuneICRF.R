#' Tune mtry to the optimal value with respect to out-of-bag error for an ICcforest model
#' 
#' Starting with the default value of mtry, search for the optimal value 
#' (with respect to Out-of-Bag error estimate) of mtry for ICcforest.
#' 
#' @param formula a formula object, with the response being a \code{\link[survival]{Surv}} 
#' object, with form 
#' 
#' 
#' \code{Surv(time1, time2, type="interval2")}.
#' @param data a data frame containing the variables named in \code{Formula}.
#' @param mtryStart starting value of \code{mtry}; default is \code{sqrt(nvar)}.
#' @param stepFactor at each iteration, \code{mtry} is inflated (or deflated) by this value.
#' @param control a list with control parameters, see \code{\link[partykit]{cforest}}. 
#' The default values correspond to those of the default values used by \code{\link{ICcforest}}.
#' @param ntreeTry number of trees used at the tuning step.
#' @param trace whether to print the progress of the search. \code{trace = TRUE} is set by default.
#' @param plot whether to plot the out-of-bag error as a function of \code{mtry}.
#' @param doBest whether to run an ICcforest using the optimal mtry found.
#' @param suppress a logical specifying whether the messages from \code{\link[icenReg]{getFitEsts}} 
#' are suppressed. If \code{FALSE}, the messages are printed. \code{suppress = TRUE} is set by default.
#' @param ... additional arguments.
#' @keywords mtry, out-of-bag errors, brier score
#' @return 
#' If \code{doBest=FALSE} (default), this returns the optimal mtry value of those searched.
#' @return 
#' If \code{doBest=TRUE}, this returns the ICcforest object produced with the optimal mtry.
#' @import partykit
#' @importFrom survival Surv
#' @importFrom graphics axis
#' @import stats
#' @import utils
#' @import icenReg
#' @seealso \code{\link{sbrier_IC}} for evaluation of model fit for interval-censored data 
#' when searching for the optimal value of \code{mtry}.
#' @examples
#' ### Example with dataset tandmob2
#' library(icenReg)
#' data(miceData)
#' 
#' ## For ICcforest to run, Inf should be set to be a large number, for example, 9999999.
#' miceData$u[miceData$u == Inf] <- 9999999.
#'
#' ## Create a new variable to be selected from
#' miceData$new = rep(1:4)
#' 
#' ## Tune mtry 
#' mtryTune <- tuneICRF(Surv(l, u, type = "interval2") ~ grp + new, data = miceData)
#' 
#' @export

tuneICRF <- function(formula, data, mtryStart = NULL, stepFactor = 1.5, ntreeTry = 100L, 
                     control = partykit::ctree_control(teststat = "quad", testtype = "Univ", 
                                                       mincriterion = 0, saveinfo = FALSE,
                                                       minsplit = nrow(data) * 0.15, 
                                                       minbucket = nrow(data) * 0.06), 
                     suppress = TRUE, trace = TRUE, plot = FALSE, doBest = FALSE) {
  # package version dependency
  if (packageVersion("partykit") < "1.2.2") {
    stop("partykit >= 1.2.2 needed for this function.", call. = FALSE)
  }
  
  if (packageVersion("icenReg") < "2.0.8") {
    stop("icenReg >= 2.0.8 needed for this function.", call. = FALSE)
  }
  
  # number of the variables
  allX <- substring(formula,1)[[3]]
  nameX <- strsplit(gsub("\\+", '', allX)," ")[[1]]
  nameX <- nameX[nameX!=""]
  nvar = length(nchar(nameX))
  if (is.null(mtryStart)){
    mtryStart = ceiling(sqrt(nvar))
  } 
  
  # integrated Brier score of out-of-bag samples for a mtry value at test 
  errorOOB_mtry <- function(formula, data, mtryTest, ntreeTry, control, suppress){
    cfOOB <- ICcforest(formula, data, mtry = mtryTest, ntree = ntreeTry, control = control, suppress = suppress)
    predOOB <- predict(object = cfOOB, newdata = NULL, OOB = TRUE, type = "prob", suppress = suppress)
    aSurv <- data[,as.character(formula[[2]][[2]])]
    bSurv <- data[,as.character(formula[[2]][[3]])]
    
    testObj <- survival::Surv(aSurv,bSurv,type = "interval2")
    errorOOB <- unname(sbrier_IC(testObj, predOOB, type = "IBS")[1])
    
    rm(cfOOB)
    rm(predOOB)
    return(errorOOB)
  }
  
  
  # errorOld
  errorOld <- errorOOB_mtry(formula, data, mtryTest=mtryStart, ntreeTry, control, suppress = suppress)

  if (errorOld < 0) stop("Initial setting gave 0 error and no room for improvement.")
  if (trace) {
    cat("mtry =", mtryStart, " OOB Brier score =",
        errorOld, "\n")
  }
  
  oobError <- list()
  oobError[[1]] <- errorOld
  names(oobError)[1] <- mtryStart  
  
  for (direction in c("left", "right")) {
    if (trace) cat("Searching", direction, "...\n")
    mtryCur <- mtryStart
    while (mtryCur != nvar) {
      mtryOld <- mtryCur
      mtryCur <- if (direction == "left") {
        max(1, ceiling(mtryCur / stepFactor))
      } else {
        min(nvar, floor(mtryCur * stepFactor))
      }
      if (mtryCur == mtryOld) break
      
      errorCur <- errorOOB_mtry(formula, data, mtryTest = mtryCur, ntreeTry, control, suppress = suppress)

      if (trace) {
        cat("mtry =",mtryCur, "\tOOB error =",errorCur, "\n")
      }
      oobError[[as.character(mtryCur)]] <- errorCur
      errorOld <- errorCur
    }
  }
  mtry <- sort(as.numeric(names(oobError)))
  res_all <- unlist(oobError[as.character(mtry)])
  res_all <- cbind(mtry = mtry, OOBError = res_all)
  res <- res_all[which.min(res_all[,2]), 1]
  
  if (plot) {
    res = res_all
    plot(res_all, xlab = expression(m[try]), ylab = "OOB Error", type = "o", log = "x", xaxt = "n")
    axis(1, at=res_all[,"mtry"])
  }
  
  if (doBest) 
    res <- ICcforest(formula, data, mtry = res, ntree = ntreeTry, control = control, suppress = suppress)
  
  return(res)
}
