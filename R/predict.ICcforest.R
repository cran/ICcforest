#' Predict from an ICcforest model
#' 
#' Compute predictions from ICcforest objects.
#' 
#' @param object an object as returned by \code{\link{ICcforest}}.
#' @param newdata an optional data frame containing test data.
#' @param type a character string denoting the type of predicted value returned. 
#' 
#' 
#' For \code{"type = response"}, the mean of a numeric response, the median survival time for 
#' the interval-censored response is returned. For \code{"type = prob"}, a list with the survival 
#' function constructed using the non-parametric maximum likelihood estimator for each observation 
#' is returned. \code{"type = weights"} returns an integer vector of prediction weights. 
#' For \code{type = "node"}, a list of terminal node ids for each of the trees in the 
#' forest is returned.
#' @param FUN a function to compute summary statistics. Predictions for each node must be 
#' computed based on arguments \code{(y, w)} where \code{y} is the response and \code{w} are case 
#' weights. 
#' @param OOB a logical specifying whether out-of-bag predictions are desired
#' 
#' 
#' (only if \code{newdata = NULL}).
#' @param simplify a logical indicating whether the resulting list of predictions should be 
#' converted to a suitable vector or matrix (if possible), see \code{\link[partykit]{cforest}}.
#' @param scale a logical indicating scaling of the nearest neighbor weights by the sum of weights 
#' in the corresponding terminal node of each tree, see \code{\link[partykit]{cforest}}.
#' @param suppress a logical specifying whether the messages from \code{\link[icenReg]{getFitEsts}} 
#' are suppressed. If \code{FALSE}, the messages are printed. \code{suppress = TRUE} is set by default.
#' @param ... additional arguments.
#' @return An object of class \code{ICcforest}, as a subclass of \code{\link[partykit]{cforest}}.
#' @import partykit
#' @import survival 
#' @import stats
#' @import utils
#' @import icenReg
#' @seealso \code{\link{sbrier_IC}} for evaluation of model fit for interval-censored data
#' @examples 
#' library(icenReg)
#' data(miceData)
#' 
#' ## For ICcforest to run, Inf should be set to be a large number, for example, 9999999.
#' miceData$u[miceData$u == Inf] <- 9999999.
#' 
#' ## First, fit an iterval-censored conditional inference forest
#' Cforest <- ICcforest(formula = Surv(l,u,type="interval2")~grp, data = miceData)
#' ## Predict the survival function constructed using the non-parametric maximum likelihood estimator
#' Pred <- predict(Cforest, type = "prob")
#' 
#' ## Out-of-bag prediction of the median survival time
#' PredOOB <- predict(Cforest, type = "response", OOB = TRUE)
#' 
#' @export
#' 
#' 

predict.ICcforest <- function(object, newdata = NULL, OOB = FALSE, suppress = TRUE,
                              type = c("response", "prob", "weights", "node"),
                              FUN = NULL, simplify = TRUE, scale = TRUE, ...){
  
  # package version dependency
  if (packageVersion("partykit") < "1.2.2") {
    stop("partykit >= 1.2.2 needed for this function.", call. = FALSE)
  }
  
  if (packageVersion("icenReg") < "2.0.8") {
    stop("icenReg >= 2.0.8 needed for this function.", call. = FALSE)
  }
  
  if (missing(type)){
    stop("The value of 'type' is not specified.")
  }
  
  pred_Surv_prob <- function(y, w) {
    if (length(y) == 0) return(NA)
    idx = which(w>0)
    y = y[idx]
    w = w[idx]
    yy = as.matrix(y[,1:2])
    ic_np(yy, weights = w)
  }
  
  pred_Surv_response <- function(y, w) {
    if (length(y) == 0) return(NA)
    getFitEsts(pred_Surv_prob(y, w))
  }
  
  if (is.null(FUN)) 
    Fun <- if (type == "response") pred_Surv_response else pred_Surv_prob
  else
    Fun <- NULL
  
  if (suppress == TRUE){
    invisible(capture.output(res <- partykit::predict.cforest(object = object, newdata = newdata, type = type, 
                                                              OOB = OOB, simplify = simplify, scale = scale,
                                                              FUN = Fun, ...)))
  } else{
    res <- partykit::predict.cforest(object = object, newdata = newdata, type = type, 
                                     OOB = OOB, simplify = simplify, scale = scale,
                                     FUN = Fun, ...)
  }
  return(res)
}


