getFitEstsN <- function(fit, newdata = NULL, p = NULL, q = NULL){
  invisible(capture.output(res <- getFitEsts(fit, newdata, p, q)))
  return(res)
}

# get survival probability at time teval for ic_np, ic_sp, ic_par or ic_bayes object
getFitEsts_Surv <- function(Curve, newData = NULL, teval){
  if (is.null(newData)){
    Shat <- 1-getFitEstsN(fit = Curve, q = teval)
  } else {
    Shat <- 1-getFitEstsN(fit = Curve, newdata = newData, q = teval)
  }
  return(Shat)
}

#' Model Fit For Interval-Censored Data
#'
#' Compute the (integrated) Brier score to evaluate the model fit for interval-censored survival data.
#'
#' @param obj an object of class \code{\link[survival]{Surv}}.
#' @param pred predicted values. This can be a matrix of survival probabilities evaluated
#' at a sequence of time points for a set of new data, a list of \code{\link[survival]{survfit}}
#' objects, a list \code{\link[icenReg]{ic_np}}
#' objects, or a list of \code{\link[icenReg]{ic_sp}} objects.
#' @param btime a vector of length two indicating the range of times that the scores are computed on.
#' The default \code{btime} is set to be the vector of the smallest and the largest values among
#' all left and right endpoints given in \code{obj}.
#' @param type a character string denoting the type of scores returned. For \code{"IBS"},
#' the integrated Brier score over the \code{btime} is returned. For \code{"BS"}, the
#' Brier score at every left and right endpoint of all censoring intervals that lie within
#' \code{btime} is returned.
#' @keywords Brier score, integrated Brier score
#' @return
#' If \code{type = "IBS"}, this returns the integrated Brier score.
#' @return
#' If \code{type = "BS"}, this returns the Brier scores.
#' @import partykit
#' @import icenReg
#' @importFrom survival Surv
#' @import stats
#' @import ipred
#' @import utils
#' @references S. Tsouprou. Measures of discrimination and predictive accuracy for interval-censored
#' data. Master thesis, Leiden University. https://www.math.leidenuniv.nl/scripties/MasterTsouprou.pdf.
#' @examples
#' ### Example with dataset miceData
#' library(survival)
#' library(icenReg)
#' data(miceData)
#'
#' ## For proper evaluation, Inf should be set to be a large number, for example, 9999999.
#' idx_inf <- (miceData$u == Inf)
#' miceData$u[idx_inf] <- 9999999.
#'
#' obj <- Surv(miceData$l, miceData$u, type = "interval2")
#'
#' ## Model fit for an NPMLE survival curve with survfit
#' pred <- survival::survfit(formula = Surv(l, u, type = "interval2") ~ 1, data = miceData)
#' # Integrated Brier score up to time = 642
#' sbrier_IC(obj, pred, btime = c(0, 642), type = "IBS")
#'
#' ## Model fit for a semi-parametric model with icenReg::ic_sp()
#' pred <- icenReg::ic_sp(formula = Surv(l, u, type = "interval2") ~ 1, data = miceData)
#' # Integrated Brier score up to the largest endpoints of all censoring intervals in the dataset
#' sbrier_IC(obj, pred, type = "IBS")
#'
#' ## Model fit for an NPMLE survival curve with icenReg::ic_np()
#' pred <- icenReg::ic_np(miceData[,c('l', 'u')])
#' # Brier score computed at every left and right endpoints of all censoring intervals in the dataset
#' sbrier_IC(obj, pred, type = "BS")
#'
#'
#'
#' @export

sbrier_IC <- function(obj, pred, btime = range(as.numeric(obj[,1:2])), type = c("IBS","BS")) {
  if (missing(type))
    stop("The value of 'type' is not specified.")

  if (!inherits(obj, "Surv"))
    stop("obj is not of class Surv")
  class(obj) <- NULL
  # number of obs
  N <- nrow(obj)

  intL <- obj[, 1]
  intR <- obj[, 2]
  # the observation time in test set
  time <- c(intL,intR)

  if (is.null(btime))
    stop("btime not given")
  if (length(btime) < 1)
    stop("btime not given")
  if (length(btime) != 2) {
    stop("btime should be a vector of two indicating the range of the times")
  } else {
    if (btime[1] < min(time))
      warning("btime[1] is smaller than min(time)")
    if (btime[2] > max(time))
      warning("btime[2] is larger than max(time)")
    # all the observation times in the range of obs times
    btime <- unique(time[time >= btime[1] & time <= btime[2]])
    btime <- sort(btime)
  }

  ptype <- class(pred)[[1]]
  if (is.null(ptype)) {
    if (is.vector(pred))
      ptype <- "vector"
    if (is.list(pred))
      ptype <- "list"
  }

  if (ptype == "numeric" && is.vector(pred))
    ptype <- "vector"

  survs <- NULL
  switch(ptype, survfit = {
    survs <- try(ipred::getsurv(pred, time), silent=T)
    if(class(survs) == "try-error"){
      stop("please check whether the package 'ipred' is installed")
    }
    survs <- matrix(rep(survs, N), nrow = length(time))
  }, ic_np = {

    survs <- try(getFitEsts_Surv(Curve = pred, teval = time), silent=T)
    if(class(survs) == "try-error"){
      stop("please check whether the package 'icenReg' is installed")
    }
    survs <- matrix(rep(survs, N), nrow = length(time))
  }, ic_ph = {
    survs <- try(getFitEsts_Surv(Curve = pred, teval = time), silent=T)
    if(class(survs) == "try-error"){
      stop("please check whether the package 'icenReg' is installed")
    }
    survs <- matrix(rep(survs, N), nrow = length(time))
  }, list = {
    if (!inherits(pred[[1]], "survfit") && !inherits(pred[[1]], "ic_ph") && !inherits(pred[[1]], "ic_np"))
      stop("pred is not a list of survfit/ic_ph/ic_np objects; \n  if pred is an ICcforest prediction object, please make sure it is created using the setting type = 'prob'")
    if (length(pred) != N) stop("pred must be of length(time)")
    if (inherits(pred[[1]], "survfit")){
      M.list <- try(lapply(pred, ipred::getsurv, times = time),silent=T)
      if (class(M.list) == "try-error"){
        stop("please check whether the package 'ipred' is installed")
      }
      survs <- matrix(unlist(lapply(M.list, function(x) x)),
                      nrow = length(time), ncol = N) ## each obs in one column
    } else if (inherits(pred[[1]], "ic_np") || (inherits(pred[[1]], "ic_ph"))){
      M.list <- try(lapply(pred, getFitEsts_Surv, teval = time), silent=T)
      if(class(M.list) == "try-error"){
        stop("please check whether the package 'icenReg' is installed")
      }
      survs <- matrix(unlist(lapply(M.list, function(x) x)),
                      nrow = length(time), ncol = N) ## each obs in one column
    }
  }, vector = {
    if (length(pred) != N) stop("pred must be of length(time)")
    if (length(time) != 1) stop("cannot compute integrated Brier score with pred; \n  if pred is an ICcforest prediction object, please make sure it is created using the setting type = 'prob'")
    survs <- pred
  }, matrix = {
    if (all(dim(pred) == c(length(time), N))) survs <- pred
    else stop("wrong dimensions of pred; \n  if pred is an ICcforest prediction object, please make sure it is created using the setting type = 'prob'")
  })
  if (is.null(survs))
    stop("unknown type of pred")

  bsc <- rep(0, length(time))

  for (i in 1:length(time)) {
    bsc_temp <- matrix(0, nrow = N, ncol = 1)
    k = 0
    for (j in 1:N){
      St <- survs[i,j]
      if (time[i] <= intL[j]){
        IY = 1
      } else if(time[i] > intR[j]){
        IY = 0
      } else {
        Sr <- survs[N+j, j]
        Sl <- survs[j, j]
        if (Sl != Sr){
          IY = (St - Sr)/(Sl - Sr)
        } else {
          IY = NULL
        }
      }
      if (is.null(IY)==0){
        k = k + 1
        bsc_temp[k] = (IY - St)^2
      }
    }
    bsc[i] = mean(bsc_temp[1:k])
  }

  # get bsc and time ordered
  ot = order(time)
  bsc = bsc[ot]
  time = time[ot]

  # get bsc and time unique
  unik <- !duplicated(time)
  bsc <- bsc[unik]
  time <- time[unik]

  # bsc now corresponds with btime
  bsc <- bsc[time<=max(btime)]

  if (type == "IBS"){
    # apply trapezoid rule
    idx <- 2:length(btime)
    RET <- diff(btime) %*% ((bsc[idx - 1] + bsc[idx])/2)
    # the largest is much larger than the second largest
    if (max(btime)/tail(btime, n=2)[1] > 1000 && tail(bsc, n=1)==0 ){
      RET <- RET - (tail(btime, n=1) - tail(btime, n=2)[1])*((tail(bsc, n=1) + tail(bsc, n=2)[1])/2)
      btime <- head(btime, n = -1)
    }
    RET <- RET/diff(range(btime))
    names(RET) <- "integrated Brier score"
    attr(RET, "time") <- range(btime)
  } else if (type == "BS"){
    RET <- bsc
    names(RET) <- btime
    attr(RET, "type") <- "Brier score"
  } else {
    stop("unknown type of results")
  }
  RET
}
