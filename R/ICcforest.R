#' Fit a conditional inference forest for interval-censored survival data
#'
#' An implementation of the random forest and bagging ensemble algorithms utilizing conditional 
#' inference trees as base learners for interval-censored survival data.
#' 
#' \code{ICcforest} returns an \code{ICcforest} object.
#' The object belongs to the class \code{ICcforest}, as a subclass of \code{\link[partykit]{cforest}}.
#' This function extends the conditional inference survival forest algorithm in 
#' \code{\link[partykit]{cforest}} to fit interval-censored survival data. 
#' 
#' 
#' @param formula a formula object, with the response being a \code{\link[survival]{Surv}} object, 
#' with form
#' 
#' 
#' \code{Surv(time1, time2, type="interval2")}.
#' @param data a data frame containing the variables named in \code{formula}.
#' @param na.action a function which indicates what should happen when the data contain 
#' missing values.
#' @param mtry number of input variables randomly sampled as candidates at each node for 
#' random forest like algorithms. The default \code{mtry} is tuned by \code{\link{tuneICRF}}.
#' @param ntree an integer, the number of the trees to grow for the forest. \code{ntree = 100L} is set 
#' by default.
#' @param perturb a list with arguments \code{replace} and \code{fraction} determining which 
#' type of resampling, with \code{replace = TRUE} referring to the n-out-of-n 
#' bootstrap and \code{replace = FALSE} referring to sample splitting. \code{fraction} is 
#' the proportion of observations to draw without replacement.
#' @param applyfun an optional \code{lapply}-style function with arguments 
#' \code{function(X, FUN, ...)}. 
#' It is used for computing the variable selection criterion. The default is to use the
#' basic \code{lapply} function unless the \code{cores} argument is specified (see below). 
#' See \code{\link[partykit]{ctree_control}}.
#' @param cores numeric. If set to an integer the \code{applyfun} is set to \code{\link[parallel]{mclapply}} 
#' with the desired number of cores. See \code{\link[partykit]{ctree_control}}.
#' @param trace whether to print the progress of the search of the optimal value of \code{mtry}
#' when \code{mtry} is not specified (see \code{\link{tuneICRF}}). \code{trace = TRUE} is set by default.
#' @param control a list of control parameters, see \code{\link[partykit]{ctree_control}}.
#' \code{control} parameters \code{minsplit}, \code{minbucket} have been adjusted from the 
#' \code{\link[partykit]{cforest}} defaults. Other default values correspond to those of the 
#' default values used by \code{\link[partykit]{ctree_control}}. 
#' @param suppress a logical specifying whether the messages from \code{\link[icenReg]{getFitEsts}} 
#' are suppressed. If \code{FALSE}, the messages are printed. \code{suppress = TRUE} is set by default.
#' @param ... additional arguments.
#' @keywords Ensemble method, conditional inference forest, interval-censored data
#' @return An object of class \code{ICcforest}, as a subclass of \code{\link[partykit]{cforest}}.
#' @import partykit
#' @import survival 
#' @import stats
#' @import utils
#' @import icenReg
#' @seealso \code{\link{predict.ICcforest}} for prediction, \code{\link{gettree.ICcforest}}
#' for individual tree extraction, and \code{\link{tuneICRF}} for \code{mtry} tuning. 
#' @examples
#' #### Example with miceData
#' library(icenReg)
#' data(miceData)
#' 
#' ## For ICcforest to run, Inf should be set to be a large number, for example, 9999999.
#' miceData$u[miceData$u == Inf] <- 9999999.
#' 
#' ## Fit an iterval-censored conditional inference forest
#' Cforest <- ICcforest(Surv(l, u, type = "interval2") ~ grp, data = miceData)
#' 
#' @export

ICcforest <- function(formula, data, mtry = NULL, ntree = 100L, applyfun = NULL, 
                      cores = NULL, na.action = na.pass, suppress = TRUE, trace = TRUE,
                      perturb = list(replace = FALSE, fraction = 0.632), 
                      control = partykit::ctree_control(teststat = "quad", testtype = "Univ", 
                                                        mincriterion = 0, saveinfo = FALSE,
                                                        minsplit = nrow(data) * 0.15, 
                                                        minbucket = nrow(data) * 0.06), ...){
  
  # package version dependency
  if (packageVersion("partykit") < "1.2.2") {
    stop("partykit >= 1.2.2 needed for this function.", call. = FALSE)
  }
  
  if (packageVersion("icenReg") < "2.0.8") {
    stop("icenReg >= 2.0.8 needed for this function.", call. = FALSE)
  }
  
  requireNamespace("inum")
  if (!requireNamespace("partykit", quietly = TRUE)) {
    stop("Package \"pkg\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  
  X <- data[,as.character(formula[[2]][[2]])]
  Y <- data[,as.character(formula[[2]][[3]])]
  
  if(sum(X==Y)){ ## add small noise to observed event case
    epsilon <- min(diff(sort(unique(c(X,Y)))))/20
    ID <- which(X==Y)
    data[ID,as.character(formula[[2]][[3]])] <- epsilon + data[ID,as.character(formula[[2]][[3]])]
  }

  ## x2 is Surv(Left,right,type="interval2") object
  .logrank_trafo <- function(x2){
    if(!(is.Surv(x2) && isTRUE(attr(x2, "type") == "interval"))){
      stop("Response must be a 'Survival' object with Surv(time1,time2,event) format")
    }
    
    # Fit IC survival curve
    Curve <- ic_np(x2[, 1:2])
    # get estimated survival
    Left <- 1 - getFitEsts(Curve, q = x2[, 1])
    Right <- 1 - getFitEsts(Curve, q = x2[, 2])
    
    Log_Left <- ifelse(Left<=0,0,Left*log(Left))
    Log_Right<- ifelse(Right<=0,0,Right*log(Right))
    result <- (Log_Left-Log_Right)/(Left-Right)
    
    return(matrix(as.double(result),ncol=1))
  }
  
  h2 <-function(y, x, start = NULL, weights, offset, estfun = TRUE, object = FALSE, ...) {
    if (all(is.na(weights))==1) weights <- rep(1, NROW(y))
    s <- .logrank_trafo(y[weights > 0,,drop = FALSE])
    r <- rep(0, length(weights))
    r[weights > 0] <- s
    list(estfun = matrix(as.double(r), ncol = 1), converged = TRUE)
  }
  
  if (is.null(mtry)) 
    mtry <- tuneICRF(formula, data, control = control, suppress = suppress, trace = trace)
  
  if (suppress == TRUE){
    invisible(capture.output(res <- partykit::cforest(formula, data, control = control, na.action = na.action, ytrafo = h2, 
                                                      mtry = mtry, ntree = ntree, applyfun = applyfun, cores = cores,
                                                      perturb = perturb, ...)))
  } else {
    res <- partykit::cforest(formula, data, control = control, na.action = na.action, ytrafo = h2, 
                             mtry = mtry, ntree = ntree, applyfun = applyfun, cores = cores,
                             perturb = perturb, ...)
  }
    
  
  # define a new class as a subclass of cforest
  class(res) <- c("ICcforest", class(res))
  return(res)
}
