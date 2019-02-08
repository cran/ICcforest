#' Extract an individual tree from an ICcforest object
#'
#' Extract the i-th individual tree from the established ICcforest. The resulting object can be
#' printed or plotted, and predictions can be made using it.
#' 
#' @param object an object as returned by \code{\link{ICcforest}}.
#' @param tree an integer, the number of the tree to extract from the forest.
#' @param ... additional arguments.
#' @return An object of class \code{\link[partykit]{party}}.
#' @import partykit
#' @import icenReg
#' @importFrom survival Surv
#' @import stats
#' @import utils
#' @examples
#' 
#' #### Example with dataset miceData
#' library(icenReg)
#' data(miceData)
#' 
#' ## For ICcforest to run, Inf should be set to be a large number, for example, 9999999.
#' idx_inf <- (miceData$u == Inf)
#' miceData$u[idx_inf] <- 9999999.
#' 
#' ## First, fit an iterval-censored conditional inference forest
#' Cforest <- ICcforest(formula = Surv(l,u,type="interval2")~grp, data = miceData, ntree = 50L)
#' ## Extract the 50-th tree from the forest
#' plot(gettree(Cforest, tree = 50L))
#' 
#' @export
#'
#'

gettree.ICcforest <- function(object, tree = 1L, ...) {
  ft <- object$fitted
  ft[["(weights)"]] <- object$weights[[tree]]
  ret <- partykit::party(object$nodes[[tree]], data = object$data, fitted = ft)
  ret$terms <- object$terms
  class(ret) <- c("constparty", class(ret))
  return(ret)
}

