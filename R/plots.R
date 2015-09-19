## plot qq plot for uniform data

qqunif <- function(y, ylim, main="Uniform Q-Q Plot", 
                   xlab="Theoretical -log10(Pvalue)", ylab="Sample -log10(Pvalue)",
                   plot.it=TRUE, datax=FALSE, ...) {
  y <- sort(y)
  ylogp <- -log10(y)
  x <- (1:length(y))/length(y)
  xlogp <- -log10(x)
  if (plot.it) {
    if (missing(ylim)) ylim <- range(ylogp)
    if (datax) {
      plot(ylogp,xlogp, main=main, xlab=ylab, ylab=xlab, xlim=ylim, ...)
    } else {
      plot(xlogp, ylogp, main=main, xlab=xlab, ylab=ylab, ylim=ylim, ...)
    }
  }
  #invisible(if (datax) list(x = y, y = x) else list(x = x, 
  #                                                  y = y))
}

#' plot the q-q plot of the pvals of the permuted results
#' @param object of class GSTTDataSet with the \code{permuted.results} slot filled
#' @export
qqpermuted <- function(object, pindex="all", plot.it=T) {
  pvals <- c()
  if (pindex=="all") {
    pvals <- unlist(sapply(object@permuted.results, function(x) {unlist(x$pvalue)}))
  } else {
    if (!is.numeric(pindex)) stop("pindex has to be 'all' or a numeric vector of permutation indices")
    pvals <- inlist(sapply(pindex, function(i) {unlist(object@permuted.results[[i]]$pvalue)}))
  }
  qqunif(pvals, main="Q-Q plot of p-values", plot.it=plot.it)
}