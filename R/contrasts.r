#' Get the numerical contrast vector for the dataset
#'
#' @param gstt an object of class DESeqDataSet 
#' @seealso \code{DESeq2:::cleanContrast}
#' @export
#' @examples
#' library(DESeq2)
#' example(DESeq)
#' getContrast(dds)
getContrast <- function(object, contrast, name) {
  
  ## the following code is copied from DESeq2 package
  if (!"results" %in% mcols(mcols(object))$type) {
    stop("results should have been computed.")
  }
  
  test <- attr(object,"test")
  
  isExpanded <- attr(object, "modelMatrixType") == "expanded"
  termsOrder <- attr(terms.formula(design(object)),"order")
  # allows use of 'name' for expanded model matrices if there are interactions
  if ((test == "Wald") & isExpanded & missing(contrast) & all(termsOrder < 2)) {
    if (missing(name)) {
      designVars <- all.vars(design(object))
      lastVarName <- designVars[length(designVars)]
      lastVar <- colData(object)[[lastVarName]]
      if (is.factor(lastVar)) {
        nlvls <- nlevels(lastVar)
        contrast <- c(lastVarName, levels(lastVar)[nlvls], levels(lastVar)[1])
        cont <- getGroupsFromCharContrast(object, contrast, isExpanded)
        return(cont)
      }
    } else {     
      stop("\n
              note: an expanded model matrix was used in fitting the model.
              
              recommendation: the 'contrast' argument should be used to extract
              log2 fold changes of levels against each other.\n")
    }
  }
  if (missing(name)) {
    name <- DESeq2:::lastCoefName(object)
  }
  altHypothesis <- match.arg(altHypothesis, choices=c("greaterAbs","lessAbs","greater","less"))
  stopifnot(lfcThreshold >= 0)
  stopifnot(length(lfcThreshold)==1)
  stopifnot(length(altHypothesis)==1)
  stopifnot(length(alpha)==1)
  stopifnot(length(pAdjustMethod)==1)
  if (length(name) != 1 | !is.character(name)) {
    stop("the argument 'name' should be a character vector of length 1")
  }
  if (lfcThreshold == 0 & altHypothesis == "lessAbs") {
    stop("when testing altHypothesis='lessAbs', set the argument lfcThreshold to a positive value")
  }
  
  # check to see at least one of these are present
  WaldResults <- paste0("WaldPvalue_",name) %in% names(mcols(object))
  LRTResults <- "LRTPvalue" %in% names(mcols(object))
  if (! ( WaldResults | LRTResults) ) {
    stop("cannot find appropriate results, for available names call 'resultsNames(object)'")
  }
  
  # if performing a contrast call the function cleanContrast()
  if (!missing(contrast)) {
    if (is.character(contrast) & length(contrast) != 3) {
      stop("contrast should either be a character vector of length 3, of the form:
contrast = c('factorName','numeratorLevel','denominatorLevel'),
or a numeric vector the same length as resultsNames(dds).
see the manual page of ?results for more information")
    }
    
    # pass down whether the model matrix type was "expanded"
    cont <- getGroupsFromCharContrast(object, contrast, expanded=isExpanded)
    return(cont)
  } else {
    stop("You've uncovered a corner case! Please contact the maintainer 
         with your dataset and reproducible steps for help.")
  }
  stop("You've uncovered a corner case! Please contact the maintainer 
         with your dataset and reproducible steps for help.")
}


#' Get the list of two groups of samples for the contrast.
#' The groups of samples can then be used to get the balanced permutations.
#' @param object of class DESeqDataset
#' @param contrast a vector c('contrastFactor','Level1','Level2') if we want do Level1 vs Level2
#' 
getGroupsFromCharContrast <- function(object, contrast, expanded) {
  if ( missing(contrast) | !is.character(contrast) | (length(contrast)!=3) ) {
    stop("should be called with a character contrast vector c('factor','level1','level2')")
  }
  # check if the appropriate columns are in the resultsNames
  if (contrast[2] == contrast[3]) {
    stop(paste(contrast[2],"and",contrast[3],"should be different level names"))
  }
  contrastFactor <- contrast[1]
  colData <- colData(object)
  groups <- list()
  groups[[1]] <- which(colData[,contrastFactor]==contrast[2])
  groups[[2]] <- which(colData[,contrastFactor]==contrast[3])
  return(groups)
}
