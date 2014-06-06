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
        cont <- getNumericFromCharContrast(object, contrast, isExpanded)
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
    name <- lastCoefName(object)
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
    cont <- getNumericFromCharContrast(object, contrast, expanded=isExpanded)
    return(cont)
  } else {
    stop("You've uncovered a corner case! Please contact the maintainer 
         with your dataset and reproducible steps for help.")
  }
  stop("You've uncovered a corner case! Please contact the maintainer 
         with your dataset and reproducible steps for help.")
}


getNumericFromCharContrast <- function(object, contrast, expanded) {
  if ( missing(contrast) | !is.character(contrast) | (length(contrast)!=3) ) {
    stop("should be called with a character contrast vector c('factor','level1','level2')")
  }
  # check if the appropriate columns are in the resultsNames
  if (contrast[2] == contrast[3]) {
    stop(paste(contrast[2],"and",contrast[3],"should be different level names"))
  }
  contrastFactor <- contrast[1]
  colData <- colData(object)
  cont <- rep(0, length(nrow(colData)))
  cont(which(colData[,contrastFactor]==contrast[2])) <- 1
  return(cont)
  if (!contrastFactor %in% names(colData(object))) {
    stop(paste(contrastFactor,"should be the name of a factor in the colData of the DESeqDataSet"))
  }
  contrastNumLevel <- contrast[2]
  contrastDenomLevel <- contrast[3]
  
  # case 1: standard model matrices: build the appropriate contrast
  # coefficients names are of the form  "factor_level_vs_baselevel"
  # output: contrastNumColumn and contrastDenomColumn
  if (!expanded) {
    
    # then we have a base level for the factor
    contrastBaseLevel <- levels(colData(object)[,contrastFactor])[1]
    
    # use make.names() so the column names are
    # the same as created by DataFrame in mcols(object).
    contrastNumColumn <- make.names(paste0(contrastFactor,"_",contrastNumLevel,"_vs_",contrastBaseLevel))
    contrastDenomColumn <- make.names(paste0(contrastFactor,"_",contrastDenomLevel,"_vs_",contrastBaseLevel))
    resNames <- resultsNames(object)
    
    # check in case the desired contrast is already
    # available in mcols(object), and then we can either
    # take it directly or multiply the log fold
    # changes and stat by -1
    stop("You've uncovered a corner case! Please contact the maintainer 
         with your dataset and reproducible steps for help.")
    if ( contrastDenomLevel == contrastBaseLevel ) {
      # the results can be pulled directly from mcols(object)
      name <- make.names(paste0(contrastFactor,"_",contrastNumLevel,"_vs_",contrastDenomLevel))
      if (!name %in% resNames) {
        stop(paste("as",contrastDenomLevel,"is the base level, was expecting",name,"to be present in 'resultsNames(object)'"))
      }
      test <- "Wald"
      log2FoldChange <- getCoef(object, name)
      lfcSE <- getCoefSE(object, name)
      stat <- getStat(object, test, name)
      pvalue <- getPvalue(object, test, name)
      res <- cbind(mcols(object)["baseMean"],
                   log2FoldChange,lfcSE,stat,
                   pvalue)
      names(res) <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue")
      return(res)
    } else if ( contrastNumLevel == contrastBaseLevel ) {
      # fetch the results for denom vs num 
      # and mutiply the log fold change and stat by -1
      cleanName <- paste(contrastFactor,contrastNumLevel,"vs",contrastDenomLevel)
      swapName <- make.names(paste0(contrastFactor,"_",contrastDenomLevel,"_vs_",contrastNumLevel))
      if (!swapName %in% resNames) {
        stop(paste("as",contrastNumLevel,"is the base level, was expecting",swapName,"to be present in 'resultsNames(object)'"))
      }
      test <- "Wald"
      log2FoldChange <- getCoef(object, swapName)
      lfcSE <- getCoefSE(object, swapName)
      stat <- getStat(object, test, swapName)
      pvalue <- getPvalue(object, test, swapName)
      res <- cbind(mcols(object)["baseMean"],
                   log2FoldChange,lfcSE,stat,
                   pvalue)
      names(res) <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue")
      res$log2FoldChange <- -1 * res$log2FoldChange
      res$stat <- -1 * res$stat
      # also need to swap the name in the mcols
      contrastDescriptions <- paste(c("log2 fold change (MAP):",
                                      "standard error:",
                                      "Wald statistic:",
                                      "Wald test p-value:"),
                                    cleanName)
      mcols(res)$description[mcols(res)$type == "results"] <- contrastDescriptions
      return(res)
    }
    
    # check for the case where neither are present
    # as comparisons against base level
    if ( ! (contrastNumColumn %in% resNames &
              contrastDenomColumn %in% resNames) ) {
      stop(paste(contrastNumLevel,"and",contrastDenomLevel,"should be levels of",contrastFactor,
                 "such that",contrastNumColumn,"and",contrastDenomColumn,
                 "are contained in 'resultsNames(object)'"))
    }
    
    # case 2: expanded model matrices: build the appropriate contrasrt
    # these coefficient names have the form "factorlevel"
    # output: contrastNumColumn and contrastDenomColumn
  } else {
    
    # else in the expanded case, we first check validity
    contrastNumColumn <- make.names(paste0(contrastFactor, contrastNumLevel))
    contrastDenomColumn <- make.names(paste0(contrastFactor, contrastDenomLevel))
    if ( ! (contrastNumColumn %in% resNames & contrastDenomColumn %in% resNames) ) {
      stop(paste("both",contrastNumLevel,"and",contrastDenomLevel,"are expected to be in
resultsNames(object), prefixed by",contrastFactor))
    }
  }
  stop("You've uncovered a corner case! Please contact the maintainer 
         with your dataset and reproducible steps for help.")
}
