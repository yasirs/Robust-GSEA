#' get a GSTTDataSet from a matrix of transcript abundance counts
#' 
#' @param countData a matrix of transcript abundance counts
#' @param colData a data frame where rows represents samples (columns of contData) and columns represent properties of samples
#' @return a GSTTDataSet object
GSTTDataSetFromMatrix <- function(countData, colData, design, ...) {
  dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=design, ...)
  gds <- new("GSTTDataSet")
  gds@dds <- dds
  # cfac <- getContrast(dds, design=design)
  # gds@col.groups <- 
  # gds@dds.contrast <- c("time",60,0)
  return(gds)
}

#' Get a GSTTDataSet from a DESeq2 data saet
#' @param dds an object of class DESeqDataSet
#' @return an object of class GSTTDataSet
#' @export
#' @examples
#' example(DESeq)
#' gds <- GSTTDataSetFromDDS(dds)
GSTTDataSetFromDDS <- function(dds) {
  gds <- new("GSTTDataSet")
  gds@dds <- dds
  return(gds)
}

#' Compute the differential expression for a GSTTDataSet object
#' 
#' \code{diff.exp} runs DESeq for the dataset
#' 
#' @param object an object of class GSTTDataSet
#' @param contrast the contrast to compare the effects or levels on the expression, see \code{DESeq2::results}
#' @param name the name of the main effect, if not using contrast
#' @return an object of class GSTTDataSet with the results and col.groups slots filled
#' @export
#' @seealso \code{DESeq2::results}
diff.exp <- function(object, contrast, name, ...) {
  if (is.null(object@dds)) {
    print("Need data set! Exiting ...")
    return(object)
  }
  object@dds <- DESeq(object@dds)
  object@dds.results = results(object@dds, contrast, name, ...)
  object@col.groups <- getContrast(object@dds, contrast, name)
  object@dds.results[,"zscore"] <- qnorm(1-object@dds.results[,"pvalue"]/2)
  object@dds.results[which(object@dds.results[,"zscore"]>10),"zscore"] = 10
  return(object)
}

#' get balanced permutations
#' @return balanced permutations
#' 
#' @export
get.balanced.permutations <- function(cases, controls) {
  getgroups <- function(ii, case1s, control1s) {
    caseInd = ii[[1]]
    controlInd = ii[[2]] 
    noOverlap <- function(vec1, vec2) {!any(vec1 %in% vec2)}
    caseCond = sapply(seq_len(caseInd), function(x) noOverlap(case1s[x],case1s[caseInd]))
    controlCond = sapply(seq_len(length(control1s)), function(x) noOverlap(control1s[x],control1s[controlInd]))
    
    case2 = case1s[1:caseInd][caseCond]
    control2 = control1s[controlCond]
    if ((length(case2)<1)||(length(control2)<1))
      return(list())
    
    gr2inds = simplify2array(apply(as.matrix(expand.grid(seq_len(length(case2)),seq_len(length(control2)))),1,list))
    
    listOfBothGroups <- lapply(gr2inds, function(x) list(group1=unlist(c(case1s[caseInd],control1s[controlInd])),group2=unlist(c(case2[x[1]],control2[x[2]]))))
    return(listOfBothGroups)
  }
  nEach <- floor(min(length(cases),length(controls))/2)
  nCases <- length(cases)
  nControls <- length(controls)
  permutations <- list()
  case1s <- combn(cases, nEach) # possible samples from the cases for group 1
  case1s <- lapply(seq_len(ncol(case1s)), function(i) {case1s[,i]}) # transform to a list of vectors
  
  control1s <- combn(controls, nEach) # possible samples from the controls for group 1
  control1s <- lapply(seq_len(ncol(control1s)), function(i) {control1s[,i]}) # transform to a list of vectors
  
  allPairs <- simplify2array(apply(as.matrix(expand.grid(seq_len(length(case1s)),seq_len(length(control1s)))), 1, list))
  listOfGroups <- sapply(allPairs, function(x) getgroups(x,case1s,control1s))
  return(unlist(listOfGroups,recursive=F))
}


#' Run permutations for the GSTTDataSet
#' 
#' @param object of class GSTTDataSet
#' @return object of class GSTTDataSet with the slot permuted.results filled
#' @examples
#' example(DESeq)
#' gstt <- GSTTDataSetFromDDS(dds)
#' gstt <- diff.exp(gstt)
#' gstt <- run.permutations(gstt)
#' @export
run.permutations <- function(object) {
  if (is.null(object@dds.results)) {
    message("Differential expression has not been computed for the object.\n
            Attempting to run the diff.exp for the default contrast.")
    object <- diff.exp(object)
  }
  perm_results = list()
  col.data <- colData(object@dds)
  i=0
  for (perm in get.balanced.permutations(object@col.groups[[1]], object@col.groups[[2]])) {
    i <- i+1
    groups = rep(0,nrow(col.data))
    for (g1 in perm$group1) groups[g1] = 1
    for (g2 in perm$group2) groups[g2] = 2
    col.data$group <- factor(groups, levels=c(1,2,0))
    pdds <- DESeqDataSetFromMatrix(countData=counts(object@dds), colData=col.data, design = ~ group)
    pdds <- DESeq(pdds)
    resPerm <- results(pdds,contrast=c("group" ,1, 2))
    resPerm[,"zscore"] <- qnorm(1-resPerm[,"pvalue"]/2)
    resPerm[which(resPerm[,"zscore"]>10),"zscore"] = 10
    perm_results[i] <- resPerm
  }
  object@permuted.results <- perm_results
  return(object)
}