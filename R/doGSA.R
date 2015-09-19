

#' @export
makeData <- function(fname='transcripts.counts.matrix') {
  m <- as.matrix(read.table(fname))
  all_genes <- row.names(m)
  colData = data.frame(row.names=colnames(m),replicate=c(11,11,11,12,12,12,13,13,13), time=rep(c(0,15,60),3))
  colData$time <- factor(colData$time, levels=c(0,15,60))
  counts <- round(m)
  dds <- DESeqDataSetFromMatrix(countData=counts, colData=colData, design = ~ time)
  gds <- new("GSTTDataSet")
  gds@dds <- dds
  return(gds)
}

#' get a vector of gene sets from a file
#' 
#' @param fname name of the input file
#' @return vector of gene sets
#' @export
#' @export
read.sets <- function(fname) {
  readOneSet <- function(line) {
    w = strsplit(line, "\t")
    list(w[[1]][1], w[[1]][-1])
  }
  f = file(fname, "r")
  gene_sets <- lapply(readLines(f, warn=FALSE), readOneSet)
  close(f)
  gnames <- sapply(gene_sets, function(x) x[[1]])
  gene_sets <- sapply(gene_sets, function(x) x[[2]])
  names(gene_sets) <- gnames
  return(gene_sets)
}



#' @export
sampleData <- function(num.genes=5000) {
	library("pasilla")
	library("Biobase")
	data("pasillaGenes")
	countData <- counts(pasillaGenes)
  if (num.genes>nrow(countData)) num.genes=nrow(countData)
	countData <- countData[1:num.genes,]
	colData <- pData(pasillaGenes)[,c("condition","type")]
	gds <- GSTTDataSetFromMatrix(countData= countData,colData= colData,design= ~ condition)
	return(gds)
}

#' @export
get.sets.table <- function(gsr, fwer)
{
	path_fwer_pval <- quantile(gsr@best.pvals.set, fwer)
	path_fwer <- sapply(which(gsr@sets.pvals < path_fwer_pval), function(ii) gsr@gene.sets[[ii]][[1]])

	gene_fwer_pval <- quantile(gsr@best.pvals.gene, fwer)
	gene_fwer <- names(which(gsr@genes.pvals<gene_fwer_pval))

	gene_sets <- lapply(gsr@gene.sets, function(x) list(x[[1]],x[[2]],length(intersect(x[[2]], gene_fwer))))
	sets_table <- outer(path_fwer, path_fwer, Vectorize(function(x,y) {length(intersect(gsr@gene.sets[[x]][[2]], gsr@gene.sets[[y]][[2]]))/length(union(gsr@gene.sets[[x]][[2]], gsr@gene.sets[[y]][[2]]))}))
	rownames(sets_table) <- path_fwer
	colnames(sets_table) <- path_fwer
	h <- hclust(as.dist(1-sets_table))
	sets_table <- sets_table[h$order, h$order]
	return(sets_table)
}
