


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



	


#' Compute the p-values for the enrichment of differentially expressed genes in the given gene sets/
#' @param object of class GSTTDataSet
#' @param gene_sets a list of gene sets
#' @return object of class GSTTResults
#' @export
get.gstt.results <- function(object, gene_sets) {
	if (is.null(object@permuted.results)) {
		object <- get.permuted.results(object)
	}

	bestGene_pvals = c()
	bestPathway_pvals = c()
	
	perm_gene_pvals = list()
	perm_path_pvals = list()
	i=0
	for (resPerm in object@permuted.results) {
		i <- i+1
		if (!("zscore" %in% colnames(resPerm))) {
			resPerm[,"zscore"] <- qnorm(1-resPerm[,"pvalue"]/2)
			resPerm[which(resPerm[,"zscore"]>10),"zscore"] = 10
		}
		all_genes <- row.names(resPerm)
		perm_gene_pvals[[i]] <- resPerm[,"pvalue"]
		bestGene_pvals <- c(bestGene_pvals, min(perm_gene_pvals[[i]], na.rm=T))
		perm_path_pvals[[i]] <- sapply(gene_sets, function(x) if (sum(!is.na(resPerm[x,"pvalue"]))>0) t.test(resPerm[x,"zscore"], resPerm[setdiff(all_genes,x),"zscore"],var.equal=T)$p.value else NA)
		bestPathway_pvals <- c(bestPathway_pvals, min(perm_path_pvals[[i]], na.rm=T))
		print("got pathway pvals")
	}

	
	Pathway_pvals <- sapply(gene_sets, function(x) if (sum(!is.na(object@dds.results[x,"pvalue"]))>0) t.test(object@dds.results[x,"zscore"], object@dds.results[setdiff(all_genes,x),"zscore"],var.equal=T)$p.value else NA)
	gene_names <- row.names(object@dds.results)

	path_fdr_1_pval <- median(bestPathway_pvals)
	path_fdr_.05_pval <- quantile(bestPathway_pvals, 0.05)
	path_fdr_1 <- sapply(which(Pathway_pvals < path_fdr_1_pval), function(ii) names(gene_sets)[[ii]])
	path_fdr_.05 <- sapply(which(Pathway_pvals < path_fdr_.05_pval), function(ii) names(gene_sets[[ii]]))

	gene_fdr_.05_pval <- quantile(bestGene_pvals, .05)
	gene_fdr_.05 <- sapply(which(object@dds.results[,"pvalue"]<gene_fdr_.05_pval), function(ii) gene_names[ii])
	gene_fdr_1_pval <- median(bestGene_pvals)
	gene_fdr_1 <- sapply(which(object@dds.results[,"pvalue"]<gene_fdr_1_pval), function(ii) gene_names[ii])

	#gene_sets_.05 <- sapply(gene_sets, function(x) list(x,length(intersect(x, gene_fdr_.05))))
	sets_table_.05 <- outer(path_fdr_.05, path_fdr_.05, Vectorize(function(x,y) {list(intersect(intersect(gene_sets[[x]],gene_sets[[y]]),gene_fdr_.05))}))
	##gene_sets_1 <- lapply(gene_sets, function(x) list(x,length(intersect(x, gene_fdr_1))))
	sets_table_1 <- outer(path_fdr_1, path_fdr_1, Vectorize(function(x,y) {list(intersect(intersect(gene_sets[[x]],gene_sets[[y]][[2]]),gene_fdr_1))}))

	gsr <- new("GSTTResults")
	gsr@best.pvals.set  <- bestPathway_pvals
	gsr@best.pvals.gene <- bestGene_pvals
	gsr@gene.sets <- gene_sets
	gsr@sets.pvals <- Pathway_pvals
	gsr@genes.pvals <- object@dds.results[,"pvalue"]
	names(gsr@genes.pvals) <- rownames(object@dds.results)
	
	gsr@sets_table_1 <- sets_table_1
	gsr@sets_table_.05 <- sets_table_.05
	return(gsr)
}

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
