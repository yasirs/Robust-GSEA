setClass("GSTTResults",
         representation(best.pvals.set="numeric",
                        best.pvals.gene="numeric",
                        gene.sets="list",
                        sets.pvals="numeric",
                        genes.pvals="numeric",
                        sets_table_1="matrix",
                        sets_table_.05="matrix"),
         prototype(best.pvals.set=NULL, best.pvals.gene=NULL, gene.sets=NULL, sets.pvals=NULL, genes.pvals=NULL))


setClass("GSTTDataSet",
         representation(dds="DESeqDataSet",
                        dds.results="DESeqResults",
                        col.groups="list",
                        permuted.results="list"),
         prototype(dds=NULL, dds.results=NULL, col.groups=NULL, permuted.results=NULL))