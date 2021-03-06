% Generated by roxygen2 (4.0.1): do not edit by hand
\name{run.permutations}
\alias{run.permutations}
\title{Run permutations for the GSTTDataSet}
\usage{
run.permutations(object)
}
\arguments{
\item{object}{of class GSTTDataSet}
}
\value{
object of class GSTTDataSet with the slot permuted.results filled
}
\description{
Run permutations for the GSTTDataSet
}
\examples{
example(DESeq)
gstt <- GSTTDataSetFromDDS(dds)
gstt <- diff.exp(gstt)
gstt <- run.permutations(gstt)
}

