% Generated by roxygen2 (4.0.1): do not edit by hand
\name{GSTTDataSetFromDDS}
\alias{GSTTDataSetFromDDS}
\title{Get a GSTTDataSet from a DESeq2 data saet}
\usage{
GSTTDataSetFromDDS(dds)
}
\arguments{
\item{dds}{an object of class DESeqDataSet}
}
\value{
an object of class GSTTDataSet
}
\description{
Get a GSTTDataSet from a DESeq2 data saet
}
\examples{
example(DESeq)
gds <- GSTTDataSetFromDDS(dds)
}

