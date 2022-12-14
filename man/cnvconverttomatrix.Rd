% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Data_preprocessing.R
\name{cnvconverttomatrix}
\alias{cnvconverttomatrix}
\title{Convert to matrix}
\usage{
cnvconverttomatrix(CNV)
}
\arguments{
\item{CNV}{Copy number variation data}
}
\value{
A matrix
}
\description{
CNV data into an expression matrix
}
\examples{
data(CNV)
CNVmatrix <- cnvconverttomatrix(CNV)

}
