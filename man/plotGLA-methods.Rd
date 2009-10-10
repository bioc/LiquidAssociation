\name{plotGLA-methods}
\docType{methods}
\alias{plotGLA-methods}
\alias{plotGLA,eSet-method}
\alias{plotGLA,matrix-method}
\alias{plotGLA}
\title{ The function plots scatter plots of two variables conditioning on the value of a third variable. }
\description{
	'plotGLA' is a function to plot the scatter plots of two variables conditioning on the value of a third variable. 
}


\arguments{
	\item{object}{An numerical matrix object with three columns or an object of ExpresionSet class with three features..}
	\item{cut}{cut==M +1. M is the number of grip points pre-specifed over the third variable .}
	\item{dim}{An index of the column for the gene to be treated as the third controller variable.}
	\item{filen}{The file name for the output graph can be specified when save=TRUE}
	\item{save}{If save=TRUE then output graphs will be save as PDF files with file name as specified by filen.}
	 \item{geneMap}{A character vector with three elements representing the mapping between gene names and feature names (optional).}
	\item{...}{Other graphical parameters can be passed to function plot.}
}

\details{The input object can be a numerical matrix with three columns with row representing observations and column representing three variables. It can also be an ExpressionSet object with three features. More detail example about the usage of geneMap is demonstrated in the vignette.
}




\keyword{methods}

\references{Yen-Yi Ho, Leslie Cope, Thomas A. Louis, and Giovanni Parmigiani, GENERALIZED LIQUID ASSOCIATION (April 2009). Johns Hopkins University, Dept. of Biostatistics Working Papers. Working Paper
183. http://www.bepress.com/jhubiostat/paper183
}

\author{Yen-Yi Ho}

\examples{
data<-matrix(rnorm(300), ncol=3)

colnames(data)<-c("Gene1", "Gene2", "Gene3")

plotGLA(data, cut=3, dim=3, pch=16, filen="GLAplot", save=FALSE)

}

