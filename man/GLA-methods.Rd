\name{GLA-methods}
\docType{methods}
\alias{GLA-methods}
\alias{GLA,eSet-method}
\alias{GLA,matrix-method}
\alias{GLA}
\title{Function to calculate GLA estimate }
\description{
  'GLA' is used to calculate the GLA estimate for a gene triplet data. 
}

\arguments{
	\item{object}{An numerical matrix object with three columns or an object of ExpresionSet class with three features.}
	\item{cut}{cut==M +1. M is the number of grip points pre-specifed over the third variable.}
	\item{dim}{An index of the column for the gene to be treated as the third controller variable. Default is dim=3}
	\item{geneMap}{A character vector with three elements representing the mapping between gene names and feature names (optional).}

}

\details{The input object can be a numerical matrix with three columns with row representing observations and column representing three variables. It can also be an ExpressionSet object with three features. If input a matrix class data, all three columns of the object representing the variables should have column names. Each variable in the object will be standardized with mean 0 and variance 1 in the function. In addition, the third variable will be quantile normalized within the function. More detail example about the usage of geneMap is demonstrated in the vignette.
}


\value{
  'GLA' returns a numerical value representing the estimated value. A more detailed interpretation of the value is illustrated in the vignette. 

}
\keyword{methods}

\references{Yen-Yi Ho, Leslie Cope, Thomas A. Louis, and Giovanni Parmigiani, GENERALIZED LIQUID ASSOCIATION (April 2009). Johns Hopkins University, Dept. of Biostatistics Working Papers. Working Paper
183. http://www.bepress.com/jhubiostat/paper183
}

\author{Yen-Yi Ho}

\seealso{\code{\link{LA-methods}}, \code{\link{getsGLA-methods}}}

\examples{

data<-matrix(rnorm(300), ncol=3)

colnames(data)<-c("Gene1", "Gene2", "Gene3")
GLAest<-GLA(data, cut=4, dim=3)

GLAest

}


