\name{LA-methods}
\docType{methods}
\alias{LA-methods}
\alias{LA,eSet-method}
\alias{LA,matrix-method}
\alias{LA}
\title{Function to calculate LA estimate }
\description{
  'LA' is used to calculate the LA estimate for a gene triplet data. 
}

\arguments{
	\item{object}{An numerical matrix object with three columns or an object of ExpresionSet class with three features.}
	\item{dim}{An index of the column for the gene to be treated as the third controller variable. Default is dim=3}
	\item{geneMap}{A character vector with three elements representing the mapping between gene names and feature names (optional).}

}

\details{The input object can be a numerical matrix with three columns with row representing observations and column representing three variables. It can also be an ExpressionSet object with three features. If input a matrix class data, all three columns of the object representing the variables should have column names. Each variable in the object will be standardized with mean 0 and variance 1 in the function. In addition, the third variable will be quantile normalized within the function. More detail example about the usage of geneMap is demonstrated in the vignette.
}



\value{
  'LA' returns a numerical value representing the estimated value. A more detailed explanation of the value is illustrated in the vignette. 

}
\keyword{methods}


\references{Ker-Chau Li, Genome-wide coexpression dynamics: theory and application (2002). PNAS 99 (26): 16875-16880. 

Yen-Yi Ho, Leslie Cope, Thomas A. Louis, and Giovanni Parmigiani, GENERALIZED LIQUID ASSOCIATION (April 2009). Johns Hopkins University, Dept. of Biostatistics Working Papers. Working Paper
183. http://www.bepress.com/jhubiostat/paper183
}


\author{Yen-Yi Ho}
\seealso{\code{\link{GLA-methods}}, \code{\link{getsLA-methods}}}

\examples{
data<-matrix(rnorm(300), ncol=3)

colnames(data)<-c("Gene1", "Gene2", "Gene3")

LAest<-LA(data)

LAest

}









