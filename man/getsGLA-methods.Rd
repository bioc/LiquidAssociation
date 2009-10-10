\name{getsGLA-methods}
\docType{methods}
\alias{getsGLA-methods}
\alias{getsGLA,eSet-method}
\alias{getsGLA,matrix-method}
\alias{getsGLA}
\title{ Function to calculate the sGLA test statistic for a given triplet data}
\description{
  'getsGLA' is used to calculate the sGLA test statistic and correponding p value.
}


\arguments{
	\item{object}{An numerical matrix object with three columns or an object of ExpresionSet class with three features.}
	\item{boots}{The number of bootstrap iterations for estimating the bootstrap standard error of sGLA. Default value is boots=30.}
	\item{perm}{The number of permutation iterations for generating the null distribution of the sGLA test statistic. Default is perm=100.}
	\item{cut}{cut==M +1. M is the number of grip points pre-specifed over the third variable.}
	\item{dim}{An index of the column for the gene to be treated as the third controller variable. Default is dim=3}
	\item{geneMap}{A character vector with three elements representing the mapping between gene names and feature names (optional).}
}


\details{The input object can be a numerical matrix with three columns with row representing observations and column representing three variables. It can also be an ExpressionSet object with three features. If input a matrix class data, all three columns of the object representing the variables should have column names. Each variable in the object will be standardized with mean 0 and variance 1 in the function. In addition, the third variable will be quantile normalized within the function. More detail example about the usage of geneMap is demonstrated in the vignette.
}

\value{
  'getsGLA' returns a vector with two elements. The first element is the value of test statistic and second element is the corresponding p value.  A more detailed interpretation of these values is illustrated in the vignette. 
 }


\keyword{methods}
\keyword{htest}

\references{Yen-Yi Ho, Leslie Cope, Thomas A. Louis, and Giovanni Parmigiani, GENERALIZED LIQUID ASSOCIATION (April 2009). Johns Hopkins University, Dept. of Biostatistics Working Papers. Working Paper
183. http://www.bepress.com/jhubiostat/paper183. 
}



\seealso{\code{\link{GLA-methods}}, \code{\link{getsLA-methods}}}
\examples{
data<-matrix(rnorm(300), ncol=3)

colnames(data)<-c("Gene1", "Gene2", "Gene3")

sGLAest<-getsGLA(data, boots=20, perm=100, cut=4, dim=3)

sGLAest

}

