\name{CNM.simple-methods}
\docType{methods}
\alias{CNM.simple-methods}
\alias{CNM.simple,eSet-method}
\alias{CNM.simple,matrix-method}
\alias{CNM.simple}
\title{The function fits a reduced conditional normal model (CNM)}
\description{
 'CNM.simple' is used to fit the reduced (correlation only) conditional normal model using GEE.}


\arguments{
	\item{object}{An numerical matrix object with three columns or an object of ExpressionSet class with three features.}
	\item{geneMap}{A character vector with three elements representing the mapping between gene names and feature names (optional).}
	\item{dim}{An index of the column for the gene to be treated as the third controller variable. The default value is dim=3. }

}

\details{The input object can be a numerical matrix with three columns with row representing observations and column representing three variables. It can also be an ExpressionSet object with three features. If input a matrix class data, all three columns of the object representing the variables should have column names. Each variable in the object will be standardized with mean 0 and variance 1 in the function. In addition, the third variable will be quantile normalized within the function. More detail example about the usage of geneMap is demonstrated in the vignette.
}


\value{
  'CNM.full' returns a object of CNM class with two Slots. The first slot describes the fitted model. The second slot is a matrix contains the CNM model fitting results. The row of this matrix represents the parameters in the CNM model. The first column, estimates, is the estimated value of the corresponding parameters. The second column, san.se, is the value of sandwich standard error estimator for the estimates. The third column, wald, is the wald test statistic as described in Ho et al (2009). The corresponding p value for the wald test statistic is represented in the fourth column. A more detailed interpretation of these values is illustrated in the vignette. 
}


\keyword{methods}
\keyword{models}

\references{Yen-Yi Ho, Leslie Cope, Thomas A. Louis, and Giovanni Parmigiani, GENERALIZED LIQUID ASSOCIATION (April 2009). Johns Hopkins University, Dept. of Biostatistics Working Papers. Working Paper
183. http://www.bepress.com/jhubiostat/paper183. 

Jun Yan and Jason Fine. Estimating equations for association structuresStatistics in Medicine. 23(6): 859--74; discussion 875-7,879-80. http://dx.doi.org/10.1002/sim.1650
}

\author{Yen-Yi Ho}

\seealso{\code{\link{CNM.full-methods}}, \code{\link{CNM-class}}}

\examples{

data<-matrix(rnorm(300), ncol=3)

colnames(data)<-c("Gene1", "Gene2", "Gene3")

FitCNM.simple<-CNM.simple(data)

FitCNM.simple


}

