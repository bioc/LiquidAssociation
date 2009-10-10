\name{LiquidAssociation-package}
\alias{LiquidAssociation-package}
\alias{LiquidAssociation}
\docType{package}
\title{
    Liquid Associaiton Package
}
\description{
   The LiquidAssociation package provide methods to examine a special kind of  three-way interaction called liquid association. The term liquid association was first proposed by 
contains functions for calculate direct and model-based estimators for liquid associaiton. 
   It also provides functions for testing the existence of liquid associaiton given a gene triplet data.  
}
\details{
\tabular{ll}{
Package: \tab LiquidAssociation\cr
Type: \tab Package\cr
Version: \tab 1.0.4\cr
Date: \tab 2009-10-05\cr
License: \tab GPL version 2 or newer\cr
LazyLoad: \tab yes\cr
}
GLA
LA
CNM.full
CNM.simple
getsGLA
getsLA
plotGLA
}
\author{
Yen-Yi Ho <yho@jhsph.edu>

Maintainer: Yen-Yi Ho <yho@jhsph.edu>
}
\references{Ker-Chau Li, Genome-wide coexpression dynamics: theory and application (2002). PNAS 99 (26): 16875-16880. 


Yen-Yi Ho, Leslie Cope, Thomas A. Louis, and Giovanni Parmigiani, GENERALIZED LIQUID ASSOCIATION
(April 2009). Johns Hopkins University, Dept. of Biostatistics Working Papers. Working Paper
183. http://www.bepress.com/jhubiostat/paper183

}
\keyword{package }



\examples{
data<-matrix(rnorm(300), ncol=3)

colnames(data)<-c("Gene1", "Gene2", "Gene3")
FitCNM.full<-CNM.full(data)

FitCNM.full

}
