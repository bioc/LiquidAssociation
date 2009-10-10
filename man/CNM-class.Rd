\name{CNM-class}
\docType{class}
\alias{CNM-class}
\alias{print,CNM-method}
\alias{show,CNM-method}

\title{Class CNM}
\description{This is a class representation for CNM model fitting results. }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CNM", ...)} or the functions \code{\link{CNM.full-methods}} and \code{\link{CNM.simple-methods}}
}
\section{Slots}{
	 \describe{
    \item{\code{Model}:}{Object of class \code{character} representing the fitted CNM model.}
    \item{\code{output}:}{Object of class \code{matrix} representing the parameter estimates from the fitted CNM model.}
  }
}
\section{Methods}{
  \describe{
    \item{print}{\code{signature(x = "CNM")}: Display CNM model fitting result.}
    \item{show}{\code{signature(object = "CNM")}: Display CNM model fitting result.}
	 }
}
\references{Yen-Yi Ho, Leslie Cope, Thomas A. Louis, and Giovanni Parmigiani, GENERALIZED LIQUID ASSOCIATION (April 2009). Johns Hopkins University, Dept. of Biostatistics Working Papers. Working Paper
183. http://www.bepress.com/jhubiostat/paper183. Yen-Yi Ho, Leslie Cope, Thomas A. Louis, and Giovanni Parmigiani, GENERALIZED LIQUID ASSOCIATION (April 2009). Johns Hopkins University, Dept. of Biostatistics Working Papers. Working Paper
183. http://www.bepress.com/jhubiostat/paper183. }
\author{ Yen-Yi Ho}
\note{ The usage of this class is demonstrated in the vignette.}
\seealso{
	 related methods \code{\link{print}}, \code{\link{show}}.
}
\examples{
showClass("CNM")
}
\keyword{classes}
