\name{plot.QStatistic}
\title{Plot method for a 'QStatistic' object}
\usage{
 \method{plot}{QStatistic}(x, what = c( "QStatistic", "NP", "INF"), log = TRUE, \dots)
}
\alias{plot.QStatistic}
\arguments{
 \item{x}{ an object of class \code{QStatistic}.}
 \item{what}{ optional character string indicating what should be plotted. \code{"QStatistic"} plots
   sample space derivatives statistic 'Q', \code{"NP"} plots the nuisance parameterization
	adjustment, \code{"INF"} plots the departures from normality adjustment.
	Defaults to \code{"Qstatistic"}.}
 \item{log}{ optional logical value indicating whether a the logarithm should be used to plot the nuisance
   parameterization adjustment. Defaults to \code{TRUE}.}
\item{\dots}{ some methods for this generic require additional
     arguments. None are used in this method.} 
%  \item{evaluate}{If \code{TRUE} evaluate the new call else return the call.}
}
\description{
  Plots the modified Q statistic and/or other quantities of interest. 
}
\value{
  Plots of the Q statistic. For every parameter present in
  the \code{QStatistic} object, a plot is produced. Alternatively, other quantities are plotted
  depending on selection.
}
\author{Sigfrido Iglesias-Gonzalez \email{sigfrido@alumni.utoronto.ca}}
\seealso{
 \code{\link{QStatistic}},
 \code{\link{plot.RStar}},
 \code{\link{plot.signedRoot}},
 \code{\link{RStar}},
 \code{\link{signedRoot}}
}
\examples{
# 'signedRoot.conslmeObjects' method 
Orth.fit <- lmeObject(distance ~ 1 + Sex + I(age - 11), random = pdSymm(~ 1 + I(age - 11)),
									data = Orthodont) # 'lmeObject'
hyp.scales <- list(fixed = list(c(3,3), c(3,3) , c(3,3)), reStruct = list(list(c(2.1,2.1), c(1.4,2), c(2,1.5)))) # scales
cons.Orth.fit <- conslmeObjects(Orth.fit, precision = 10, scale = hyp.scales) # 'conslmeObjects'
qs.Orth.fit <- QStatistic(Orth.fit, cons.Orth.fit) # QStatistic object
plot(qs.Orth.fit, "Q") # plots the Q statistic
plot(qs.Orth.fit, "I") # plots the departure from normality adjustment
plot(qs.Orth.fit, "N") # plots nuisance parameterization adjustment
}
\keyword{models}
