\name{plot.signedRoot}
\title{Plot method for a 'signedRoot' object}
\usage{
 \method{plot}{signedRoot}(x, natural = TRUE, \dots)
}
\alias{plot.signedRoot}
\arguments{
 \item{x}{ an object of class \code{signedRoot}.}
 \item{natural}{ a logical value indicating whether the test values should be plotted against
   the natural parameterization of the linear mixed model, or the unconstrained paratemterization
	should be used. Defaults to the natural paremeterization.} 
\item{\dots}{ some methods for this generic might require additional
     arguments. None are used in these methods.} 
%  \item{evaluate}{If \code{TRUE} evaluate the new call else return the call.}
}
\description{
  Plots for the signed root of the (log) likelihood ratio test in the context of a
  linear mixed effects model.
}
\value{
  Plots of the signed root of the (log) likelihood ratio test. For every parameter present in
  the \code{signedRoot} object, a plot is produced, with lines indicating the 95% and 99%
  critical levels based on the standard normal approximation.
}
\author{Sigfrido Iglesias-Gonzalez \email{sigfrido@alumni.utoronto.ca}}
\seealso{
 \code{\link{conslmeObjects}},
 \code{\link{hypotheses}},
 \code{\link{intervals.signedRoot}},
 \code{\link{signedRoot}}
}
\examples{
Orth.fit <- lmeObject(distance ~ 1 + Sex + I(age - 11), random = pdSymm(~ 1 + I(age - 11)),
	data = Orthodont) # 'lmeObject'
hyp.scales <- list(fixed = list(c(3,3), c(3,3) , c(3,3)), reStruct = list(list(c(2.1,2.1), c(1.4,2), c(2,1.5)))) # scales
cons.Orth.fit<- conslmeObjects(Orth.fit, precision = 10, scale = hyp.scales) # 'conslmeObjects'
sr.Orth.fit <- signedRoot(cons.Orth.fit) # computing the likelihood ratio
plot(sr.Orth.fit, FALSE) # plots (unconstrained parameterization)
plot(sr.Orth.fit, TRUE) # plots (natural parameterization)

}
\keyword{models}
