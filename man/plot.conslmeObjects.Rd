\name{plot.conslmeObjects}
\title{Plot method for 'conslmeObjects' objects}
\usage{
 \method{plot}{conslmeObjects}(x, uncons = TRUE, \dots)
}
\alias{plot.conslmeObjects}
\arguments{
 \item{x}{ an object of class \code{conslmeObjects}.}
 \item{uncons}{ a logical value indicating whether the (log) likelihood should be plotted against
   the unconstrained parameterization of the linear mixed model, or against its natural
  parameterization. Defaults to the unconstrained parameterization.} 
\item{\dots}{ some methods for the generic might require additional
     arguments. None are used in these methods.} 
}
\description{
  Plots for the (log) profle likelihood of a linear mixed effects model.
}
\value{
  Plots of the (log) profile likelihood. For every parameter present in
  the \code{conslmeObjects} object, a plot is produced, with lines indicating the 90%, 95% and 99%
  critical levels based on the standard normal approximation.
}
\author{Sigfrido Iglesias-Gonzalez \email{sigfrido@alumni.utoronto.ca}}
\seealso{
 \code{\link{conslmeObjects}},
 \code{\link{hypotheses}},
 \code{\link{intervals.signedRoot}},
 \code{\link{lmeObject}},
 \code{\link{plot.signedRoot}}
 \code{\link{signedRoot}}
}
\examples{
Orth.fit <- lmeObject(distance ~ 1 + Sex + I(age - 11), random = pdSymm(~ 1 + I(age - 11)),
	data = Orthodont) # 'lmeObject'
hyp.scales <- list(fixed = list(c(3,3), c(3,3) , c(3,3)), reStruct = list(list(c(2.1,2.1), c(1.4,2), c(2,1.5)))) # scales
cons.Orth.fit<- conslmeObjects(Orth.fit, precision = 10, scale = hyp.scales) # 'conslmeObjects'
plot(cons.Orth.fit, TRUE) # plots (unconstrained parameterization)
plot(cons.Orth.fit, FALSE) # plots (natural parameterization)
}
