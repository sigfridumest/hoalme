\name{signedRoot}
\title{Signed Root of the (log) Likelihood Ratio Test}
\usage{
 signedRoot(object, \dots)
 \method{signedRoot}{conslmeObjects}(object, \dots)
 \method{signedRoot}{RStar}(object, \dots)
}
\alias{signedRoot}
\alias{signedRoot.conslmeObjects}
\alias{signedRoot.RStar}
\alias{signedRoot.lmeObject}
\arguments{
 \item{object}{ an object from which a signed root of the (log) likelihood ratio is to be computed).}
 \item{\dots}{ some methods for this generic might require additional
     arguments. None are used in these methods.} 
}
\description{
  This generic function computes the signed root of the log likelihood ratio test.
  The methods here provided compute the test for every parameter present in the
  context of the linear mixed model on which a constrained optimization of its likelihood 
  is performed via \code{conslmeObjects} (usually over a normal region
  around the maximum likelihood estimator).
}
\value{
  An object of class \code{signedRoot}, which is a matrix. Use \code{hypotheses},
  \code{plot}, \code{intervals} methods to extract further information.
}
\author{Sigfrido Iglesias-Gonzalez \email{sigfrido@alumni.utoronto.ca}}
\seealso{
  \code{\link{conslmeObjects}},
  \code{\link{hypotheses}},
  \code{\link{intervals.signedRoot}},
  \code{\link{plot.signedRoot}},
  \code{\link{RStar}}
}
\examples{
# 'signedRoot.conslmeObjects' method 
Orth.fit <- lmeObject(distance ~ 1 + Sex + I(age - 11), random = pdSymm(~ 1 + I(age - 11)),
		data = Orthodont) # 'lmeObject'
hyp.scales <- list(fixed = list(c(3,3), c(3,3) , c(3,3)), reStruct = list(list(c(2.1,2.1), c(1.4,2), c(2,1.5)))) # scales
cons.Orth.fit <- conslmeObjects(Orth.fit, precision = 10, scale = hyp.scales) # 'conslmeObjects'
sr.Orth.fit <- signedRoot(cons.Orth.fit) # computing the likelihood ratio
plot(sr.Orth.fit, FALSE) # plots (unconstrained parameters) and intervals (natural parameterization)
intervals(sr.Orth.fit)

#extracting from an 'RStar' object
rs.Orth.fit <- RStar(Orth.fit, cons.Orth.fit) # 'RStar' object
signedRoot(rs.Orth.fit) # 'signedRoot'

}
