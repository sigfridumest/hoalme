\name{intervals.signedRoot}
\title{Confidence Intervals using the Signed Root of the (log) Likelihood Ratio Test}
\usage{
\method{intervals}{signedRoot}(object, level = 0.95, \dots)
}
\alias{intervals.signedRoot}
\arguments{
 \item{object}{an \code{signedRoot} object.}
 \item{level}{an optional numeric value for the interval confidence
   level. Defaults to 0.95.}
 \item{\dots}{some methods for the generic may require additional
   arguments.}
}
\description{
  Confidence intervals on the parameters associated with the linear mixed effects model fitted by
  based on the signed root of the (log) likelihood ratio test are obtained.
}
\value{
  A matrix containing estimators and confidence limits.
}

\author{Sigfrido Iglesias-Gonzalez, \email{sigfrido@alumni.utoronto.ca}}
\seealso{
 \code{\link{intervals}},
 \code{\link{intervals.signedRoot}}
}

\examples{
Orth.fit <- lmeObject(distance ~ 1 + Sex + I(age - 11), random = pdSymm(~ 1 + I(age - 11)),
		data = Orthodont) # 'lmeObject'
hyp.scales <- list(fixed = list(c(3,3), c(3,3) , c(3,3)), reStruct = list(list(c(2.1,2.1), c(1.4,2), c(2,1.5)))) # scales
cons.Orth.fit <- conslmeObjects(Orth.fit, precision = 10, scale = hyp.scales) # 'conslmeObjects'
sr.Orth.fit <- signedRoot(cons.Orth.fit) # computing the likelihood ratio
intervals(sr.Orth.fit) # intervals 
}