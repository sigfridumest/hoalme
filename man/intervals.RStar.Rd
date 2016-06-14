\name{intervals.RStar}
\title{Confidence Intervals using the Modified Ratio "R Star"}
\usage{
\method{intervals}{RStar}(object, level = 0.95, \dots)
}
\alias{intervals.RStar}
\arguments{
 \item{object}{an \code{RStar} object.}
 \item{level}{an optional numeric value for the interval confidence
   level. Defaults to 0.95.}
 \item{\dots}{some methods for the generic may require additional
   arguments. None are used in this method.}
}
\description{
  Confidence intervals on the parameters associated with the linear mixed effects model fitted by
  \code{lmeObject} (used in the computation of the \code{RStar}) object are obtained. The intervals
  are obtained by using the standard normal approximation to
  the modified ratio.
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
cons.Orth.fit<- conslmeObjects(Orth.fit, precision = 10, scale = hyp.scales) # 'conslmeObjects'
rs.Orth.fit <- RStar(Orth.fit, cons.Orth.fit) # 'RStar' object
intervals(rs.Orth.fit) # intervals (natural parameterization)
}
