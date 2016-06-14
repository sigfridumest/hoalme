\name{QStatistic}
\title{Sample Space Derivatives Statistic 'Q'}
\usage{
 QStatistic(object, \dots)
 \method{QStatistic}{lmeObject}(object, conslmeObjs, \dots)
 \method{QStatistic}{RStar}(object, \dots)
}
\alias{QStatistic}
\alias{QStatistic.lmeObject}
\alias{QStatistic.RStar}
\arguments{
 \item{object}{ an object from which a Q statistic can be computed or extracted.}
 \item{conslmeObjs}{ an object of class \code{conslmeObjects}.}
 \item{\dots}{ some methods for this generic require additional
     arguments. None are used in these methods.} 
}
\description{
  This function computes or extract sample space derivatives statistic (Q statistic) using the techniques
  by Fraser and Reid. This statistic is a necessary ingredient for computing the modified ratio
  test for scalar hypotheses.
}
\value{
  An object of class \code{QStatistic}, which is a matrix, containing Q statistic
  computed for every parameter present in \code{object}, within the
  ranges specified in \code{conslmeObjects}. Use \code{plot} method for visualizing information
  regarding the nuisance parameter and departures from normality adjustments.
}
\details{
  A necessary ingredient for computing a modified ratio test is the calculation of
  a sample space derivatives statistic, in this case called "Q statistic". The statistic
  is composed by a combination of an adjustment for departures of normality
  of the standard likelihood ratio (INF), and an adjustement accounting for the presence
  of nuisance parameters (NP). The correction term resulting from the Q statistic provides
  a third order normal approximation to the modified ratio test.
}
\author{Sigfrido Iglesias-Gonzalez \email{sigfrido@alumni.utoronto.ca}}
\seealso{
  \code{\link{conslmeObjects}},
  \code{\link{hypotheses}},
  \code{\link{intervals.RStar}},
  \code{\link{plot.RStar}},
  \code{\link{QStatistic}}
}
\examples{
# direct computation
Orth.fit <- lmeObject(distance ~ 1 + Sex + I(age - 11), random = pdSymm(~ 1 + I(age - 11)),
									data = Orthodont) # 'lmeObject'
hyp.scales <- list(fixed = list(c(3,3), c(3,3) , c(3,3)), reStruct = list(list(c(2.1,2.1), c(1.4,2), c(2,1.5)))) # scales
cons.Orth.fit<- conslmeObjects(Orth.fit, precision = 10, scale = hyp.scales) # 'conslmeObjects'
qs.Orth.fit <- QStatistic(Orth.fit, cons.Orth.fit) # 'QStatistic' object
plot(qs.Orth.fit) # plots

# via 'RStar'
rs.Orth.fit <- RStar(Orth.fit, cons.Orth.fit) # 'RStar' object
QStatistic(rs.Orth.fit)
}
