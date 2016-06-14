\name{RStar}
\title{Modified Likelihood Ratio Test for Scalar Hypotheses}
\usage{
 RStar(object, \dots)
 \method{RStar}{lmeObject}(object, conslmeObjs, delta = 0.7, \dots)
 \method{RStar}{hoalme}(object, \dots)
}
\alias{RStar}
\alias{RStar.lmeObject}
\alias{RStar.hoalme}
\arguments{
 \item{object}{ an object from which a modified ratio "R Star" is to be computed.}
 \item{conslmeObjs}{ an object of class \code{conslmeObjects}.}
 \item{delta}{ a positive numeric value, preferabily between .5 and 1, indicating the region around zero to be left out
   in the computations; see details. Defaults to 0.7 .}
 \item{\dots}{ some methods for this generic may require additional
     arguments. None are used in these methods.}
}
\description{
  This generic function computes or extracts the modified likelihood ratio test using the techniques
  by Fraser and Reid. This modified ratio for scalar hypotheses has an asymptotic normal distribution to third order,
  which represents an improvement to the normal approximation to the signed root of the likelihood ratio. The
  improvement is based on the computation of the Fraser and Reid's Q statistic. The methods here provided
  computes the modified ratio in the context of the linear mixed model.
}
\value{
  An object of class \code{RStar}, which is a matrix, containing the modified
  ratio tests for every parameter present in \code{lmeObject}, within the
  ranges specified in \code{conslmeObjects}. Use \code{hypotheses},
  \code{plot}, \code{intervals} \code{QStatstic} methods to extract further information.
}
\details{
  A modified ratio test for a scalar hypothesis is a signed root of the likelihood ratio plus
  an adjustment for derpartures of normality and for the presence of nuisance parameters.
  In hypotheses close to the maximum likelihood estimator, the adjustment has numerical
  singularities. The points where these singularities are likely to be present are removed and
  substituted by a smoothed version based on the standard likelihood ratio and using the
  remaining points on the tails. As these issues arise, they might
  contaminate valid values on the tails. In general, the region where these numerical issues occur
  is hard to predict and visual checking is desirable. Use the \code{plot} method to asses the specification
  of the region via the \code{delta} argument.
}
\author{Sigfrido Iglesias-Gonzalez \email{sigfrido@alumni.utoronto.ca}}
\seealso{
  \code{\link{conslmeObjects}},
  \code{\link{hypotheses}},
  \code{\link{intervals.RStar}},
  \code{\link{plot.RStar}},
  \code{\link{QStatistic}}
  \code{\link{update.RStar}}  
}
\examples{
# usual procedure
Orth.fit <- lmeObject(distance ~ 1 + Sex + I(age - 11), random = pdSymm(~ 1 + I(age - 11)),
		data = Orthodont) # 'lmeObject'
hyp.scales <- list(fixed = list(c(3,3), c(3,3) , c(3,3)), reStruct = list(list(c(2.1,2.1), c(1.4,2), c(2,1.5)))) # scales
cons.Orth.fit<- conslmeObjects(Orth.fit, precision = 10, scale = hyp.scales) # 'conslmeObjects'
rs.Orth.fit <- RStar(Orth.fit, cons.Orth.fit) # 'RStar' object
plot(rs.Orth.fit) # plot
intervals(rs.Orth.fit) #intervals
hypotheses(rs.Orth.fit) # set of hypotheses used to compute the test curve

# using a 'hoalme' object'
hl.Orth.fit <- hoalme(Orth.fit, cons.Orth.fit) # 'hoalme' object
RStar(hl.Orth.fit) # 'RStar' object
}
