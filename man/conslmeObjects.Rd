
\name{conslmeObjects}
\title{Constrained Optimization in the context of the Linear Mixed-Effects Model}
\usage{
conslmeObjects(object, scale = 2, precision = 10, control = list(returnObject = TRUE))
}
\alias{conslmeObjects}
\alias{conslmeObjects.lmeObject}
\arguments{
 \item{object}{ an object of class \code{lmeObject}.}
%  \item{hypotheses}{ a list containing the set of parameter values in their unconstrained
%   form representing the hypotheses at which the profile likelihood is calculated. If present,
%   the \code{scales} argument is irrelevant. The list must contain the components \code{fixed}
%   and \code{reStruct}, which in turn are lists. If p is the number of coefficients describing
%   the fixed effects, then \code{fixed} is a list of p components containing a   
%   Defaults to no hypotheses.}
 \item{scale}{ an optional list with components \code{fixed} and \code{reStruct}
  containing pairs of positive numeric values representing the limits of a normalized
  interval around the maximum likelihood estimator for each parameter present in the
  model fitted by \code{lmeObject}. Alternatively, a single positive
  numeric value to be used in the construction of these limits. Common scale values
  are between 1 and 2. See Details below. Both the \code{fixed} and \code{reStruct}
  components are in turn lists. If p is the number of coefficients describing the fixed
  effects, then \code{fixed} is a list of p pairs of scales defining the inferior and
  superior limits. The \code{reStruct} list must be specified in the same way,
  except that the specification should be done in a list per each level present in the model;
  \code{reStruct} is a list of lists. If there are q levels, then \code{reStruct} should
  contain q lists, each of them with a number of pairs equal to the number of
  variance parameters at the level. Further, if a \code{pdBlocked} structure is
  present in th model fit, a similar approach should be followed for each block at each level;
  see examples below. Defaults to a general scale equal to 2.}
 \item{precision}{ an optional even integer value specifying the number of points at which the parameters
  values are to be fixed during the constrained optimization. Too few points lead to imprecise results when
  computing higher order methods. Defaults to 10.}
 \item{control}{see \code{lmeControl} for details. Defaults to \code{lmeControl} default, except that the routine
  is not terminated if convergence fails.}
}

\description{
  Constrained optimization of the likelihood function of a linear mixed effects model as specified
  in a \code{lmeObject} fit. The constrained optimization is performed fixing hypothetical values
  of the parameters present in the model. This gives the profile likelihood at the hypothetical values.
  This is a necessary step for the computation of higher order inference methods. Routines
  in \code{nlme} package for fitting linear mixed effects models use an "unconstrained"
  parameterization as opposed to the "natural" parameterization in which the model is formulated.
  The constrained optimization is performed at the unconstrained formulation.
}

\details{
 The constrained optimization is performed by fixing a parameter value and optimizing
 on the remaining parameters. This is done for all the paremeters present in the model, except
 for the residual variance; if the model includes only one coefficient describing the fixed effects,
 the constrained optimization is performed only on the variance parameters. A set of normalized
 deviations around the maximum likelihood estimator is constructed using the scales provided
 in the \code{scale} argument. For example, a symmetric interval of normalized length equal to
 4 is constructed by adding and substracting two times the standard error to the maximum likelihood
 estimator; a pair (2, 2) of scale values is thus used. Asymmetric intervals are allowed. For this
 implementation to succeed it is necessary to specify parameter values not too far from the overall
 maximum likelihood estimator. This is accomplished by building a normalized interval around
 the estimator, using its asymptotic standard error and solving. The procedure is equivalent to
 constructing a confidence interval for a parameter using the central limit theorem, where the scale
 plays the role of the standard normal quantile required to achieve the confidence level.
}

\value{
 A list with components components \code{logProfLiks} , (matrix with profile likelihoods
 for each parameter), \code{parameters} (array with estimated parameteres),
 \code{hypotheses} (matrix with hypotheses in unconstrained form),
 \code{nat.hypotheses} (matrix with hypotheses in natural form), \code{objects} (list with optimization
 results), \code{npar} (number of parameters on which hypotheses are formulated),
 \code{npar.fixed} (number of parameters describing fixed effects on which hypotheses are formulated).
}

\author{Sigfrido Iglesias-Gonzalez \email{sigfrido@alumni.utoronto.ca}}
\seealso{
 \code{\link{hypotheses}},
 \code{\link{hoalme}},
 \code{\link{lmeObject}},
 \code{\link{QStatistic}},
 \code{\link{plot.conslmeObjects}},
 \code{\link{RStar}},
 \code{\link{signedRoot}}
}

\examples{
# 'lmeObject' fit:
Orth.fit <- lmeObject(distance ~ 1 + Sex + I(age - 11), random = pdSymm(~ 1 + I(age - 11)), data = Orthodont)
# specyfing hypotheses
hyp.scales <- list(fixed = list(c(3, 3), c(3, 3) , c(3, 3)), reStruct = list( list(c(2.1, 2.1), c(1.4, 2), c(2, 1.5)) ) )
# alternatively a general scale value for the fixed effect part
hyp.scales <- list(fixed = 3, reStruct = list( list(c(2.1, 2.1), c(1.4, 2), c(2, 1.5)) ) )
#constrained optimization
cons.Orth.fit <- conslmeObjects(Orth.fit, scale = hyp.scales)
# plotting the profile likelihood
plot(cons.Orth.fit)

#two levels
WaferSubset.fit <- lmeObject(current ~ poly(voltage, 2), random = list(Wafer = pdDiag(~ 1 + poly(voltage, 1)),
Site = pdDiag(~ 1 + poly(voltage, 1))), data = WaferSubset)
# hypotheses: a list per level
hyp.scales <- list(fixed = 2.7, reStruct = list( list(c(2.2, 2.2), c(1.2, .5)), list(c(2.4, 2.3), c(2.4, 2.2) ) ) )
cons.WaferSubset.fit <- conslmeObjects(WaferSubset.fit, scale = hyp.scales, precision = 6)
plot(cons.WaferSubset.fit) # constrained parameterization
plot(cons.WaferSubset.fit, FALSE) # natural paramaterization

# 'pdBlocked' structure
# one level, two blocks, one parameter per block
# only one coefficient describing the fixed effects, no profile likelihood for this parameter
vc.fit <- lmeObject( y ~ 1, VarComp, random = pdBlocked( list(pdIdent(~ block1 - 1), pdIdent(~ block2 - 1)) ) ) 
hyp.scales <- list(reStruct = list( list(  list(c(3.5, 1.8)), list(c(3.2, 1.8)) )) ) # one list for one level, and one list per block
cons.vc.fit <- conslmeObjects(vc.fit, scale = hyp.scales)
plot(cons.vc.fit)

}
