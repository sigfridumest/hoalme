## file hoalme/R/logLik.R, v 0.0-1 2009/11/18
##
##  Copyright (C) 2009 Sigfrido Iglesias-Gonzalez
##
##  This file is part of the "hoalme" package for R.  This program is 
##  free software; you can redistribute it and/or modify it under the 
##  terms of the GNU General Public License as published by the Free 
##  Software Foundation; either version 2 of the License, or (at your 
##  option) any later version.
##
##  This library is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program; if not, write to the Free Software
##  Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
##  MA 02111-1307 USA or look up the web page 
##  http://www.gnu.org/copyleft/gpl.html.
##
##  Please send any comments, suggestions or errors found to:
##  Sigfrido Iglesias-Gonzalez
##  Email: sigfrido@alumni.utoronto.ca.  


##logLik methods

logLik.lmeObject <-
  function(object, ...){
	respSt <- respStruct(object)
	modelSt <- object$modelStruct
	secondDer <- secondDeriv(object)
	secondDer.det <- abs(determinant(secondDer))
	nuissecondDer.dets <- rep(0, ncol(secondDer))
	parnams <- dimnames(secondDer)[[ 1 ]]
	for(i in 1:ncol(secondDer)){
		nuissecondDer.dets[ i ] <- abs(det(secondDer[ -i, -i ]))
	}
	names(nuissecondDer.dets)  <- parnams
	obsDer <- attr(secondDer, "W.inv.res")
	class(obsDer) <- "groupsList"
	fixedMixedDer <- attr(secondDer, "W.inv.X")
	varMixedDer <- attr(secondDer, "W.inv.W.dot.W.inv.res")
	Vvectors <- V(respSt)
	dirDer <- -Sum(Crossprod(Vvectors, obsDer))
	fixedMixedDer.dir <- Sum(Crossprod(Vvectors, fixedMixedDer))
	varMixedDer.dir <- Sum(Crossprod(Vvectors, varMixedDer))
	nams<- dimnames(varMixedDer.dir)[[ 3 ]]
	varMixedDer.dir <- matrix(unlist(varMixedDer.dir), ncol = dim(varMixedDer.dir)[3], dimnames = list(dimnames(fixedMixedDer.dir)[[ 1 ]], nams))
	dimnames(varMixedDer.dir) <- list(dimnames(varMixedDer.dir)[[ 1 ]], nams)
	MixedDer.dir<- cbind(fixedMixedDer.dir, varMixedDer.dir)
	fact <- -crossprod(secondDer, solve(MixedDer.dir))
 	SSDer <- list(dirDeriv = dirDer, MixedDeriv = MixedDer.dir, fact = fact)
	result <- list(respStruct = respSt, 
		modelStruct = modelSt,
		secondDeriv = secondDer,
 		secondDeriv.det = secondDer.det,
		nuissecondDeriv.dets = nuissecondDer.dets,
		varPar = Coef(object, F),
		meanPars = object$coefficients,
		fixedPars = object$coefficients$fixed,
		randomPars = object$coefficients$random,
		secondDer.uncons = object$Hessian,
		unconsvarPar = object$unconsPar,
		sigmasq = object$sigma^2,
		logLik = object$logLik,
		SSDeriv = SSDer)
	class(result) <- "lmeLogLik"
	result
}

logLik.conslmeObject <-
  function(object, lmeObj, ...){
	respSt <- respStruct(object, lmeObj)
	modelSt <- object$modelStruct
	secondDer <- secondDeriv(object, lmeObj)
	obsDer <- attr(secondDer, "W.inv.res")
	class(obsDer) <- "groupsList"
	fixedMixedDer <- attr(secondDer, "W.inv.X")
	varMixedDer <- attr(secondDer, "W.inv.W.dot.W.inv.res")
	Vvectors <- V(lmeObj)
	dirDer <- -Sum(Crossprod(Vvectors, obsDer))
	fixedMixedDer.dir <- Sum(Crossprod(Vvectors, fixedMixedDer))
	varMixedDer.dir <- Sum(Crossprod(Vvectors, varMixedDer))
	nams<- dimnames(varMixedDer.dir)[[ 3 ]]
	varMixedDer.dir <- matrix(unlist(varMixedDer.dir), ncol = dim(varMixedDer.dir)[3], dimnames = list(dimnames(fixedMixedDer.dir)[[ 1 ]], nams))
	dimnames(varMixedDer.dir) <- list(dimnames(varMixedDer.dir)[[ 1 ]], nams)
	MixedDer.dir<- cbind(fixedMixedDer.dir, varMixedDer.dir)
 	SSDer <- list(dirDeriv = dirDer, MixedDeriv = MixedDer.dir)
	if(!is.null(object$fixedPars))
		fixedPars <- object$fixedPars
		else
			fixedPars <- object$coefficients$fixed
	result <- list(respStruct = respSt, 
		modelStruct = modelSt,
		secondDeriv = secondDer,
		varPar = Coef(object, F),
		meanPars = object$coefficients,
		fixedPars = fixedPars,
		randomPars = object$coefficients$random,
		secondDer.uncons = object$Hessian,
		unconsvarPar = object$unconsPar,
		sigmasq = object$sigma^2,
		logLik = object$logLik,
		SSDeriv = SSDer)
	class(result) <- "conslmeLogLik"
 	attr(result, "hypothesis") <- attr(object, "hypothesis")
	invisible(result)
}

logLik.conslmeObjects <-
  function(object, lmeObj = NULL, ...){
 	objects <- object$objects
	hypotheses <- object$nat.hypotheses
	nams <- dimnames(hypotheses)
	npar <- nrow(hypotheses)
	if(is.null(lmeObj))
		lmeObj <- attr(object, "lmeObject")
	result <- lapply(objects, logLik, lmeObj)
	attr(result, "npar") <- npar
	MLE <- hypotheses[ , 1]
	names(MLE) <- nams[[ 1 ]]
	attr(result, "MLE") <- MLE
	hypotheses <- matrix(c(hypotheses[ , -1]), nrow = npar)
	dimnames(hypotheses) <- list(nams[[ 1 ]], nams[[ 2 ]][ -1 ])
	attr(result, "hypotheses") <- hypotheses
	attr(result, "lmeObject") <-  lmeObj
	class(result) <- "conslmeLogLiks"
	result
}

