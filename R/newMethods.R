
## file hoalme/R/newMethods.R, v 0.0-1 2009/11/18
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


##New methods for generics from other packages

model.matrix.grps.model.frame <-
 function(object, random = FALSE, ...){
	lapply(object, function(object, random){
	cols <- attributes(object)
	Q <- attr(object, "Q")
	if(random)
		result <- as.data.frame(object[, cols[["Zcols"]] ])
		else
			result <- as.data.frame(object[, cols[["Xcols"]] ])
	result
	}, random)
}

residuals.grps.model.frame <-
  function(object, ...){
	lapply(object, function(object){
		cols <- attributes(object)
		residuals <- object[, cols$rescol ]
		as.matrix(residuals)
		})
}

##chol methods

chol.pdDiag <-
  function(x, ...){
	x <- as.matrix(x)
	diag(x) <- sqrt(diag(x))
	x
}

chol.pdSymm <-
  function(x, ...){
	x <- as.matrix(x)
	x <- chol(x)
	x
}

chol.pdLogChol <-
  function(x, ...){
	x <- as.matrix(x)
	x <- chol(x)
	x
}

chol.pdCompSymm <-
  function(x, ...){
	x <- as.matrix(x)
	x <- chol(x)
	x
}

chol.pdIdent <-
  function(x, ...){
	x <- as.matrix(x)
	x <- sqrt(x)
	x
}

chol.pdBlocked <-
  function(x, ...){
	result <- lapply(x, chol)
	result
}

chol.reStruct <-
  function(x, ...){
	result <- lapply(x, chol)
	result
}

chol.lmeStructInt <-
  function(x, ...){
	result <- lapply(x, chol)
	result
}

##intervals methods

intervals.RStar <-
  function(object, level = 0.95, ...){
  	lmeObj <- attr(object, "lmeObject")
  	lmeSt <- lmeObj$modelStruct
 	non.neg.limits <- Diag(lmeSt)
 	npar.fixed <- attr(object, "npar.fixed")
 	if(npar.fixed > 0)
 	  	mles <- c(lmeObj$coefficients$fixed, Coef(lmeObj, F))
 	  	else
 	  		mles <- Coef(lmeObj, F)
  	p <- (1 - level)/2
  	quant <- qnorm(p)
  	hypotheses <- hypotheses(object)
  	npar <- nrow(object)
  	result <- matrix(0, nrow = npar, ncol = 3)
  	for(i in 1:npar){
		spline.obj <- smooth.spline(y = hypotheses[ i, ], x = object[ i, ], all.knots = TRUE, spar = 1e-07)
		lims <- predict(spline.obj, x = c(-quant, quant))$y
		if(i > npar.fixed)
			if( lims[1] < 0 && non.neg.limits[ i - npar.fixed] ) lims[1] <- 0
		result[ i , c(1, 3)] <- lims
		result[ i , 2 ] <- mles[ i ]
	}
	dimnames(result) <- list(dimnames(hypotheses)[[1]], c("lower", "estimate", "upper"))
	attr(result, "npar.fixed") <- npar.fixed
	attr(result, "level") <- level
	attr(result, "lmeObject") <- lmeObj
	class(result) <- "intervals.RStar"
	result
}

intervals.signedRoot <-
  function(object, level = 0.95, ...){
  	lmeObj <- attr(object, "lmeObject")
  	lmeSt <- lmeObj$modelStruct
 	non.neg.limits <- Diag(lmeSt)
 	npar.fixed <- attr(object, "npar.fixed")
 	if(npar.fixed > 0)
 	  	mles <- c(lmeObj$coefficients$fixed, Coef(lmeObj, F))
 	  	else
 	  		mles <- Coef(lmeObj, F)
  	p <- (1 - level)/2
  	quant <- qnorm(p)
  	hypotheses <- hypotheses(object)
  	npar <- nrow(object)
  	result <- matrix(0, nrow = npar, ncol = 3)
  	for(i in 1:npar){
		spline.obj <- smooth.spline(y = hypotheses[ i, ], x = object[ i, ], all.knots = TRUE)
		lims <- predict(spline.obj, x = c(-quant, quant))$y
		if(i > npar.fixed)
			if( lims[1] < 0 && non.neg.limits[ i - npar.fixed] ) lims[1] <- 0
		result[ i , c(1, 3)] <- lims
		result[ i , 2 ] <- mles[ i ]
	}
	dimnames(result) <- list(dimnames(hypotheses)[[1]], c("lower", "estimate", "upper"))
	attr(result, "npar.fixed") <- npar.fixed
	attr(result, "level") <- level
	attr(result, "lmeObject") <- lmeObj
	class(result) <- "intervals.signedRoot"
	result
}

intervals.lmeObject <-
  function(object, level = 0.95, ...){
	class(object) <- c("lme", "Object")
	intervals(object, level, ...)
}


