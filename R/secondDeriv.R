## file hoalme/R/secondDeriv.R, v 0.0-1 2009/11/18
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


##secondDeriv methods

secondDeriv <-
  function(object, ...)
	UseMethod("secondDeriv")

secondDeriv.lmeObject <-
  function(object, ...){
	respSt <- respStruct(object)
	W.inv <- Var(respSt, T)
	X <- model.matrix(respSt)
	res <- residuals(respSt)
# 	U.inv <- lapply(W.inv, function(x) t(attr(x, "cholvar.inv")))
	W.dot <- pdDeriv(respSt)
# 	class(U.inv) <- "groupsList"
	W.inv.X <- Crossprod(W.inv, X)
# 	secDer.fixed <- -Sum(Crossprod(X, W.inv.X))
	secDer.fixed <- -solve(object$varFix)
	W.dot.W.inv <- Crossprod(W.dot, W.inv) #will be used transposed
	W.inv.res <- Crossprod(W.inv, res)
	W.inv.W.dot.W.inv.res <- Crossprod(W.dot.W.inv, W.inv.res)
	secDer.fixed.var <- -Sum(Crossprod(X, W.inv.W.dot.W.inv.res))
	secDer.fixed.var <- matrix(c(secDer.fixed.var), ncol = dim(secDer.fixed.var)[ 3 ])
	secDer.var <- lmeObjectHessian(object)
	nams <- dimnames(W.dot[[ 1 ]])[[3]]
	dimnames(secDer.fixed.var) <- list(dimnames(X[[ 1 ]])[[2]], nams)
	dimnames(secDer.var) <- list(nams, nams)
	A <- cbind(secDer.fixed, secDer.fixed.var)
	B <- cbind(t(secDer.fixed.var), secDer.var)
	secDer <- rbind(A, B)
	attr(secDer, "secDer.fixed") <- secDer.fixed
	attr(secDer, "secDer.fixed.var") <- secDer.fixed.var
	attr(secDer, "secDer.var") <- secDer.var
	attr(secDer, "W.inv.res") <- W.inv.res
	attr(secDer, "W.inv.W.dot.W.inv.res") <- W.inv.W.dot.W.inv.res
	attr(secDer, "W.inv.X") <- W.inv.X
	attr(secDer, "variance estimates") <- Coef(object, F)
	class(secDer) <- c("secondDeriv")
	secDer
}

print.secondDeriv <-
  function(x, ...){
	print(x[ , ], ...)
}

secondDeriv.conslmeObjects <-
   function(object, ...){
	lmeObj <- attr(object, "lmeObject")
	lapply(object$objects, function(x, lmeObj){
			secondDeriv(x, lmeObj)
			}, lmeObj)
}

secondDeriv.conslmeObject <-
  function(object, lmeObj, ...){
	respSt <- respStruct(object, lmeObj)
	W.inv <- Var(respSt, T)
	X <- model.matrix(respSt)
	res <- residuals(respSt)
# 	U.inv <- lapply(W.inv, function(x) t(attr(x, "cholvar.inv")))
	W.dot <- pdDeriv(respSt)
# 	class(U.inv) <- "groupsList"
	W.inv.X <- Crossprod(W.inv, X)
	p <- ncol(X[[1]])
	if(ncol(object$varFix) == p)
 		secDer.fixed <- -solve(object$varFix)
		else
			secDer.fixed <- -Sum(Crossprod(X, W.inv.X))
	W.dot.W.inv <- Crossprod(W.dot, W.inv) #will be used transposed
	W.inv.res <- Crossprod(W.inv, res)
	W.inv.W.dot.W.inv.res <- Crossprod(W.dot.W.inv, W.inv.res)
	secDer.fixed.var <- -Sum(Crossprod(X, W.inv.W.dot.W.inv.res))
	secDer.fixed.var <- matrix(c(secDer.fixed.var), ncol = dim(secDer.fixed.var)[ 3 ])
	secDer.var <- lmeObjectHessian(object)
	nams <- dimnames(W.dot[[ 1 ]])[[3]]
	dimnames(secDer.fixed.var) <- list(dimnames(X[[ 1 ]])[[2]], nams)
	dimnames(secDer.var) <- list(nams, nams)
	A <- cbind(secDer.fixed, secDer.fixed.var)
	B <- cbind(t(secDer.fixed.var), secDer.var)
	secDer <- rbind(A, B)
	attr(secDer, "secDer.fixed") <- secDer.fixed
	attr(secDer, "secDer.fixed.var") <- secDer.fixed.var
	attr(secDer, "secDer.var") <- secDer.var
	attr(secDer, "W.inv.res") <- W.inv.res
	attr(secDer, "W.inv.W.dot.W.inv.res") <- W.inv.W.dot.W.inv.res
	attr(secDer, "W.inv.X") <- W.inv.X
	attr(secDer, "variance estimates") <- Coef(object, F)
	class(secDer) <- c("secondDeriv")
	secDer
}

determinant.secondDeriv <-
  function(x, logarithm = TRUE, ...){
	det(x[ , ])
}
