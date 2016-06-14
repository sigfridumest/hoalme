## file hoalme/R/Hessian.R, v 0.0-1 2009/11/18
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


##Hessian methods

Hessian <- 
  function(object)
	UseMethod("Hessian")

Hessian.lmeObject <-
  function(object){
	hess <- object$Hessian
	result <- hess$Hessian
	if(all(eigen(result)$values < 0)) attr(result, "Negative-definite") <- T #comes negative valued
		else attr(result, "Negative-definite") <- F
	attr(result, "logLik") <- hess$mean
	attr(result, "gradient") <- hess$gradient
	class(result) <- "Hessian"
	result
}

Hessian.conslmeObjects <-
  function(object){
	result <- lapply(object$objects, Hessian)
	class(result) <- "Hessian.conslmeObjects"
	result
}

print.Hessian <-
  function(x, ...){
	loglik <- attr(x, "logLik")
	grad <- attr(x, "gradient")
	nd <- attr(x, "Negative-definite")
	cat("\n   LogLikelihood:\n", loglik, "\n")
	cat("\n   Gradient:\n", grad, "\n")
	cat("\n   Hessian:\n")
	print(x[ , ])
	cat("\n   Negative-definite:\n", nd, "\n\n")
}

