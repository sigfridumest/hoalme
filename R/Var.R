## file hoalme/R/Var.R, v 0.0-1 2009/11/18
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



##Var methods

Var <-
  function(object, ...)
	UseMethod("Var")

Var.reStruct <-
  function(object, respFr, cholSt = NULL, m = length(object), sigmasq = 1, ...){

	if(is.null(cholSt)) cholSt <- chol(object) #getting cholesky factor of variance (reRtruct) if not given

	##extracting Z matrix at the current level
	grps <- as.data.frame(respFr[ , attr(respFr, "gcols")]) #groups
	grps.thislevel <- factor(grps[ , m], levels = unique(grps[ , m]))
	ngrps <- length(attr(grps.thislevel, "levels")) #how many groups at the current level

	if(ngrps!=1){
		result <- lapply(split(respFr, grps.thislevel), FUN = function(x, reSt, cholSt, m) Var(reSt, x, cholSt, m), reSt = object,
				cholSt = cholSt, m = m)
	}
	else {
		#Z at the current grouping factor
		Z <- as.data.frame(respFr[ , attr(respFr, "Zcols")$cols])
		starts <- attr(respFr, "Zcols")$starts
		ends <- attr(respFr, "Zcols")$ends
		Z.thislevel <- as.matrix(Z[ , starts[m]:ends[m] ])
		n <- nrow(Z.thislevel)

		##variance component at the current level
		if( inherits(object[[ m ]], "pdBlocked") ){
			chol.var.list <- cholSt[[ m ]]
			chol.var <- bdMatrices(chol.var.list)
			result <- tcrossprod(tcrossprod(Z.thislevel, chol.var))
		}
		else
			result <- tcrossprod(tcrossprod(Z.thislevel, cholSt[[ m ]]))

		##going for the variance components at lower levels
		if(m-1 > 0){
  			grps.factor <- factor(grps[ , m -1], levels = unique(grps[ , m -1]))
 			result.prev <- lapply(split(respFr, grps.factor), FUN = function(x, reSt, cholSt, m) Var(reSt, x, cholSt, m), reSt = object, cholSt = cholSt, m = m -1)
			dims <- unlist(lapply(result.prev, nrow))
			dims <- rep(cumsum(dims), dims)
			entries <- row(result) <= matrix(dims, ncol = n, nrow = n, byrow = T) & col(result) <= dims
			result[entries] <- result[entries] + unlist(result.prev)
		}
	}
	if(sigmasq!=1) 
		result <- lapply(result, function(x, sigmasq) x*sigmasq, sigmasq)
	result
}

Var.respStruct <-
  function(object, inverse = FALSE, ...){
	if(inverse){
		result <- object$respVar.inv
		class(result) <- "groupsList"
	}
	else{
		result <- object$respVar
		class(result) <- "groupsList"
	}
	result
}

Var.lmeStruct <-
  function(object, meanEff, sigmasq = 1, ...){
	##reStruct
	reSt <- object$reStruct
	cholreSt <- chol(reSt)
	result <- lapply(meanEff, function(x, reSt, cholreSt, m)
						Var(reSt, x, cholreSt, m), reSt, cholreSt, length(reSt))

	##corStruct and varStruct (missing)
	if(sigmasq !=1 ) 
		result <- lapply(result, function(x, sig2) x*sig2, sig2 = sigmasq)
	
	result
}

Var.lmeObject <-
  function(object, inverse = FALSE, ...){
  respSt <- respStruct(object)
  Var(respSt, inverse)
}