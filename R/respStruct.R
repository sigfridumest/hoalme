## file hoalme/R/respStruct.R, v 0.0-1 2009/11/18
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


##respStruct generics and methods, and methods related to 'respStruct' class

respStruct <-
function(object, ...)
	UseMethod("respStruct")

respStruct.lmeObject <-
  function(object, ...){

	if(!is.null(object$respStruct)) return(object$respStruct) #avoids computations
	
	sigmasq <- object$sigma^2
	lmeSt <- object$modelStruct
  	names.lmeSt <- names(Coef(lmeSt))

	##getting the mean effects (random and fixed), together with grouping structure in a data.frame
	resp.frame <- model.matrix(object) #object of class grps.model.frame
	grps <- as.data.frame(resp.frame[ , attr(resp.frame, "gcols")])
	Q <- attr(resp.frame, "Q") #number of grouping levels

	##grouping the main effects frame at the highest grouping level
  	meanEff <- split( resp.frame, factor(grps[ , Q], levels = unique(grps[ , Q])) )
	class(meanEff) <- c("groupsList", "grps.model.frame")

	###
	###Constructing variance and variance derivatives, and storing them in a list representing the highest level of grouping
	###

	##getting the variance
	respVar <- Var(lmeSt, meanEff)
	#for now, not assuming any other structure than reStruct
	respVar <- lapply(respVar, function(x, sigmasq){
				diag(x) <- diag(x) + 1
				x*sigmasq
			}, sigmasq)
	
	##getting the variance derivatives at each structure and level
	varDeriv <- pdDeriv(lmeSt, meanEff)
# print(varDeriv)
	#>>>>>> not including pdBlock structures
	reStDeriv <- varDeriv$reStruct
 	
 	#getting the residual variance derivative
	varDeriv$residual <- lapply(reStDeriv, function(x, Q){
			x <- x[[ Q ]]
			dimax <- length(dim(x))
			if(dimax < 3) result <- diag(nrow(x))
				else  result <- diag(nrow(x[ , , 1]))
			result
		}, Q)
	
	##getting response variance derivatives (constructing block diagonal matrices for the lower levels)
	if( Q > 1){ 		
		#getting the derivative variances matrices as block diagonal matrices
		deriv.matrices <- lapply(varDeriv$reStruct, function(x, Q, nams){
		##( grand lapply) discloses the ***lower levels @ each group*** of the highest level
			result <- list()
			result[[ Q ]] <- x[[ Q ]]
			for(i in 1:(Q -1)){
				result[[ i ]] <- bdMatrices(x[[ i ]])
			}
			names(result) <- nams
			result
		}, Q, names(lmeSt))

		varDeriv <- groupsList( deriv.matrices, varDeriv$residual, nams = c("reStruct", "residual") )
	}
 		else varDeriv <- groupsList(reStDeriv, varDeriv$residual, nams = c(names(lmeSt), "residual") )
	##*

	#putting everything in a list of arrays
	listn <- length(varDeriv[[ 1 ]])
	deriv.unlist <- lapply(varDeriv, function(x, listn) {
			n <- nrow(x[[ listn ]])	
			x <- unlist(x)
			d <- as.integer( length(x)/(n*n) )
			x <- array( x, c(n, n, d) )
			x
			}, listn )
	var.deriv <- groupsList(respVar, deriv.unlist, nams = c("respVar", "derivatives"))

	rm(deriv.unlist)
 
	## assuming no other Structures (**)
	nams <- c(names.lmeSt, "residual") #Coef has names!!! rules....

	var.deriv <- lapply(var.deriv, function(x, nams){
				dimnames(x[[ 2 ]]) <- list(NULL, NULL, nams)
				x
			}, nams)
	rm(nams)
	
	##cholesky factorization and derivatives
	chol.objects <- lapply(var.deriv, function(x){
			var <- x[[ 1 ]]
			der <- x[[ 2 ]]
			nams <- names(der)
			cholderiv <- cholDeriv(var, der)
			cholderiv
		})
	cholVar <- lapply(chol.objects, function(x) x[ , 1]) #cholesky decomposition of the response variance
	cholDeriv <- lapply(chol.objects, function(x) x[ , -1]) #cholesky decomposition of the derivatives of the response variance

	##getting the V vectors
	X <- model.matrix(meanEff)
	res <- residuals(meanEff)

	V.vectors <- V(X, res, cholVar, cholDeriv)

 	V.vectors <- lapply(V.vectors, function(x, X) {
 			xnams <- dimnames(X[[ 1 ]])[[ 2 ]]
 			dimnames(x) <- list(NULL, c(xnams, names.lmeSt, "residual"))
 			x
  		}, X )
	class(V.vectors ) <- c("V", "groupsList")

	##getting the inverse of the response variance
	respVar.inv <- lapply(V.vectors, function(x) {
			chinv <- attr(x, "cholvar.inv")
			varinv <- crossprod(chinv)
			attr(varinv, "cholvar.inv") <- chinv
			varinv
		} )
	class(respVar.inv) <- "groupsList"

	stdRes <- lapply(V.vectors, function(x) attr(x, "std.res") )
	class(stdRes) <- "groupsList"

	var.deriv <- lapply(var.deriv, function(x) x[[ - 1 ]])

	## Second Derivatives
	##Will change with other variance structures *()
	##**
	varsecDer <- NULL
	##*
	
	result <- list(meanEffects = meanEff,
		respVar = respVar,
		varDeriv = var.deriv,
		cholVar = cholVar,
		cholDeriv = cholDeriv,
		V = V.vectors,
		respVar.inv = respVar.inv,
		Residuals = res,
		stdResiduals = stdRes,
		varsecDer = varsecDer,
		sigmasq = sigmasq)

	class(result) <- "respStruct"
	object$respStruct <- result

	result
}

respStruct.conslmeObject <-
  function(object, Struct = NULL, ...){
	#some needed objects
	if(!is.null(object$respStruct)) return(object$respStruct) #avoids computations
	
	sigmasq <- object$sigma^2
	lmeSt <- object$modelStruct
  	names.lmeSt <- names(Coef(lmeSt))

	##getting the mean effects (random and fixed), together with grouping structure in a data.frame
	resp.frame <- model.matrix(object) #object of class grps.model.frame
	grps <- as.data.frame(resp.frame[ , attr(resp.frame, "gcols")])
	Q <- attr(resp.frame, "Q") #number of grouping levels

	##grouping the main effects frame at the highest grouping level
  	meanEff <- split( resp.frame, factor(grps[ , Q], levels = unique(grps[ , Q])) )
	class(meanEff) <- c("groupsList", "grps.model.frame")

	###
	###Constructing variance and variance derivatives, and storing them in a list representing the highest level of grouping
	###

	##getting the variance
	respVar <- Var(lmeSt, meanEff)
	#for now, not assuming any other structure than reStruct
	respVar <- lapply(respVar, function(x, sigmasq){
				diag(x) <- diag(x) + 1
				x*sigmasq
			}, sigmasq)
	
	#################################################################################################3
	 varDeriv <- pdDeriv(Struct)
	
	##################################################################################################
	
	res <- residuals(meanEff)

	##getting the inverse of the response variance

	respVar.inv <- lapply(respVar, function(x) {
			ch <- chol(x)
			varinv <- chol2inv(ch)
			attr(varinv, "cholvar.inv") <- t(solve(ch))
			varinv
		} )
	class(respVar.inv) <- "groupsList"

	##Will change with other variance structures *()
	##**
	varsecDer <- NULL
	##*

	result <- list(meanEffects = meanEff,
		respVar = respVar,
		varDeriv = varDeriv,
		respVar.inv = respVar.inv,
		Residuals = res,
		varsecDer = varsecDer,
		sigmasq = sigmasq)

	class(result) <- "respStruct"
# 	object$respStruct <- result

	result
}

model.matrix.respStruct <-
 function(object, random = FALSE, ...){
	meanEff <- object$meanEffects
	if(random){
		result <- lapply(meanEff, function(x){
			Zcols <- attr(x, "Zcols")
			result <- x[ , Zcols$cols]
			result <- as.matrix(result) 
			Zcols$cols <- NULL
			attr(result, "Zcols") <- Zcols
			result
		} )
	}
	else{
		result <- lapply(meanEff, function(x){
			cols <- attr(x, "Xcols")
			X <- x[ , cols]
			result <- as.matrix(X)
			if(is.vector(X) && unique(X) == 1){
				xname <- "(Intercept)"
				dimnames(result) <- list(NULL, xname)
			}
			result
		} )
	}
	class(result) <- "groupsList"
	result
}

residuals.respStruct <-
  function(object, std = FALSE, ...){
	if(std) result <- object$stdResiduals
		else
			result <- object$Residuals

	class(result) <- "groupsList"
	result
}


# getResponse.respStruct <-
#    function(object){
# 	x <- object$meanEffects
#  	result <- lapply(x, function(x) x$y)
# 	class(result) <- "groupsList"
# 	result
# }
# 
# getGroups.respStruct <-
#    function(object){
# 	x <- object$meanEffects
#  	result <- lapply(x, function(x){
# 		gcols <- attr(x, "gcols")
# 		result <- x[ , gcols]
# 		attr(result, "Q") <- attr(x, "Q")
# 		result
# 	} )
# 	class(result) <- "groupsList"
# 	result
# }
