## file hoalme/R/groupsList.R, v 0.0-1 2009/11/18
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



groupsList <-
  function(object, ...)
	UseMethod("groupsList")

groupsList.default <- 
  function(object, ..., nams = NULL){
	x <- list(object, ...)
	if(length(x) == 1) {
		y <- x[[ 1 ]]
		class(y) <- "groupsList"
		return(y)
	}
	#initial checking
	if( !all(unlist(lapply(x, function(x) is.list(x)))) ) stop("Not all arguments are lists")
	lengths <- unlist(lapply(x, function(x) length(x)))
	lengths.1 <- lengths[ 1 ]
	if( any(ifelse(lengths != lengths.1, TRUE, FALSE)) )  stop("Lists not of the same length.")

	#getting the names
	if(is.null(nams)){
		mc <- match.call(expand.dots = F)
		nams <- as.character(mc[[ 2 ]])
	}
	#constructing indexes
	llen <- length(x)
	indexes <- matrix(1:(lengths.1*llen), ncol = llen)

	#reindexing
	result <- apply(indexes, 1, function(indexes, lists, nams) {
					y <- unlist(lists, FALSE)
					y <- y[indexes]
  					names(y) <- nams
					y
				}, lists = x, nams)
	class(result) <- "groupsList"
	result
}

solve.groupsList <-
  function(a, b = NULL, ...){
	result <- lapply(a, solve)
	class(result) <- "groupsList"
	result
}

names.groupsList <-
  function(x){
	names(x[[ 1 ]])
}

"names<-.groupsList" <-
  function(x, value){
	result <- lapply(x, function(x) { names(x) <- value; x } )
	class(result) <- "groupsList"
	result
}

bdMatrices <-
  function(x, ...)
	UseMethod("bdMatrices")

bdMatrices.default <-
  function(x, ...){
	##generating the matrix to contain the submatrices and the block diagonal entries
	dims <- unlist(lapply(x, function(x) dim(x)[1]))
	dims <- rep(cumsum(dims), dims)
	last.mat <- length(dims)
	n <- dims[last.mat]
	bigmat <- diag(dims[last.mat])
	entries <- row(bigmat) <= matrix(dims, ncol = n, nrow = n, byrow = T) & col(bigmat) <= dims
	dimx <- dim(x[[ 1 ]])

	if(length(dimx) == 3){ 
		bigarray <- array(bigmat, c(dim(bigmat), dimx[3]))

		##Getting the matrices across the list and not in each list(per matrix): l1: mat1, ..., matm  ...  ln:mat1, ..., matm
		for(i in 1:dimx[3]){ 
				result <- lapply(x, function(x, i) x[ , , i ], i  )
				bigarray[ , , i ][entries] <- unlist(result)
		}

	}
	else{ #when only one matrix is present as an element of the list 
		bigarray <- bigmat
		result <- lapply(x, function(x) x)
		bigarray[entries] <- unlist(result)
	}
	bigarray
}


Crossprod <-
  function(object, ...)
	UseMethod("Crossprod")

Crossprod.default <-
  function(object, y = NULL, ...){
	crossprod(object, y)
}

Crossprod.list <-
  function(object, y = NULL, ...){
	if(is.null(y)) result <- lapply(object, Crossprod)
		else{
			object <- groupsList(object, y)
			result <- Crossprod(object)
		}
	class(result) <- "groupsList"
	result
}

Crossprod.groupsList <-
  function(object, ...){
	object <- groupsList(object, ...)
	result <- lapply(object, function(object){
				if(is.list(object)){
					if(length(object) > 2) stop("No methods defined for more than two lists per group.")
						else{
							obj1 <- object[[ 1 ]]
							obj2 <- object[[ 2 ]]
							if(is.vector(obj1)) obj1 <- as.matrix(obj1)
							if(is.vector(obj2)) obj2 <- as.matrix(obj2)
							result <- Crossprod(obj1, obj2)
						}
				}
				else 
					result <- Crossprod(object)
			result
			})
	class(result) <- "groupsList"
	result
}


Crossprod.matrix <-
 function(object, y = NULL, ...){
	if(is.null(y)){
		result <- crossprod(object)
	}
	else{
		if(is.matrix(y)) result <- crossprod(object, y)
			else{
				if(is.array(y)){
					result <- c(apply(y, 3, function(y, object) crossprod(object, y), object))
					result <- array(result, c(dim(object)[2], dim(y)[2:3]))
					dimnames(result) <- list(NULL, NULL, dimnames(y)[[ 3 ]])
				}
				else stop("No methods defined for other objects than matrices or arrays.")
		}
	}
	result
}


Crossprod.array <-
 function(object, y = NULL, ...){
	nams <- dimnames(object)[[ 3 ]]
	dimobj <- dim(object)
	if(is.null(y)){
				result <- c(apply(object, 3, Crossprod))
				result <- array(result, c(dimobj[2], dimobj[2], dimobj[3]))
				dimnames(result) <- list(NULL, NULL, nams)
	}
	else{
		dimy <- dim(y)
		if(is.matrix(y)){
			result <- c(apply(object, 3, function(object, y) crossprod(object, y), y))
			result <- array(result, c(dimobj[ 2 ], dimy[ 2 ], dimobj[ 3 ]))
			dimnames(result) <- list(NULL, NULL, nams)
		}
		else{
			if(is.array(y)){
				result <- c( apply(object, 3, function(object, y){
							apply(y, 3, function(y, object) crossprod(object, y), object)
							}, y) )
				result <- array(result, c(dimobj[2 ], dimy[2], dimobj[ 3 ]*dimy[ 3 ]))
			}
			else
				stop("No methods defined for other objects than matrices or arrays.")
		}
	}
	result
}


Sum <-
  function(object, ...)
	UseMethod("Sum")

Sum.groupsList <-
 function(object, ...){
	Reduce("+", object)
}

Sum.list <-
 function(object, ...){
	Reduce("+", object)
}

# Trace <-
#   function(object, ...)
# 	UseMethod("Trace")
# 
# tr <-
#   function(x, y = NULL){
# 	if(is.null(y)){
# 		result <- sum(diag(x))
# 	}
# 	else{
# 		n <- as.integer(nrow(x))
# 		x <- as.double(x)
# 		y <- as.double(y)
# 		w <- as.double(0)
# 		result <- .C("trace", x, y, w, n)[[ 3 ]]
# 	}
# 	return(result)
# }
# 
# Trace.default <-
#   function(object, ...){
# 	lapply(object, tr)
# }
# 
# Trace.groupsList <-
#   function(object, y = NULL, t = FALSE, ...){
# 	if(t) object <- t(object)
# 
# 	if(is.null(y)){
# 		result <- lapply(object, function(object) Trace(object))
# 	}
# 	else{
# 		object <- groupsList(object, y)
# 		result <- lapply(object, function(object){
# 					if(is.list(object)){
# 						if(length(object) > 2) stop("No methods defined for more than two (inner) lists.")
# 						else{
# 							result <- Trace(object[[ 1 ]], object[[ 2 ]])
# 						}
# 					}
# 				})
# 	}
# 	class(result) <- "groupsList"
# 	result
# }
# 
# 
# Trace.matrix <-
#   function(obejct, y = NULL, ...){
# 	if(is.null(y)){
# 		result <- tr(obejct)
# 	}
# 	else{
# 		if(is.matrix(y)) result <- tr(object, y)
# 			else{
# 				if(is.array(y)){
# 					result <- apply(y, 3, function(y, object) tr(object, y), object)
# 					result <- matrix(result, nrow = 1)
# 				}
# 				else stop("No methods defined for other objects than matrices or arrays.")
# 			}
# 	}
# 	result
# }
# 
# Trace.array <-
#  function(object, y = NULL, ...){
# 	dimobj <- dim(object)
# 	if(is.null(y)){
# 		result <- apply(object, 3, tr)
# 		result <- matrix(result, nrow = 1)
# 	}
# 	else{
# 		dimy <- dim(y)
# 		if(is.matrix(y)){
# 			result <- apply(object, 3, function(object, y) tr(object, y), y)
# 			result <- matrix(result, nrow = 1)
# 		}
# 		else{
# 			if(is.array(y)){
# 				result <- apply(object, 3, function(object, y){
# 							result <- apply(y, 3, function(y, object){ tr(object, y)}, object)
# 							}, y)
# 				result <- matrix(result, nrow = 1)
# 			}			
# 			else
# 				stop("No methods defined for other objects than matrices or arrays.")
# 		}
# 	}
# 	result
# }

t.groupsList <-
function(x){
	result <- lapply(x, function(x){
			if(is.list(list) & length(x) > 1) stop("No methods for handling more than one (inner) list")
			else{
				if(is.vector(x)) result <- t(as.matrix(x))
				else{
					if(is.matrix(x)) result <- t(x)
					else{
						if(is.array(x)){
							result <- c(apply(x, 3, t)) 
							result <- array(result, dim(x)[c(2,1,3)])
						}
						else stop("No methods for handling objects of class: ' ", class(x), " '")
					}
				}
			}
		result
		})
	class(result) <- "groupsList"
	result
}

Tcrossprod <-
  function(object, ...)
	UseMethod("Tcrossprod")

Tcrossprod.default <-
  function(object, y = NULL){
	crossprod(object, y)
}

Tcrossprod.list <-
  function(object, y = NULL){
	if(is.null(y)) result <- lapply(object, Tcrossprod)
		else{
			object <- groupsList(object, y)
			result <- Tcrossprod(object)
		}
	class(result) <- "groupsList"
	result
}

Tcrossprod.groupsList <-
  function(object, ...){
	object <- groupsList(object, ...)
	result <- lapply(object, function(object){
				if(is.list(object)){
					if(length(object) > 2) stop("No methods defined for more than two lists per group.")
						else{
							obj1 <- object[[ 1 ]]
							obj2 <- object[[ 2 ]]
							if(is.vector(obj1)) obj1 <- as.matrix(obj1)
							if(is.vector(obj2)) obj2 <- as.matrix(obj2)
							result <- Tcrossprod(obj1, obj2)
						}
				}
				else 
					result <- Tcrossprod(object)
			result
			})
	class(result) <- "groupsList"
	result
}

Tcrossprod.matrix <-
 function(object, y = NULL){
	if(is.null(y)){
		result <- tcrossprod(object)
	}
	else{
		if(is.matrix(y)) result <- tcrossprod(object, y)
			else{
				if(is.array(y)){
					result <- c(apply(y, 3, function(y, object) tcrossprod(object, y), object))
					result <- array(result, c(dim(object)[2], dim(y)[2:3]))
					dimnames(result) <- list(NULL, NULL, dimnames(y)[[ 3 ]])
				}
				else stop("No methods defined for other objects than matrices or arrays.")
		}
	}
	result
}


Tcrossprod.array <-
 function(object, y = NULL){
	nams <- dimnames(object)[[ 3 ]]
	dimobj <- dim(object)
	if(is.null(y)){
				result <- c(apply(object, 3, Tcrossprod))
				result <- array(result, c(dimobj[2], dimobj[2], dimobj[3]))
				dimnames(result) <- list(NULL, NULL, nams)
	}
	else{
		dimy <- dim(y)
		if(is.matrix(y)){
			result <- c(apply(object, 3, function(object, y) tcrossprod(object, y), y))
			result <- array(result, c(dimobj[ 2 ], dimy[ 2 ], dimobj[ 3 ]))
			dimnames(result) <- list(NULL, NULL, nams)
		}
		else{
			if(is.array(y)){
				result <- c( apply(object, 3, function(object, y){
							apply(y, 3, function(y, object) tcrossprod(object, y), object)
							}, y) )
				result <- array(result, c(dimobj[2 ], dimy[2], dimobj[ 3 ]*dimy[ 3 ]))
			}
			else
				stop("No methods defined for other objects than matrices or arrays.")
		}
	}
	result
}
