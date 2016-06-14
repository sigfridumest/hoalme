## file hoalme/R/V.R, v 0.0-1 2009/11/18
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


##V methods

V <-
  function(object, ...)
	UseMethod("V")

V.default <-
  function(object, res, chol.var, chol.deriv, ...){
	X <- object
	grouped.args <- groupsList(X, res, chol.var, chol.deriv)
	result <- lapply(grouped.args,
		function(x){
			#extracting the relevant matrices
			X <- as.matrix(x[[ 1 ]])
			res <- as.matrix(x[[ 2 ]])
			chol.var <- x[[ 3 ]]
			chol.deriv <- x[[ 4 ]]
			chol.var <- lwtr(chol.var, nrow(res))

			cholvar.inv <- solve(chol.var)
			std.res <- cholvar.inv%*%res
			#performing the matrix multiplications
			v.var <- apply(chol.deriv, 2, function(chol.deriv, std.res, cholvar.inv){
								chol.deriv <- lwtr(chol.deriv, nrow(std.res))
								v <- chol.deriv%*%std.res
								v
							}, std.res, cholvar.inv)
			v.vectors <- cbind(X, v.var)
			attr(v.vectors, "cholvar.inv") <- cholvar.inv
			attr(v.vectors, "std.res") <- std.res
			v.vectors
		} )
	class(result) <- c("V", "groupsList")
	result
}

V.lmeObject <-
  function(object, ...){
	V(respStruct(object))
}

V.respStruct  <-
  function(object, ...){
	object$V
}

print.V <-
  function(x, ...){
	result <- lapply(x, function(x){
			attr(x, "std.res") <- NULL
			attr(x, "cholvar.inv") <- NULL
			x
		})
	print(result)
}
