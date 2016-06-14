## file hoalme/R/cholDeriv.R, v 0.0-1 2009/11/18
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


##cholDeriv generics and methods

cholDeriv <-
  function(x, y, ...)
	UseMethod("cholDeriv")

cholDeriv.default <-
  function(x, y, ...){
	n <- nrow(x)
	d <- dim(y)[ 3 ]
	x. <- x[ row(x) <= col(x) ]
	y. <- c(apply(y, 3, function(x) x[ row(x) <= col(x) ]))
 	matrices <- .C( 'cholderiv', as.double(x.), as.double(y.), as.integer(n), as.integer(d) )
	x[ row(x) <= col(x) ] <- matrices[[ 1 ]]
	x <- t(x)
	x <- x[ row(x) >= col(x) ]
	matrices[[ 2 ]] <- matrix(matrices[[ 2 ]], ncol = d)
	deriv <- apply(matrices[[ 2 ]], 2, function(x, n){
				result <- diag(n)
				result[row(result) <= col(result)] <- x
				result <- t(result)
				result[row(result) >= col(result)]
			}, n)
	result <- matrix( c(x, deriv), ncol = d +1, nrow = length(x.) )
	dimnames(result) <- list(NULL, c("respVar", dimnames(y)[[ 3 ]]))
	result
}

