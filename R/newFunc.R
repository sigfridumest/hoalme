## file hoalme/R/newFunc.R, v 0.0-1 2009/11/18
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


##general functions (invisible to users)

lwtr <-
  function(x, n){
	x <- as.vector(x)
	n <- as.integer(n)
	lower.t <- diag(n)
	lower.t[col(lower.t) <= row(lower.t)] <- x
	lower.t
}

is.even <-
 function(x){
	x <- as.integer(x)
	y <- floor(x/2)
	if(2*y == x) TRUE
		else FALSE
}

rawpar <-
 function(object){
	parameters <- object$unconsPar
	parameters
}
