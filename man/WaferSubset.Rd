
\name{WaferSubset}
\alias{WaferSubset}
\title{Modeling of Analog MOS Circuits}
\description{
    The \code{WaferSubset} data frame has 200 rows and 4 columns, which is a
	 subset of the \code{Wafer} dataset from the \code{nlme} package by Bates
	 and Pinheiro.
}
\format{
    This data frame contains the following columns:
    \describe{
	\item{Wafer}{
	    a factor with levels
	    \code{1} 
	    \code{2} 
	    \code{3} 
	    \code{4} 
	    \code{5} 
	    \code{6} 
	    \code{7} 
	    \code{8} 
	    \code{9} 
	    \code{10} 
	}
	\item{Site}{
	    a factor with levels
	    \code{1} 
	    \code{2} 
	    \code{3} 
	    \code{4} 
	}
	\item{voltage}{
	    a numeric vector
	}
	\item{current}{
	    a numeric vector
	}
    }
}
\source{
    Pinheiro, J. C. and Bates, D. M. (2000), \emph{Mixed-Effects Models in S
	and S-PLUS}, Springer, New York.
}
\keyword{datasets}
