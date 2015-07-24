\name{DetPass}
\alias{BMhyd}
\title{
	Verification of the determinant of the matrix 
	}
	
\description{
	Calculate the determinant of the matrix(see details).
	}

\usage{
	DetPass(phy)
	}

\arguments{
  \item{phy}{
 an object of class 'phylo'
 } 
}

\details{
	The DetPass functon first calculates variance covariance matrix C for the tree. It then calculates the determinant of the C, 0.0001*C, and 1000*C and verifies whether all the three determinant are positive and finite. If yes, the C matrix is likely to be the well-conditioned.  
	}

\value{
 a TRUE/FALSE value for determining whether the matrix is likely to be well-conditioned. 
 }

\author{
	Brian O'Meara <omeara.brian@gmail.com>
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2015. \emph{Studying trait evolution in hybrid species on phylogenetic networks}. Submitted. 
	}

\examples{
	library(ape)
	##simulate a tree of 5 taxa
	phy<-rtree(5)
	DetPass(phy)
	}