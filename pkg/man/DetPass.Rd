\name{DetPass}
\alias{BMhyd}
\title{
 Determinant of the matrix 
}
\description{
Verify the determinant of the matrix to have the positive and finite determinant.
}

\usage{
DetPass(phy)
}
\arguments{
  \item{phy}{
 a phylogenetic tree of class 'phylo'
 } 
}
\details{
	The DetPass functon first calculate variance covariance matrix C for the tree phylo. It then calculates the determinant of the C, 0.0001*C, and 1000*C and verify whether all the three determinant are positive and finite. If yes, the C matrix is likely to be the well-conditioned.  
}
\value{
 \item{det.pass}{a TRUE/FALSE value}
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