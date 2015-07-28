\name{AdjustForDet}
\alias{BMhyd}
\title{
 Adjust the phylogeny 
 
\description{
	A function adjusts the phylogeny when the determinant of the variance covariance matrix of the the phylogeny is negative. 
	}

\usage{
	AdjustForDet(phy, max.attempts=100)
	}

\arguments{
  \item{phy}{
	 	an object of class 'phylo'
 		} 
 		
 	\item{max.attepmts}{
 		the maximum numer of attempts to adjust the phylogeny. The default is 100.
 		}
	}

\details{
	This functon first calculates the determiant of the phylogeny using the function \code{DetPass}. When the determinant of the variance covariance matrix for the phylogeny is negative, this function will adjust the phylogeny by slightly lengthening the edge lengths of the tree. This procedure repeats until a good adjusted phylogeny is returned.}

\value{
 an object of class 'phylo'
 }

\author{
	Brian O'Meara
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2015. \emph{Studying trait evolution in hybrid species on phylogenetic networks}. Submitted. 
	}

\examples{
	library(ape)
	#simulate a tree of 5 taxa
	phy<-rtree(5)
	#check and adjust the phylogeny
	AdjustForDet(phy,max.attempts=100)
	}