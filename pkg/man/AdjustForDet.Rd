\name{AdjustForDet}
\alias{BMhyd}
\title{
 Adjust the phylogeny 
 
\description{
	Adjust the phylogeny when the determinant of the variance covariance matrix of the phylogeny is ill-conditioned. 
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
	The AdjustForDet functon first calculate the determiant of the phylogeny using the DetPass function. If the determinant of the phylogeny is difficult to calculate, the function will adjust the phylogeny by slightly lengthening the edge lengths of the tree. This procedure repeats until a good adjusted phylogeny is returned.}

\value{
 an object of class 'phylo'
 }

\author{
	Brian O'Meara <omeara.brian@gmail.com>,  Dwueng-Chwuan Jhwueng <tony.jhwueng@gmail.com>
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2015. \emph{Studying trait evolution in hybrid species on phylogenetic networks}. Submitted. 
	}

\examples{
	library(ape)
	phy<-rtree(5)
	AdjustForDet(phy,max.attempts=100)
	}