\name{AdjustForDet}
\alias{BMhyd}
\title{
 Adjust the phylogeny 
 
\description{
Adjust the phylogeny due to the difficulty of calculation of the determinant. 
}

\usage{
AdjustForDet(phy, max.attempts=100)
}
\arguments{
  \item{phy}{
 a phylogenetic tree of class 'phylo'
 } 
 \item{max.attepmts}{
 	the maximum numer of attempts to adjust the phylogeny. The default is 100.
 }
}
\details{
	The AdjustForDet functon first calculate the determiant of the phylogeny using the DetPass function, if the determinant of the phylogeny is difficult to calculate, the function will adjust the phylogeny by slightly lengthening the edge lengths of the tree. It repeats this until a phylogeny with good determinant is returned.}
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
phy<-rtree(5)
AdjustForDet(phy,max.attempts=100)
}