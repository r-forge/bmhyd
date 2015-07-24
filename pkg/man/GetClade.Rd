\name{GetClade}
\alias{BMhyd}
\title{
	Get clade
	}

\description{
	Search and report the nodes of the tree that have the desired clade size.
	}

\usage{
	GetClade(phy, clade.size)
	}

\arguments{
  \item{phy}{an object of class 'phylo'}
  \item{clade.size}{size of the clade}
	}

\details{
	This function uses the 'phylo' class where the edges(ancestor-decedant relationship) of the tree are used to identify the interior node with desired number of descedent (clade.size) on the tips. It applies the findMRCA function in phytools to search the ancestor.   
	} 

\value{
	the interior node that has descedants of size clade.size on the tip.   
 	}

\author{
	Brian O'Meara <omeara.brian@gmail.com>, Dwueng-Chwuan Jhwueng <tony.jhwueng@gmail.com>
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2014. \emph{Studying trait evolution in hybrid species on phylogenetic networks}. Submitted. 
	}

\examples{
	library(ape)
	library(phytools)
 	phy<-rtree(3)
 	plot(phy) # plot the tree
 	clade.size<-2
 	GetClade(phy,clade.size) 
 	}