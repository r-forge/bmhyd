\name{GetAncestor}
\alias{BMhyd}
\title{
	Ancestor identification
	}

\description{
	Identify the ancestor node for the given descedent note.
	}

\usage{
	GetAncestor(phy, node)
	}

\arguments{
  \item{phy}{an object of class 'phylo'}
  \item{node}{descedant node}
	}

\details{
	Given the descedent node for the tree, the function will return its ancestor. 
	} 

\value{
	the ancestor node 
 	}

\author{
	Brian O'Meara <omeara.brian@gmail.com>, Dwueng-Chwuan Jhwueng <tony.jhwueng@gmail.com>
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2015. \emph{Studying trait evolution in hybrid species on phylogenetic networks}. Submitted. 
	}

\examples{
	library(ape)
 	phy<-rtree(3)
 	plot(phy) # plot the tree
 	node<-1
 	GetAncestor(phy,node) 
 	}