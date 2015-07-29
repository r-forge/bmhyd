\name{AttachHybridsToDonor}
\alias{BMhyd}
\title{
	Attach hybrids to the donor species
	}

\description{
	This function attachs the hybrid to the donor species.  
	}

\usage{
	AttachHybridsToDonor(phy, flow, suffix="_DUPLICATE")
    }
    
\arguments{
  \item{phy}{a phylogeny of phylo class.}
  \item{flow}{structure of gene flow.}
  \item{suffix}{a duplicated species that used for hybrid information in advanced.}
  }
  
\details{
	This function attachs hybrid to the donor and return the new phylogeny. It first identifies the flow clades using the function \code{LumpIntoClades}. Thenfor each clade containingthe flow, it searches the donor taxa and then attach the hybrid into the clade ancestor. For empirical analysis, users shall have the correct format for the flow strcutre to ensure the hybrid information. 
	}
	 
\value{
  a merged tree that identifies the hybrid species from the original tree.
  }

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2015. \emph{Trait evolution on phylogenetic networks}. Submitted. 
	}

\examples{
	#set the number of non hybrid species
	ntax.nonhybrid<-4
	#set the number of hybrid species
	ntax.hybrid<-1
	#simulate a network
	network<-SimulateNetwork(ntax.nonhybrid=ntax.nonhybrid, ntax.hybrid=ntax.hybrid, flow.proportion=0.5, origin.type='clade', birth = 1, death = 0.5, sample.f = 0.5, tree.height = 1, allow.ghost=FALSE)
	#show the tree. 
	network$phy
	#show the flow information
	network$flow
 	#Attach hybrid species to the donor and show the new tree. 
 	#The name of the hybrid is displayed with suffix  "_DUPLICATE"
	AttachHybridsToDonor(network$phy, network$flow, suffix="_DUPLICATE")
	}