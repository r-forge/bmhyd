\name{AttachHybridsToDonor}
\alias{BMhyd}
\title{
Attach hybrids to the donor
}
\description{
This function attach the hybrid to the donor
}

\usage{
AttachHybridsToDonor(phy, flow, suffix="_DUPLICATE")
\arguments{
  \item{phy}{a phylogeny of phylo class}
  \item{flow}{structure of gene flow}
  \item{suffix}{a duplicate species}
  }
\details{
	This function attachs hybrid to the donor. It first pull the clade and merge ...........
	} 
\value{
 \item{phy.merged}{the merged tree with network}
 }

\author{
	Brian O'Meara <omeara.brian@gmail.com>, Dwueng-Chwuan Jhwueng <tony.jhwueng@gmail.com>
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2015. \emph{Trait evolution on phylogenetic networks}. Submitted. 
	}

\examples{
	ntax.nonhybrid<-4
	ntax.hybrid<-1
	network<-SimulateNetwork(ntax.nonhybrid=ntax.nonhybrid, ntax.hybrid=ntax.hybrid, flow.proportion=0.5, origin.type='clade', birth = 1, death = 0.5, sample.f = 0.5, tree.height = 1, allow.ghost=FALSE)
	AttachHybridsToDonor(network$phy, network$flow, suffix="_DUPLICATE")
	}