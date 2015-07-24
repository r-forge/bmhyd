\name{LumpIntoClade}
\alias{BMhyd}
\title{
Lump into clade
}
\description{
This function Lump the .....
}

\usage{
LumpIntoClades(phy,flow)
\arguments{
  \item{phy}{a phylogeny of phylo class}
  \item{flow}{structure of gene flow}
   }
   
\details{
 The function returns the donor clade, recipient clade, time from root donor, time from root recipient as a dataframe.   
    } 

\value{
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
	LumpIntoClades(network$phy, network$flow) 
	}