\name{SimulateNetwork}
\alias{BMhyd}
\title{
Simulate Network
}
\description{
Simulate network with desired number of taxa and number of hybrid.
}

\usage{
SimulateNetwork(ntax.nonhybrid=100, ntax.hybrid=10, flow.proportion=0.5, origin.type=c("clade", "individual"), birth = 1, death = 0.5, sample.f = 0.5, tree.height = 1, allow.ghost=FALSE)
}
\arguments{
  \item{ntax.nonhybrid}{
	number of non hybrid taxa.
	}
  \item{network}{
	number of hybrid.
	}
  \item{flow.proportion}{
  	the gene flow proportion from the parent
  }
  \item{origin.type}{
  	the type where the hybrids were formed. clade: hybrid formed and then speciate; individial: hybrid species formed and evolved without speciation.
  }
  
  \item{birth}{
  	birth rate for the tree speciation
  	}
  
  \item{death}{
  	death rate (extinction) of the tree 
  	}
  	
  \item{sample.f}{
  	sampling frequency
  	}
  
  \item{tree.height}{
  	the height of the tree
  	}
  	
  \item{allow.ghost}{
  	allow.ghost allows ghost lineage: something that persists for awhile, hybridizes, goes extinct. Otherwise, hybridization event 	
  	must between coeval with extant descedants
  }
}

\details{	
	The SimulateNetwork generates tree with relevant hybridization information(the structure flow information) where it takes a 		recipient clade. The gene flow must happen on its stem edge, which starts at shortest from root and goes up to longest from root. 	Gene flow cannot go back in time. In general, gene flow can go forward in time via ghost lineages. If no ghost lineages, then there must be temporal overlap between the donor and recipient lineages. So the tipward end of the donor edge must be later than the rootward end of the recipient edge. When the non of hybris is 1 and clade sampling, this will do individual sampling instead. 
	} 


\value{
	a list that contains
 \item{phy}{a birth death tree with number of taxa of ntax.nonhybrid + ntax.hybrid}
 \item{flow}{The flow structure of the gene flow}
 }

\author{
	Brian O'Meara <omeara.brian@gmail.com>, Dwueng-Chwuan Jhwueng <tony.jhwueng@gmail.com>
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2015. \emph{Trait evolution on phylogenetic networks}. Submitted. 
	}
	
\examples{
	library(TreeSim)
	##simulate a network of 4 taxa where 1 taxa is the hybrid
	non.hybrid <- 3
	hybrid <- 1
	flow<- 0.5
	origins<- "individual"
	network<-SimulateNetwork(ntax.nonhybrid=non.hybrid, ntax.hybrid=hybrid, flow.proportion=flow, origin.type="individual", birth = 1, death = 0.5, sample.f = 0.5, tree.height = 1, allow.ghost=FALSE)
(network)
	}