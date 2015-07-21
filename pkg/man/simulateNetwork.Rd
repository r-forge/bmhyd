\name{SimulateNetwork}
\alias{BMhyd}
\title{
Simulate Network
}
\description{
Randomly simulate network with desired number of taxa and number of hybrid.
}

\usage{
SimulateNetwork(ntax.nonhybrid=100, ntax.hybrid=10, flow.proportion=0.5, origin.type=c("clade", "individual"), birth = 1, death = 0.5, sample.f = 0.5, tree.height = 1, allow.ghost=FALSE)
}
\arguments{
  \item{ntax.nonhybrid}{
number of non hybrid taxa (see Details).
}
  \item{network}{
number of hybrid  (see Details).
}
  \item{flow.proportion}{
  	the flow proportion from the parent
  }
  \item{origin.type}{
  	the type where the hybrids were formed. clade: hybrid formed and then speciate; individial: hybrid form without speciation.
  }
  \item{birth}{
  	birth rate for the tree speciation
  }
  \item{death}{
  	death rate (extinction)
  }
  \item{sample.f}{
  	sampling frequency
  }
  
  \item{tree.height}{
  	the height of the tree
  }
  \item{allow.ghost}{
  	allow.ghost allows ghost lineage: something that persists for awhile, hybridizes, goes extinct. Otherwise, hybridization event must between coeval with extant descedants
  }
}
\details{
	
	It take a recipient clade. The gene flow must happen on its stem edge, which starts at shortest from root and goes up to longest from root. Gene flow cannot go back in time.
} 

\value{
 \item{phy}{a birth death tree with number of taxa of ntax.nonhybrid + ntax.hybrid}
 \item{flow}{The flow structure of the gene flow}
 }

\author{
Brian O'Meara <omeara.brian@gmail.com>
}

\references{
Jhwueng D.C. and O'Meara B.C. 2014. \emph{Studying trait evolution in hybrid species on phylogenetic networks}. Submitted. 
}
\examples{
##simulate a network of 4 taxa where 1 taxa is the hybrid
non.hybrid <- c(3)
hybrid <- c(1)
flow<- 0.5
origins<- c("clade")
sigma.sq <- c(0.01)
mu <- c(1)
bt  <- c(0.1)
vh <- c(1)
SE <- c(0.5)
tree.height<- c(50)
nreps <- 50
SimulateNetwork(ntax.nonhybrid=non.hybrid, ntax.hybrid=hybrid, flow.proportion=flow.vector, origin.type= origin, sigma.sq=sigma.sq, mu=mu, bt=bt, vh=vh, SE=SE, tree.height=tree.height.v, allow.ghost=FALSE)
}