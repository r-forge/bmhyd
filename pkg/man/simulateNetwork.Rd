\name{SimulateNetwork}
\alias{BMhyd}
\title{
	Simulate phylogenetic network
	}

\description{
	Simulate network with desired number of taxa and number of hybrid.
	}

\usage{
	SimulateNetwork(ntax.nonhybrid=100, ntax.hybrid=10, flow.proportion=0.5, origin.type=c("clade", "individual"), birth = 1, death = 	0.5, sample.f = 0.5, tree.height = 1, allow.ghost=FALSE)
	}

\arguments{
  \item{ntax.nonhybrid}{
	number of non hybrid taxa.
	}
	
  \item{network}{
	number of hybrid.
	}
	
  \item{flow.proportion}{
  	the gene flow proportion from the parent.
  	}
  
  \item{origin.type}{
  	the type where the hybrids were formed. clade: hybrid formed and then speciate; individial: hybrid species formed and evolved 		without speciation.
  	}
  
  \item{birth}{
  	birth rate for the tree speciation.
  	}
  
  \item{death}{
  	death rate (extinction) of the tree. 
  	}
  	
  \item{sample.f}{
  	sampling frequency.
  	}
  
  \item{tree.height}{
  	the height of the tree.
  	}
  	
  \item{allow.ghost}{
  	allows ghost lineage persists hybridization for a while and then goes extinct. Otherwise, hybridization event 	
  	will be between coeval with extant descedants.
  	}
	}

\details{	
	This function generates tree with relevant hybridization information (the structure of gene flow information). The gene flow must happen on its stem edge, which starts at shortest from root and goes up to longest from root. 	Gene flow cannot go back in time. In general, gene flow can go forward in time via ghost lineages. If no ghost lineages, then there must be temporal overlap between the donor and recipient lineages. So the tipward end of the donor edge must be later than the rootward end of the recipient edge. When the non of hybris is 1 and clade sampling, this will do individual sampling instead. 
	} 


\value{
	a list that contains
	 \item{phy}{a birth death tree with number of taxa of ntax.nonhybrid + ntax.hybrid}
 	\item{flow}{The flow structure of the gene flow}
 	}

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2015. \emph{Trait evolution on phylogenetic networks}. Submitted. 
	}
	
\examples{
	library(TreeSim)
	#set up the number for non hybrid species
	non.hybrid <- 3
	#set up the number for hybrid species
	hybrid <- 1
	#set up the gene flow porportion
	flow<- 0.5
	#set pu the hybridization type to original
	origins<- "individual"
	#start to simulate the network
	network<-SimulateNetwork(ntax.nonhybrid=non.hybrid, ntax.hybrid=hybrid, flow.proportion=flow, origin.type="individual", birth = 1, death = 0.5, sample.f = 0.5, tree.height = 1, allow.ghost=FALSE)
	#print out the result
	(network)
	}