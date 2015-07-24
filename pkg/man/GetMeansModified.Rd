\name{GetMeanModified}
\alias{BMhyd}
\title{
 	The species means
	}
	
\description{
	The GetMeansModified function returns the species mean for the network model. 
	}

\usage{
	GetMeansModified(x, phy, flow, actual.params)
	}
	
\arguments{
  \item{x}{
   vector of the parameters: including rate of evolution \eqn{\sigma^2}, overall mean \eqn{\mu}, measurement error SE, hybrid vigor \eqn{\beta}, variation at the burst of hybridization \eqn{v_H}.  
 	} 
 \item{phy}{
 	a class of phylo tree
 	}
 \item{flow}{
 	a flow struture of gene flow
 	}
 \item {actual.params}{
 	The parameters that are related to the hybridzation including the hybrid vigor \eqn{\beta} and variation at the burst of hybridization \eqn{v_H}. 
 	}
	}

\details{
	The GetMeansModified functon return the means for each species in vector form. For the non hybrid species, it returns as a parameter \eqn{\mu} and for hybrid species it returns \eqn{\mu+ \log \beta}.   
	}
	
\value{
	the species mean
 	}

\author{
	Brian O'Meara <omeara.brian@gmail.com>, Dwueng-Chwuan Jhwueng <dcjhwueng@fcu.edu.tw>
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2015. \emph{Trait evolution on phylogenetic networks}. Submitted. 
	}
	
\examples{
	ntax.nonhybrid<-4	
	ntax.hybrid<-1
	network<-SimulateNetwork(ntax.nonhybrid=ntax.nonhybrid, ntax.hybrid=ntax.hybrid, flow.proportion=0.5, origin.type='clade', birth = 1, death = 0.5, sample.f = 0.5, tree.height = 1, allow.ghost=FALSE)
sigma.sq <- 0.01;mu <- 1;SE <- 0
	GetMeansModified(c(sigma.sq,mu,SE), network$phy, network$flow, "bt")
	}