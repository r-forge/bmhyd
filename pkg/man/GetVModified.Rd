\name{GetVModified}
\alias{BMhyd}
\title{
	 The variance covariance matrix for the network model
	}

\description{
	The GetVModified function returns the variance covariance function for the network model. 
	}

\usage{
	GetVModified(x, phy, flow, actual.params)
	}

\arguments{
  \item{x}{
	  vector of the parameter. 
 	} 
 \item{phy}{
 	a class of phylo tree
 	}
 \item{flow}{
 	a flow struture of gene flow
 	}
 \item {actual.params}{
 	The parameters that are related to the hybridzation including the hybrid vigor \eqn{\beta} and variation at the 	burst of hybridization \eqn{v_H}. 	
 	}
}

\details{
	The GetVModified functon calculate the variance covariance matrix under the network model.  The function first read tree in 'phylo' class and convert the variance covariance matrix use the Brownian motion model. It then uses the structure of the gene flow to identify the hybrid species, parents and the direction of the flow. The covariance of between the hybrid and non hybrid species is the weighted sum of the covariance from evolution along the tree plus evolution along the migration path.
	}

\value{
 	the variance covariance matrix will be returned and is a function of parameters including the hybrid vigor \eqn{\beta} , the variation at the burst of hybridization the \eqn{v_H}, the rate of evolution \eqn{\sigma^2} and the measurement error SE.
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
	network<-SimulateNetwork(ntax.nonhybrid=ntax.nonhybrid, ntax.hybrid=ntax.hybrid, flow.proportion=0.5, 		origin.type='clade', birth = 1, death = 0.5, sample.f = 0.5, tree.height = 1, allow.ghost=FALSE)
	sigma.sq <- 0.01;mu <- 1;SE <- 0
	GetVModified(c(sigma.sq,mu,SE), network$phy, network$flow, "bt")
}