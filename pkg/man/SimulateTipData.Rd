\name{SimulateTipData}
\alias{SimulateTipData}
\title{
	Simulate data on a network
	}

\description{
	Simulate data on a network given its flow parameters
	}

\usage{
	SimulateTipData(phy, flow, params, suffix="_DUPLICATE")
	}

\arguments{
  \item{phy}{
	phylo object
	}
	
  \item{flow}{
	data.frame of the flow parameters
	}
	
  \item{params}{
  	named vector of BMhyd parameter values
  	}
  
  \item{suffix}{
	taxa are duplicated on the network; this suffix distinguishes them
	 }
}
 \details{	
	This function simulates tip data on a network, under a BMhyd model. 
	} 


\value{
	a named vector of tip values
 	}

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2015. \emph{Trait evolution on phylogenetic networks}. Submitted. 
	}
	
\examples{
#add example later. Perhaps using cichlid

	}