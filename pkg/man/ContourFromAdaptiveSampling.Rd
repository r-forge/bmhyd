\name{ContourFromAdaptiveSampling}
\alias{BMhyd}
\title{
	Contour maps 
	}
	
\description{
	This function generates the countour maps for a pair of parameters of interest.  
	}

\usage{
	ContourFromAdaptiveSampling(sims, params.of.interest=NULL)
    }
    
\arguments{
  \item{sims}{a data frame with the simulation results}
  \item{params.of.interest}{name of parameter that is of interest to be plotted}
   }

\details{
	This function generates the countour plots for the parameter of interest. It uses the simulated data from generating the confidence interval under adapative sampling technique where the likelihood value and parameters values are stored. Then it uses the function \code{Plot_ConvexHull} to generate the contour maps. 
	} 

\value{
	It returns the plot of the contour map for a pair of parameters.
 	}

\author{
	Brian O'Meara
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2015. \emph{Studying trait evolution in hybrid species on phylogenetic networks}. Submitted. 
	}
	
\examples{
	#simulate a data frame for two parameters
	sims<-matrix(rnorm(300),ncol=3)
	sims<-as.data.frame(sims)
	colnames(sims)<-c("neglnL","param1","param2")
	#plot the contour map for the parameters
	ContourFromAdaptiveSampling(sims, params.of.interest=NULL)	
	}