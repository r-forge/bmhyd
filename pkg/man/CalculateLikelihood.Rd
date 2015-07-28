\name{CalculateLikelihood}
\alias{BMhyd}
\title{
	Calculate the likelihood value for the model
	}
	
\description{
	This function calculates the likelihood value for the model
	}

\usage{
	CalculateLikelihood(x, data, phy, flow, actual.params, precision=2, proportion.mix.with.diag=0, allow.extrapolation=FALSE) 
	}
	
\arguments{
  \item{x}{parameter of interest}
  \item{data}{the trait values}
  \item{phy}{a tree of phylo class}
  \item{flow}{strcuture of gene flow}
  \item{actual.params}{the free parameters for hybridzation}
  \item{precision}{a number to verify the condition of the variance covariation for the network model} 
  \item{proportion.mix.with.diag}{the proportion value that applies to the diagonal of the vcv matrix}
  \item{allow.extrapolation}{a TRUE/FALSE argument}
	}
	
\details{
	This function calculates the likelihood value. As described in the argument, the \eqn{x} contains the parameter of interests including over all mean \eqn{\mu}, the rate of evolution \eqn{\sigma^2}, the measurement error SE and the hybrid vigor \eqn{\beta} and the variation at the burst of hybridization \eqn{v_H}. Names of the taxa in data vector and the phylogeny must be match for further analysis, otherwise it would terminate immediately. The negative log likelihood function is calculated and a reasonable likelihood value is retruned by checking the variance covariance matrix where a precision number is used for verifying the log conditioned number of the variance covaraince of the network model. When the matrix is ill conditioned, we appraopriately adjust the matrix by the shrink the value of the off diagonal matrix using a proportion array. Then the likelihood is calculated using a spline approximation.  
	} 

\value{
	negative log likelihood value
	 }

\author{
	Brian O'Meara
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2014. \emph{Studying trait evolution in hybrid species on phylogenetic networks}. Submitted. 
	}
	
\examples{
	#set the number of hybrid 
	ntax.nonhybrid<-5	
	#set the number of hybrid
	ntax.hybrid<-1
	#simulate a network
	network<-SimulateNetwork(ntax.nonhybrid=ntax.nonhybrid, ntax.hybrid=ntax.hybrid, flow.proportion=0.5, origin.type='clade', birth = 1, death = 0.5, sample.f = 0.5, tree.height = 1, allow.ghost=FALSE)
sigma.sq <- 0.01;mu <- 1;SE <- 0
	#simulatedata
	data<-rnorm(ntax.nonhybrid+ntax.hybrid)
	names(data)<-paste("t",(1:(ntax.nonhybrid+ntax.hybrid)),sep="")
	#calculate the likelihood value
	CalculateLikelihood(c(sigma.sq,mu,SE), data, network$phy, network$flow, "vh", precision=2, proportion.mix.with.diag=0, allow.extrapolation=TRUE) 
	}