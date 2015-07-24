\name{AdaptiveConfidenceIntervalSampling}
\alias{BMhyd}
\title{
Adaptive Confidence Interval Sampling
}
\description{
This function useadaptive confidence interval sampling technique to generate the confidence interval for parameter of interest.
}

\usage{
AdaptiveConfidenceIntervalSampling(par, fn, lower=-Inf, upper=Inf, desired.delta = 2, n.points=5000, verbose=TRUE, ...)
}
\arguments{
  \item{par}{parameter of interest for sampling in the model.}
  \item{fn}{the likelihood function}
  \item{lower}{the lower bound for the values allowed for sampling} 
  \item{upper}{the upper bound for the values allowed for sampling}
  \item{desired.delta}{a numeric value. The default is 2 meaning that the log likelihood shall no more than 2 unit away from he maximum.}
  \item{n.points}{number of required points for calculating likeilihood value}
}
\details{
The function starts with a set of parameters value, then the points are simulated by the GeneratedValues function. The likelihood is calculated. We set the width.ratio defined by by the ratio of the sampled ranged over the total range of the points to ensure the sampling is widely enough.  
} 

\value{
 \item{results}{a three column matrix where the first two columns represent the grid points of the parameter and the third column is  the calculated likelihood.}
 }

\author{
	Brian O'Meara <omeara.brian@gmail.com>, Dwueng-Chwuan Jhwueng <tony.jhwueng@gmail.com>
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2014. \emph{Studying trait evolution in hybrid species on phylogenetic networks}. Submitted. 
	}

\examples{
	ntax.nonhybrid<-4
	ntax.hybrid<-1
	network<-SimulateNetwork(ntax.nonhybrid=ntax.nonhybrid, ntax.hybrid=ntax.hybrid, flow.proportion=0.5, origin.type='clade', birth = 1, death = 0.5, sample.f = 0.5, tree.height = 1, allow.ghost=FALSE)

	tips<-rnorm(ntax.nonhybrid+ntax.hybrid)
	names(tips)<-paste("t",(1:(ntax.nonhybrid+ntax.hybrid)),sep="")
		
	free.parameters<-rep(TRUE, 5)
	names(free.parameters) <- c("sigma.sq", "mu", "bt", "vh", "SE")

	par<-c(0.01,1,1,0,0)

	interval.results <- AdaptiveConfidenceIntervalSampling(par, fn=CalculateLikelihood, lower=c(0, -Inf, 0, 0, 0)[which(free.parameters)], data=tips, phy=network$phy, flow=network$flow, actual.params=free.parameters[which(free.parameters)], allow.extrapolation=allow.extrapolation)

	interval.results.in <- interval.results[which(interval.results[,1]-min(interval.results[,1])<=2),]
	interval.results.out <- interval.results[which(interval.results[,1]-min(interval.results[,1])>2),]		
	}