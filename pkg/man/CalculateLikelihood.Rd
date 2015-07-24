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
  \item{n.points}{number of required simulated likeilihood value}
}
\details{
The function starts to sample by proposing a set of parameters and estimating the likelihood. the likelihood value is then plotted versus each of the parameters used in that set in a different subplot as a gray or black dots. Dots in black represent the desired likelihood value and the width of the adaptive confidence intervals are measured by the left most and right model black dots.
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
ntax.hybrid<-2
network<-SimulateNetwork(ntax.nonhybrid=ntax.nonhybrid, ntax.hybrid=ntax.hybrid, flow.proportion=0.5, origin.type='clade', birth = 1, death = 0.5, sample.f = 0.5, tree.height = 1, allow.ghost=FALSE)
sigma.sq <- 0.01;mu <- 1;SE <- 0

data<-rnorm(ntax.nonhybrid+ntax.hybrid)
names(tips)<-paste("t",(1:(ntax.nonhybrid+ntax.hybrid)),sep="")

CalculateLikelihood(c(sigma.sq,mu,SE), data, network$phy, network$flow, "vh", precision=2, proportion.mix.with.diag=0, allow.extrapolation=TRUE) 
}