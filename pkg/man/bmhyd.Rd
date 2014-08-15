\name{BMhyd}
\alias{BMhyd}
\title{
Comparative method for hybridization (include gene flow)
}
\description{
Fit Brownian motion model of continuous character to investigate hybrid species through the hybrid vigor and hybrid variation. 
}

\usage{
BMhyd(data, phy, flow,...)
}
\arguments{
  \item{data}{
continuous trait data containing species information in vector format (see Details).
}
  \item{phy}{
a rooted phylogenetic tree in Newick format (see Details).
}
  \item{flow}{
  a data frame for parental information of the hybrids 
  }
}
\details{Function BMhyd fits likelihood models for continuous characters through incorporating phylogenetic network and comparative data. 
The full likelihood model includes four parameters: the ancenstral state \eqn{\mu}, the overall variation \eqn{\sigma^2}, 
the hybrid vigor \eqn{\beta}, and hybrid burst variation at formation \eqn{v_H}. We use four submodels that allows some fixed values of parameters 
and treat others as free parameters: (1) \eqn{\beta=1, v_H=0}; (2) \eqn{\beta} is free parameter, \eqn{v_H=0};   
(3) \eqn{\beta=1, v_H} is free parameter (4)\eqn{\beta} and \eqn{v_H} are free parameters. 

BMhyd fits the model through maximum likelihood technique. It returns MLEs and standard errors for above four models. The standard errors from adaptive sampling could be be optionally returned upon users request.  

As input, \pkg{BMhyd} requires (a) a phylogenetic tree with known topology and branch length and (b) a set of comparative data (c) a flow data that includes taxon names, proportion (m) inherited from the parents, and hybridization time (gene flow time). 

The variance covariance matrices can be poorly conditioned, making them hard to deal with numerically. If this is detected, we have an error checking and an approximation approach built in.
In each case, we convert the matrix to a weighted sum of matrices: the original one and the diagonals of the original one, which has the effect of decreasing the magnitude of the off-diagonal elements
We have found that in cases where we can estimate the likelihood reliably, as the weight changes so that the weighted matrix approaches the original matrix, the likelihood smoothly moves towards the 
likelihood of the original matrix. In cases with numerical problems, it moves smoothly for awhile and then jumps around. There are two approaches implemented. By default, with allow.extrapolation=FALSE,
we return a very high (poor) likelihood value if this jumping is observed. This serves to tell the search algorithm that the likelihood is very bad at this point (this is how bounds are commonly
implemented in programs). The problem with this approach is that the MLE may truly be in this region, it's just that numerically dealing with the matrices is hard for R (biology is not so convenient as to only
occur at parameter values for which the math is easy). The other approach,
allow.extrapolation=TRUE, uses this weighting and when the likelihood change stops being smooth, extrapolates from the points with the least amount of transformation from the original matrix. This might
give a more reasonable estimate of the likelihood if just a small change to the matrix is all that is needed, but if not can extrapolate to ridiculous likelihoods (sometimes ridiculously good ones, 
so the estimates will come from here).
} 

\value{
 \item{param.est}{an array that contains the maximum likelihood estimates.}
 \item{param.se}{The approximated standard error for the MLE estimates.}
 \item{loglik}{the maximum log-likelihood.}
 \item{AICc}{Akaike information criterion with sample size correction.}
 \item{Akaike.weight}{The Akaike weight \eqn{\exp(- \Delta AICc/2)}.}
 \item{wgt.param.est.se}{The weighted parameters and standard errors calculated through four models.}
 \item{convergence}{An integer code. 0 indicates successful completion. see optim}
 }

\author{
Dwueng-Chwuan Jhwueng <djhwueng@umail.iu.edu>, Brian O'Meara <bomeara@utk.edu>
}

\references{
Jhwueng D.C. and O'Meara B.C. 2014. \emph{Studying trait evolution in hybrid species on phylogenetic networks}. Submitted. 

Burnham, K.P., and D.R. Anderson. 2004. \emph{Model selection and inference: a practical information-theoretic approach}. Sec. Ed. Springer, New York.

Gilbert, P.  2011. numDeriv: Accurate Numerical Derivatives. R package version 
2010.11-1. http://CRAN.R-project.org/package=numDeriv

Nash, J.C. and R. Varadhan (2011). Unifying Optimization Algorithms to Aid Software System Users: optimx for R. Journal of Statistical Software, 43(9), 1-14. URL http://www.jstatsoft.org/v43/i09/.
}
\examples{

## A 15 taxon network, with three hybrids
\dontrun{network <- SimulateNetwork(ntax.nonhybrid = 12, ntax.hybrid=3, origin.type="individual")

#What does the flow data.frame look like:
print(network$flow)

params<-c(0.1, 0, 1, 0, 0)
names(params) <- c("sigma.sq", "mu", "bt", "vh", "SE")

data <- SimulateTipData(network$phy, network$flow, params)

results <- BMhyd(data, network$phy, network$flow, plot.se=FALSE)

#Now look at results
print(results)}

}