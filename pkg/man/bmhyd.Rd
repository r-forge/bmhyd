\name{BMhyd}
\alias{BMhyd}
\title{
Comparative method for hybridzation
}
\description{
Fit Brownian motion model of continuous character to investigate hybrid species through the  hybrid vigor, 
heritability factor, and hybrid variation. 
}

\usage{
BMhyd(data, network)
}
\arguments{
  \item{data}{
continuous trait data containing species information in vector format (see Details).
}
  \item{network}{
a rooted phylogenetic network input in eNewick format (see Details).
}
}
\details{Function BMhyd fits likelihood models for continuous characters through incorporating phylogenetic network and comparative data. 
The full likelihood model includes five parameters: the ancenstral state \eqn{\mu}, the overall variation \eqn{\sigma^2}, 
the hybrid vigor \eqn{\beta}, proportion of heritability \eqn{0\leq h \leq 1}, and hybrid burst variation at formation \eqn{v_H}. 
Current version fixes heritability parameter, \eqn{h}, to 0.5. We use four submodels that allows some fixed values of parameters 
and treat others as free parameters: (1) \eqn{\beta=1, v_H=0}; (2) \eqn{\beta} is free parameter, \eqn{v_H=0};   
(3) \eqn{\beta=1, v_H} is free paramter (4)\eqn{\beta} and \eqn{v_H} are free parameters. 

BMhyd fits the model through maximum likelihood technique. It returns MLEs and their approximated standard errors for above four models. 
Model averaged parameter estimates and model averaged parameter standard errors for the four models are reported. 

As input, \pkg{BMhyd} requires (a) a network N with known topology and branch length and (b) a set of comparative data. 
The network must be entered in \emph{eNewick format} (Cardona et al. 2008). 
The comparative data format is set to the \eqn{n} by \eqn{1} vector that contains the continuous trait of \eqn{n} species of interest.  

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
Dwueng-Chwuan Jhwueng <djhwueng@umail.iu.edu>
}

\references{
Jhwueng D.C. and O'Meara B.C. 2014. \emph{Studying trait evolution in hybrid species on phylogenetic networks}. Submitted. 

Cardona, G., F. Rossello, and G. Valiente. 2008. \emph{Extended Newick: it is time for a standard representation of phylogenetic networks}. BMC Bioinformatics 9:532.

Burnham, K.P., and D.R. Anderson. 2004. \emph{Model selection and inference: a practical information-theoretic approach}. Sec. Ed. Springer, New York.

Gilbert, P.  2011. numDeriv: Accurate Numerical Derivatives. R package version 
2010.11-1. http://CRAN.R-project.org/package=numDeriv

Nash, J.C. and R. Varadhan (2011). Unifying Optimization Algorithms to Aid Software System Users: optimx for R. Journal of Statistical Software, 43(9), 1-14. URL http://www.jstatsoft.org/v43/i09/.
}
\examples{
## A 5 taxon phylogenetic network in eNewick format. It has four normal species and
#  one hybrid (species 2).
network<-c("(((1:0.5,(2:0.5)7:0)6:0.1,(7:0,3:0.5)8:0.1)10:0.4,(4:0.3,5:0.3)9:0.7)11:0;")
## Comparative data
data<-matrix(c( 1.2, 3.5, 1.6, 2.7, 4.3  ),ncol=1)
BMhyd(data,network)

}