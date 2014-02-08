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
Jhwueng D.C. and O'Meara B.C. 2013. \emph{Studying trait evolution in hybrid species on phylogenetic networks}. submitted. 

Cardona, G., F. Rossello, and G. Valiente. 2008. \emph{Extended Newick: it is time for a standard representation of phylogenetic networks}. BMC Bioinformatics 9:532.

Burnham, K.P., and D.R. Anderson. 2004. \emph{Model selection and inference: a practical information-theoretic approach}. Sec. Ed. Springer, New York.

Gilbert, P.  2011. numDeriv: Accurate Numerical Derivatives. R package version 
2010.11-1. http://CRAN.R-project.org/package=numDeriv

Nash, J.C. and R. Varadhan (2011). Unifying Optimization Algorithms to Aid Software System Users: optimx for R. Journal of Statistical Software, 43(9), 1-14. URL http://www.jstatsoft.org/v43/i09/.
}
\examples{
## A 9 taxa phylogenetic network in eNewick format. 7 normal species and 2 hybrid.
network<-c("(((1:0.35,2:0.35)10:0.5,(3:0.7,(4:0.45,(5:0.1,(6:0.1)13:0)12:0.35)11:0.25)16:0.15)18:0.15,((13:0,7:0.1)14:0.55,(8:0.2,9:0.2)15:0.45)17:0.35)19:0;")
## Comparative data
data<-matrix(c( 16.00,12.00,16.00,23.00,18.00,18.00,25.00,20.00,14.00  ),ncol=1)
BMhyd(data,network)


## A 22 taxa phylogenetic network in eNewick format. Nodes a,b,c represent the extinct taxa. Not run due to speed issues
#data(cichlid)
#network<-c("((((((1:0.495,((((2:0.215,(5:0.215)27:0)26:0.1,(14:0.315)34:0)33:0.055,(a:0.1,(3:0.1)21:0)20:0.27)37:0.06,((21:0,4:0.1)22:0.26,(27:0,(b:0.17,(6:0.17)24:0)23:0.045)28:0.145)36:0.07)38:0.065)39:0.11,7:0.605)41:0.045,((8:0.575,9:0.575)40:0.06,10:0.635)42:0.015)43:0.125,((((24:0,11:0.17)25:0.075,12:0.245)29:0.475,(13:0.68,(34:0,(15:0.275,(16:0.275)31:0)30:0.04)35:0.365)44:0.04)45:0.02,(31:0,c:0.275)32:0.465)46:0.035)47:0.11,(17:0.835,18:0.835)48:0.05)49:0.115,19:1)50:0;")
## Comparative data. NA represent the data for extinct species.
#data<-matrix(c( 13.5,12.4,NA,7.0,5.8,5.5,NA,16.0,15.0,25.0,6.1,6.5,5.5,11.0,7.8,15.0,4.3,4.0,NA,8.6,10.3,8.5),ncol=1)
#BMhyd(data,network)

}