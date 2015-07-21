\name{AICc}
\alias{BMhyd}
\title{
Computing AICc
}
\description{
Function to compute the second order Akaike's information criterion 
}

\usage{
BMhyd(n, k,LogLik)
}
\arguments{
  \item{n}{
number of taxa
}
  \item{k}{
number of parameters in the model
}
  \item{LogLik}{
   negative log likelihood value  	
}


}
\details{
	'AICc' is a function to compute the AICc values and is valid to select among different models. 
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
Dwueng-Chwuan Jhwueng <dcjhwueng@fcu.edu.tw>
}

\references{
Burnham, K.P., and D.R. Anderson. 2004. \emph{Model selection and inference: a practical information-theoretic approach}. Sec. Ed. Springer, New York.
}
\examples{
## A 5 taxon phylogenetic network in eNewick format. It has four normal species and
n<-5
k<-3
LogLik<--2
AICc(n,k,LogLik)
}