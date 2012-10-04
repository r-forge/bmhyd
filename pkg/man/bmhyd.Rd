\name{bmhyd}
\alias{bmhyd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
comparative method for hybridzation
}
\description{
Fit Brownian motion model of conitnuous character to investigate hybrid species through the  hybrid vigor \eqn{\beta}, heritability factor \eqn{h}, and hybrid variation \eqn{v_H}. 
}

\usage{
bmhyd(data, network)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
a data vector of size n containing species information (see Details).
}
  \item{network}{
a rooted phylogenetic network (see Details).
}
}
\details{Function bmhyd fits likelihood models for continuous characters through incorporating phylogenetic network and comparative data. The full likelihood model includes five parameters: two of them are the ancenstral state \eqn{\mu}, and the overall variation \eqn{\sigma^2};  where three new paramters are hybrid vigor \eqn{\beta \geq 0}, proportion of heritability \eqn{ 0\leq h \leq 1}, and hybrid variation \eqn{v_H \geq 0}. Current version set \eqn{h} to 0.5 as a fixed parameter. In addition, several possible models are considered where models  allow some fixed values of parameters and treat others as free parameters: (1) \eqn{\beta=1, v_H=0} (2) \eqn{v_H=0, \beta} as free parameter  (3) \eqn{\beta=1, v_H} as free paramter (4) both \eqn{\beta} and \eqn{v_H} are free parameters. bmhyd analyzes data through maximum likelihood technique and returns MLEs estimates and their approximated standard errors for the four models. Model averaged parameters and model averaged standard errors will also be reported. The R package \pkg{optimx} provides a useful and open-source library for multinational optimization. As input, \pkg{bmhyd} requires a network \eqn{\aleph} with known topology and branch length and a comparative data. The network must be entered in \emph{eNewick format} (Cardona et al. 2008) . The comparative data is the \eqn{n} by \eqn{1} vector that contains the continuous trait of \eqn{n} species of interest.  

Mathematical Model:

Let \eqn{\bold{Y}=(y_1,y_2, ... ,y_n )} be the trait values for \eqn{n} species some of which are, possibly old,hybrids. According to the assumption of Brownian motion for trait evolution, \eqn{\bold{Y}} will follow the multivariate normal distribution \eqn{\bold{Y} \sim MVN(\bold{\mu},\bold{V}_R)} where \eqn{\bold{\mu}=(\mu_1,\mu_2,...,\mu_n )} and \eqn{\mu_i=\mu} or \eqn{\beta \mu} is the mean for usual species orthe mean for hybrid species; \eqn{\bold{V}_R = \bold{V}_R (\beta,h,\sigma^2,v_H|\bold{\aleph})} is the variance-covariance matrix, and \eqn{\bold{\aleph}} is a known network topology given with the branch lengths. The negative log likelihood function is 
\deqn{-\log L(\bold{\mu},\beta,h,\sigma^2,v_H |\bold{Y},\bold{\aleph}) = \frac{n}{2}\log 2\pi +\frac{1}{2} \log |\bold{V}d_R| 
+ \frac{1}{2} (\bold{Y}-\bold{\mu})^t \bold{V}_R^{-1}(\bold{Y}-\bold{\mu})}
where   \eqn{|\bold{V}_R|}  is the determinant of \eqn{\bold{V}_R} and \eqn{(\bold{Y}-\bold{\mu})^t} is the transpose of \eqn{(\bold{Y}-\bold{\mu})}
For estimation of standard error of the parameter estimate, the Hessian matrix is used. The variance-covaraince matrix of the estimate of \eqn{\mu, \beta, \sigma^2, v[H]} is obtained by evaluating the inverse of the Hessian matrix. The standard error for each paramter estimate is equal to the squares root of the corresponding diagonal element of the Hessian matrix. 

The model average parameter \eqn{\hat{\bar{\theta}}} and model average standard error \eqn{se(\hat{\bar{\theta}})} is calculated by (Burham and Anderson 2004) in the following. 

 \eqn{\hat{\bar{\theta}}}=\sum_{i=1}^4 w_i\hat{\theta}_i where \eqn{w_i, i=1,1,\cdots,4} is the Akakei weight which is just \eqn{\exp( -0.5 * \Delta AICc_i ). 

and

 \eqn{se(\hat{\bar{\theta}})= \sum_{i=1}^4\sqrt{\hat{var}(\hat{\theta}_i|g_i)+(\hat{\theta}_i-\hat{\bar{\theta}})^2}} where \hat{var(\hat{\theta_i|g_i})} is the variation of MLE \eqn{\hat{\theta}_i}  in
the ith model and is calculated by taking the diagonal element of the inverse of the Hessian matrix through R package numDeriv (Gilbert 2011).

}
\value{
 \item{loglik}{the maximum log-likelihood.}
 \item{AICc}{Akaike information criterion with sample size correction.}
 \item{param.est.}{an array that contains the maximum likelihood estimates.}
 \item{param.se}{The approximated standard error for the MLE estimates.}
 \item{Akaike.weight}{The Akaike weight calculated by \eqn{\exp(- \Delta AICc/2)}.}
 \item{wgt.param.est}{The weighted parameter calculated through four models.}
 }

\author{
Dwueng-Chwuan Jhwueng <djhwueng@umail.iu.edu>
}

\references{
Jhwueng D.C. and O'Meara B.C. 2012. \emph{Studying trait evolution in hybrid species on phylogenetic networks}. submitted. 

Cardona, G., F. Rossello, and G. Valiente. 2008. \emph{Extended Newick: it is time for a standard representation of phylogenetic networks}. BMC Bioinformatics 9:532.

Burnham, K.P., and D.R. Anderson. 2004. \emph{Model selection and inference: a practical information—theoretic approach}. Sec. Ed. Springer, New York.

Gilbert, P.  2011. numDeriv: Accurate Numerical Derivatives. R package version 
2010.11-1. http://CRAN.R-project.org/package=numDeriv

Nash, J.C. and R. Varadhan (2011). Unifying Optimization Algorithms to Aid Software System Users: optimx for R. Journal of Statistical Software, 43(9), 1-14. URL http://www.jstatsoft.org/v43/i09/.
}

\examples{
## A 9 taxa phylogenetic network in eNewick format. 7 normal species and 2 hybrid. 
network<-"(((1:0.35,2:0.35)10:0.5,(3:0.7,(4:0.45,(5:0.1,(6:0.1)13:0)12:0.35) 
11:0.25)16:0.15)18:0.15,((13:0,7:0.1)14:0.55,(8:0.2,9:0.2)15:0.45)17:0.35)19:0;"
## Comparative data
data<-matrix(c( 16.00,12.00,16.00,23.00,18.00,18.00,25.00,20.00,14.00  ),ncol=1)
bmhyd(data,network)
}
