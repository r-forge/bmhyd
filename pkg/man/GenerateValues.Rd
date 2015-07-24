\name{GenerateValues}
\alias{BMhyd}
\title{
	Generate parameter values 
	}

\description{
	Generate valid parameter values 
	}

\usage{
	GenerateValues(par, lower, upper, max.tries=100, expand.prob=0, examined.max, examined.min) 
	}

\arguments{
  \item{par}{parameters of interest}
  \item{lower}{the lower bound for the parameters}
  \item{upper}{the upper bound for the parameters} 
  \item{max.tries}{maximum number of attempts tp generate the parameter value. The default is set to 100.}
  \item{expand.prob}{a probability value to expand. Currently not used.}
  \item{examined.max}{the allowed maximum for the parameters values.}
  \item{examined.min}{the allowed minimum for the parameters values}
  }
 
\details{
	It generate new parameter values using the uniform distribution on the interval of 
	\eqn{ (L, U )}  where \eqn{ L=\max (lower, 0.9*examined.min )} and  \eqn{ U =\min (upper, 1.1*examined.max )}. 
	}
	
\value{
	The generated parameter values.
 	}

\author{
	Brian O'Meara <omeara.brian@gmail.com>, Dwueng-Chwuan Jhwueng <tony.jhwueng@gmail.com>
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2014. \emph{Studying trait evolution in hybrid species on phylogenetic networks}. Submitted. 
	}

\examples{
	mu<-1;sigma.sq<-0.1
	par<-c(mu,sigma.sq)
	lower=c(-2,0)
	upper<-c(2,1)
	examined.max<-c(1,1)
	examined.min<-c(0.001,0.001)
	GenerateValues(par, lower, upper, max.tries=100, expand.prob=0, examined.max, examined.min)
	}