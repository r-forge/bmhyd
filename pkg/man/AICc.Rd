\name{AICc}
\alias{BMhyd}
\title{
	AICc
	}

\description{
	Calculates the second order Akaike's information criterion 
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
	'AICc' is a function to compute the AICc values and is valid to select among different models. \ean{ AICc = 2*n*k/(n-k-1) -2 log L} where \eqn{L} is the maximum likelihood for the model.  
	} 

\value{
	The AICc values.
 	}

\author{
	Brian O'Meara <omeara.brian@gmail.com>, Dwueng-Chwuan Jhwueng <dcjhwueng@fcu.edu.tw>
	}

\references{
	Burnham, K.P., and D.R. Anderson. 2004. \emph{Model selection and inference: a practical information-theoretic approach}. Sec. Ed. Springer, New York.
	}
\examples{
	n<-5
	k<-3
	LogLik<- -2
	AICc(n,k,LogLik)
	}