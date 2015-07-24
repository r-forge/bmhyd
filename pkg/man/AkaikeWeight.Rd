\name{AkaikeWeight}
\alias{BMhyd}
\title{
	Akaike weight 
	}

\description{
	Calculate Akaike weight for the models 
	}

\usage{
	AkaikeWeight(Delta.AICc.Array)
	}

\arguments{
  \item{Delta.AICc.Array}{
	Delta.AICc.Array is defined as the difference between the AICc value and the minimum AICc value among candidate models. 
 	}
	}

\details{	
	For \eqn{n} models of interest, the Akaikei weight for the \eqn{i} th model is defined as \eqn{w = exp(-0.5* dAICc_i)/ sum_i dAICc_i  } where dAICc_i is the delta AICc value for the \eqn{i} th model. 
	} 

\value{
 	Akaike weights that appear in vector format.
 	}

\author{
	Dwueng-Chwuan Jhwueng <dcjhwueng@fcu.edu.tw>
	}

\references{ 
	Burham, K.P., Anderson, D.R. (2002) Model selection and multimodel inference: a practical information-theoretic approach. Second edition. Springer. New York.
	}
\examples{
	#simulate 4 AICc values for 4 models, respecitvely
	AICc_Array<-rnorm(4, mean=10,sd=1)
	#calculate the delta AICc
	Delta.AICc.Array<-AICc_Array-min(AICc_Array)
	#calculate the Akaike weight
	AkaikeWeight(Delta.AICc.Array)
	}