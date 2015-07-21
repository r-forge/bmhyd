\name{AkaikeWeight}
\alias{BMhyd}
\title{
Akaike weight 
}
\description{
Calculate Akaike weight 
}

\usage{
AkaikeWeight(Delta.AICc.Array)
}
\arguments{
  \item{Delta.AICc.Array}{
The delta AICc value which is defined as the difference between the AICc value and the minimum AICc value among a candidate model. 
 }
}

\details{	
For n models of interest, the Akaikei weight for the i th model is defined as \eqn{w = exp(-0.5* dAICc_i)/ sum_i dAICc_i  } where dAICc_i is the delta AICc value for the ith model. 
} 

\value{
 \item{w}{Akaike weight}
 }

\author{
Brian O'Meara <omeara.brian@gmail.com>
}

\references{ 
	Burham, K.P., Anderson, D.R. (2002) Model selection and multimodel inference: a practical information-theoretic approach. Second edition. Springer. New York.
	}
\examples{
#simulate 4 AICc values for 4 models, respecitvely
AICc_array<-rnorm(4, mean=10,sd=1)
#calculate the delta AICc
Delta.AICc.array<-AICc_array-min(AICc_array)
#calculate the Akaike weight
AkaikeWeight(Delta.AICc.Array)
}