\name{BMhyd}
\alias{BMhyd}
\title{
 	Comparative method for studying hyridization using Brownian motion for trait evolution
	}

\description{
	This function fits the Brownian motion model of continuous character to investigate hybrid species through the hybrid vigor \eqn{\beta}, and the variation at the burst of hybridization \eqn{v_H}. Measurement error SE is also considered as well as the  parameters including the over all mean \eqn{\mu} and the overall variance \eqn{\sigma^2} in the typical Brownian motion model. 
	}

\usage{
	BMhyd(data, phy, flow, opt.method="Nelder-Mead", models=c(1,2,3,4), verbose=TRUE, 
	get.se=TRUE, plot.se=TRUE, store.sims=FALSE, precision=2, auto.adjust=FALSE, 
	likelihood.precision=0.001, allow.extrapolation=FALSE)
	}

\arguments{
  \item{data}{continuous trait data containing species information in vector format}
  \item{phy}{a tree in phylo class}
  \item{flow}{a struture of gene flow} 
  \item{opt.method}{the method for used for optimization. The default is Nelder-Mead}
  \item{models}{the model used for analysis (see details)}
  \item{verbose}{a TRUE/FALSE argument to start optimization}
  \item{get.se}{a TRUE/FALSE argument estimation for doing simulation to estimate parameter uncertainty (see details)}
  \item{plot.se}{a TRUE/FALSE argument for output the uncertainty plot for the model (see details)}
  \item{store.sims}{a TRUE/FALSE argument to record the the parameter estimates and relevant values.}
  \item{precision}{a numeric value to present the cutoff at which the user thinks the estimates become unreliable due to ill conditioned matrix.}
  \item{auto.adjust}{a TRUE/FALSE argument to adjust the the phylogeny }
  \item{likelihood.precision}{a numerical value used for verifying the convergent of the estimation.}
  \item{allow.extrapolation}{a TRUE/FALSE argument. If TRUE, the VCV matrix was ill-conditioned, so used splines to estimate its likelihood.}
	}

\details{Function BMhyd fits likelihood models for continuous characters. It incorporates phylogenetic tree, structure of gene flow 	and comparative data. The full likelihood model includes several parameters: the ancenstral state \eqn{\mu}, the overall variation 	\eqn{\sigma^2}, the hybrid vigor \eqn{\beta}, the hybrid burst variation at formation \eqn{v_H} and the measurement error SE. The 	struture of gene flow is a five column table where the first and the second column contain the donor and recipient information. 	The thrid column is the information about the heritibaility factor m which is a fraction of the recipient trait that comes from 	the source. The fourth column, recipient time, records time from the root of the recipient that counting forward from the root 		when the gene flow happened from the donor. The fifth column, recipient time, records time from the root of the recipient that 		counting forward from the root when the gene flow happened from the donor. For detail modeling, see the manuscript in Jhwueng and 
	O'Meara (2015). The function allows some fixed values of parameters and treat others as free parameters: model 1 fixes \eqn{\beta} 
	at 1 but 	allow \eqn{v_H} to vary; model 2 allows \eqn{\beta} to vary but fixes \eqn{v_H} at 0; model 3 fixes \eqn{\beta} at 1 
	and \eqn{v_H} 	at 0; and model 4 allows both to vary.  BMhyd fits the model through maximum likelihood technique, it returns 
	MLEs. When setting 	the get.se to TRUE, the BMhyd will do simulation to estimate parameter uncertainty and the confidence 
	intervals by adaptive 		confidence interval sampling method(see manuscript in Jhwueng and O'Meara). It setting plot.se to be 
	TRUE, BMhyd will be saved the 	CI in PDF file once the simulations are finished.  Model averaged parameter estimates are 
	calculated by the Akaike weight. As an 	input, \pkg{BMhyd} requires a phylogenetic tree of the phylo class, a structure of gene 
	flow and a comparative data. Currently the 	method is develop for univariate analysis where the comparative data includes a single 
	trait for analyses. 
 	} 

\value{
	A summaried table including the type of model, the corresponding number of parameter, the parameter estimates, the likelihood values, the upper bound and lower bound of 	the parameters, and the Akaike weights for model averaging. 
 	}

\author{
	Brian O'Meara, Dwueng-Chwuan Jhwueng.
	}

\references{
	Jhwueng D.C. and O'Meara B.C. 2015. \emph{Studying trait evolution in hybrid species on phylogenetic networks}. Submitted. 

	Burnham, K.P., and D.R. Anderson. 2004. \emph{Model selection and inference: a practical information-theoretic approach}. Sec. Ed. 	Springer, New York.
	}
	
\examples{
	#set up the number of non hybrid
	ntax.nonhybrid<-2
	#set up the number of non hybrid
	ntax.hybrid<-1
	#simulate a network
	network<-SimulateNetwork(ntax.nonhybrid=ntax.nonhybrid, ntax.hybrid=ntax.hybrid, 
	flow.proportion=0.5, origin.type='clade', birth = 1, death = 0.5, sample.f = 0.5,
	tree.height = 1, allow.ghost=FALSE)
    #simulate the tips data
	tips<-rnorm(ntax.nonhybrid+ntax.hybrid)
	names(tips)<-paste("t",(1:(ntax.nonhybrid+ntax.hybrid)),sep="")
	#run the analysis uses model 3
	\donttest{
	BMhyd(tips,network$phy,network$flow, opt.method="Nelder-Mead", models=3, verbose=TRUE,
	get.se=FALSE, plot.se=FALSE, store.sims=FALSE, precision=2, auto.adjust=FALSE, 
	likelihood.precision=0.001, allow.extrapolation=FALSE) 
	}
	}