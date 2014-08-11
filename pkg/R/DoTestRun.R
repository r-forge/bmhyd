library(BMhyd)

non.hybrid.vector <- c(30, 100)
hybrid.vector <- c(1, 5, 10)
flow.vector <-c(0, .1, .5)
origins.vector <- c("clade", "individual")
sigma.vector <- c(0.01) #maybe?
mu.vector <- c(1)
bt.vector <- c(0.1, 1, 10)
vh.vector <- c(0, 1, 10)
SE.vector <- c(0)
tree.height.vector <- c(50)
nreps <- 50
possible.sets <- expand.grid(ntax.nonhybrid=non.hybrid.vector, ntax.hybrid=hybrid.vector, flow.proportion=flow.vector, origin.type=origins.vector, sigma.sq=sigma.vector, mu=mu.vector, bt=bt.vector, vh=vh.vector, SE=SE.vector, tree.height=tree.height.vector)


#names(free.parameters) <- c("sigma.sq", "mu", "bt", "vh", "SE")
#GetMeansModified <- function(x, phy, flow, actual.params) {
#params must be named vector
#SimulateTipData <- function(phy, flow, params) {

#SimulateNetwork <- function(ntax.nonhybrid=100, ntax.hybrid=10, flow.proportion=0.5, origin.type=c("clade", "individual"), birth = 1, death = 0.5, sample.f = 0.5, tree.height = 1, allow.ghost=FALSE) {

RunFromSet<-function(x, id) {
	params<-rep(NA,5)
	names(params) <- c("sigma.sq", "mu", "bt", "vh", "SE")
	for (i in sequence(length(params))) {
		params[i]<-x[which(names(x)==names(params)[i])]	
	}
	params<-unlist(params)
	network<-SimulateNetwork(ntax.nonhybrid=as.numeric(x['ntax.nonhybrid']), ntax.hybrid=as.numeric(x['ntax.hybrid']), flow.proportion=as.numeric(x['flow.proportion']), origin.type=as.character(x['origin.type']), birth = 1, death = 0.5, sample.f = 0.5, tree.height = as.numeric(x['tree.height']), allow.ghost=FALSE)
	tips<-SimulateTipData(network$phy, network$flow, params)
	results <- BMhyd(tips, network$phy, network$flow, get.se=FALSE, precision=5, verbose=TRUE)
	save(list=ls(), file=paste(id, ".RData"))
}

#RunFromSet(possible.sets[1,], id="1")
library(foreach)
library(doMC)
options(cores=12)
foreach(i=sequence(dim(possible.sets)[1])) %dopar% RunFromSet(possible.sets[i,], id=as.character(i))
