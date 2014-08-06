
AICc<-function(n,k,LogLik){
  return(2*n*k/(n-k-1)+2*LogLik)
}

AkaikeWeight<-function(Delta.AICc.Array){
  return(exp(-Delta.AICc.Array/2) /sum(exp(-Delta.AICc.Array/2) ))
}

se.function<-function(cov.matrix,var.name){
  name.Index<-which(rownames(cov.matrix)==var.name)
  if( length(name.Index)==1){
    return( cov.matrix[name.Index,name.Index])    
  }else{return(0)}
}

var.model.Index.function<-
  function(cov.matrix,var.name){
    name.Index<-which(rownames(cov.matrix)==var.name)
    if( length(name.Index)==1){
      return( cov.matrix[name.Index,name.Index])    
    }else{return(0)}
  }

weight.para.value<-
  function(para.vect,weight){
    return(sum(para.vect*weight))
  }


####################################################################
###################### MAIN PROGRAM ################################
####################################################################

#data is vector with names() matching taxon names
#phy is an ape phylo object
#flow is a data.frame with four columns
#donor = the taxon that is the gene flow source
#recipient = the taxon that is the gene flow recipient
#m = the fraction of the recipient trait that comes from the source. In the case of an equal hybridization between the recipient's sister on the tree and the donor, this is 0.5. In other cases where only, say, 10% of the recipient's quantitative trait
#	loci come from the donor, it would be 0.1
#time.from.root.donor = the time, counting forward FROM THE ROOT, when the gene flow happened from the donor. It may not be the same as time.from.root.recipient, as it may have spent time in a now extinct ghost lineage first (though time.from.root.donor <= time.from.root.recipient). It's treated as a one time event, which makes sense in the case of a single allopolyploid speciation event, probably less so in the case
#	of ongoing gene flow. Too bad. 
#time.from.root.recipient = the time, counting forward FROM THE ROOT, when the gene flow happened from the donor
#If the gene flow happened to or from a lineage with multiple descendant species, use one row for each pair. For example, if lineage (A,B) had 20% of their genes coming in from lineage (C,D,E) at 14.5 MY since the root (not back in time), you would have
#	a flow data.frame of
#donor	recipient	m	time.from.root.donor	time.from.root.recipient
#C		A			0.2	14.5					14.5
#D		A			0.2	14.5					14.5
#E		A			0.2	14.5					14.5
#C		B			0.2	14.5					14.5
#D		B			0.2	14.5					14.5
#E		B			0.2	14.5					14.5
#We may write a utility function for dealing with this case in the future.
#Note the use of all updates of V.modified based on V.original; we don't want to add v_h to A three different times, for example, for one migration event (so we replace the variance three times based on transformations of the original variance)
#Note that we do not assume an ultrametric tree
BMhyd <- function(data, phy, flow, opt.method="Nelder-Mead", models=c(1,2,3,4), verbose=TRUE, get.se=TRUE, plot.se=TRUE, store.sims=FALSE, precision=2) {
	if(min(flow$m)<0) {
		stop("Min value of flow is too low; should be between zero and one")	
	}
	if(max(flow$m)>1) {
		stop("Min value of flow is too high; should be between zero and one")	
	}
	results<-list()
	#hessians <- list()
	results.summary <-data.frame()
	phy.geiger.friendly <- phy #geiger can't handle branch lengths near zero. Let's lengthen them if needed
	if(min(phy.geiger.friendly$edge.length)<0.00001) {
		phy.geiger.friendly$edge.length <- phy.geiger.friendly$edge.length + 0.00001
	}
	phy <- AdjustForDet(phy)
	all.sims<-list()
	if(verbose) {
		print("Getting starting values from Geiger")	
	}

	starting.from.geiger<-fitContinuous(phy.geiger.friendly, data, model="BM", SE=data*NA)$opt
	starting.values <- c(starting.from.geiger$sigsq, starting.from.geiger$z0, 1,  0, starting.from.geiger$SE) #sigma.sq, mu, beta, vh, SE
	if(verbose) {
		print("Done getting starting values")
	}
	for (model.index in sequence(length(models))) {
		step.count <- 0
		if(verbose) {
			print(paste("Starting model", model.index, "of", length(models)))
		}
		free.parameters<-rep(TRUE, 5)
		names(free.parameters) <- c("sigma.sq", "mu", "bt", "vh", "SE")
		model <- models[model.index]
		if(model==1) {
			free.parameters[which(names(free.parameters)=="bt")]<-FALSE
		}
		if(model==2) {
			free.parameters[which(names(free.parameters)=="vh")]<-FALSE
		}
		if(model==3) {
			free.parameters[which(names(free.parameters)=="bt")]<-FALSE
			free.parameters[which(names(free.parameters)=="vh")]<-FALSE
		}
		best.run <- optim(par=starting.values[free.parameters], fn=CalculateLikelihood, method=opt.method, hessian=FALSE, data=data, phy=phy, flow=flow, actual.params=free.parameters[which(free.parameters)], precision=precision)
		if(verbose) {
			results.vector<-c(step.count, best.run$value, best.run$par)
			names(results.vector) <- c("step","negloglik", names(free.parameters[which(free.parameters)]))
			print(results.vector)
		}
		#this is to continue optimizing; we find that optim is too lenient about when it accepts convergence
		times.without.improvement <- 0
		while(times.without.improvement<10) {
			times.without.improvement <- times.without.improvement+1
			new.run <- best.run
			if(times.without.improvement%%4==0) {
				new.run <- optim(par=best.run$par, fn=CalculateLikelihood, method=opt.method, hessian=FALSE, data=data, phy=phy, flow=flow, actual.params=free.parameters[which(free.parameters)], precision=precision)
			} else {
				new.run <- optim(par=GenerateValues(best.run$par, lower=c(0, -Inf, 0, 0, 0)[which(free.parameters)], upper=rep(Inf, sum(free.parameters)), examined.max=10*best.run$par, examined.min=0.1*best.run$par), fn=CalculateLikelihood, method=opt.method, hessian=FALSE, data=data, phy=phy, flow=flow, actual.params=free.parameters[which(free.parameters)], precision=precision)					
			}
			#print("new.run best.run")
			#print(c(new.run$value, best.run$value))
			if(new.run$value<best.run$value) {
				best.run <- new.run
				times.without.improvement <- 0
				if(verbose) {
					print("New improvement found")	
				}
			}
			if(verbose) {
				step.count <- step.count+1
				results.vector<-c(step.count, times.without.improvement, best.run$value, best.run$par)
				names(results.vector) <- c("step", "steps.without.improvement","negloglik", names(free.parameters[which(free.parameters)]))
				print(results.vector)
			}
		}
		results[[model.index]] <- best.run
		#try(hessians[[model.index]] <- hessian(func=CalculateLikelihood, x=new.run$par, data=data, phy=phy, flow=flow, actual.params=free.parameters[which(free.parameters)]))
		results.vector.full <- c(NA, NA, 1, 0, 0)
		names(results.vector.full) <- names(free.parameters)
		names(best.run$par) <- names(free.parameters[which(free.parameters)])
		for (i in sequence(length(best.run$par))) {
			results.vector.full[which(names(results.vector.full)==names(best.run$par)[i])] <- best.run$par[i]
		}
		#print(hessians[[model.index]])
		#try(print(solve(hessians[[model.index]])))
		ci.vector<-rep(NA,10)
		for(parameter in sequence(5)) {
			names(ci.vector)[1+2*(parameter-1)] <- paste(names(free.parameters)[parameter],"lower", sep=".")
			names(ci.vector)[2+2*(parameter-1)] <- paste(names(free.parameters)[parameter],"upper", sep=".")
		}

		if(get.se) {
			if(verbose) {
				print("Now doing simulation to estimate SE")	
			}
			interval.results <- AdaptiveConfidenceIntervalSampling(best.run$par, fn=CalculateLikelihood, lower=c(0, -Inf, 0, 0, 0)[which(free.parameters)], data=data, phy=phy, flow=flow, actual.params=free.parameters[which(free.parameters)])
			interval.results.in <- interval.results[which(interval.results[,1]-min(interval.results[,1])<=2),]
			interval.results.out <- interval.results[which(interval.results[,1]-min(interval.results[,1])>2),]
			if(plot.se) {
				pdf(file=paste("Model",models[model.index], "_SE_plot.pdf", sep=""), height=5, width=5*sum(free.parameters))
				par(mfcol=c(1, sum(free.parameters)))
				for(parameter in sequence(sum(free.parameters))) {
					plot(x=interval.results[,parameter+1], y=interval.results[,1], type="n", xlab=names(free.parameters[which(free.parameters)])[parameter], ylab="NegLnL", bty="n", ylim=c(min(interval.results[,1]), min(interval.results[,1])+10))
					points(x=interval.results.in[,parameter+1], y=interval.results.in[,1], pch=16, col="black")
					points(x=interval.results.out[,parameter+1], y=interval.results.out[,1], pch=16, col="gray")
					points(x= best.run$par[parameter], y= best.run$value, pch=1, col="red", cex=1.5)
				}
				dev.off()
				if(verbose) {
					print(paste("SE plot has been saved in Model",models[model.index], "_SE_plot.pdf in ", getwd(), sep=""))
				}
			}
			if(store.sims) {
				all.sims<-append(all.sims, interval.results)
			}
			free.index=0
			for(parameter in sequence(5)) {
				
				if(free.parameters[parameter]) { #is estimated
					free.index <- free.index + 1
					ci.vector[1+2*(parameter-1)] <- min(interval.results.in[,free.index+1])
					ci.vector[2+2*(parameter-1)] <- max(interval.results.in[,free.index+1])
				} else {
					ci.vector[1+2*(parameter-1)] <- results.vector.full[parameter]
					ci.vector[2+2*(parameter-1)] <- results.vector.full[parameter]	
				}
			}
		}
		local.df <- data.frame(matrix(c(model.index, results.vector.full, AICc(Ntip(phy),k=length(free.parameters[which(free.parameters)]), best.run$value), best.run$value, length(free.parameters[which(free.parameters)]), ci.vector), nrow=1))
		colnames(local.df) <- c("Model", names(results.vector.full), "AICc", "NegLogL", "K", names(ci.vector))
		print(local.df)
		results.summary <- rbind(results.summary, local.df)

	}
	results.summary <- cbind(results.summary, deltaAICc=results.summary$AICc-min(results.summary$AICc))
	results.summary<-cbind(results.summary, AkaikeWeight = AkaikeWeight(results.summary$deltaAICc))
	if(store.sims) {
		return(list(results=results.summary, sims=all.sims))	
	}
	return(results.summary)
}

DetPass <- function(phy) {
	det.pass <- TRUE
	vcv.result <- vcv.phylo(phy)
	det.tries <- c(det(vcv.result), det(1000*vcv.result), det(0.0001*vcv.result))
	if(min(det.tries)<0) {
		det.pass <- FALSE
	}
	if(sum(is.finite(det.tries))!=length(det.tries)) {
		det.pass <- FALSE
	}
	return(det.pass)
}

AdjustForDet <- function(phy, max.attempts=100) {
	attempts<-0
	if(!DetPass(phy)) {
		warning("Determininant of the phylogeny was difficulty to calculate, so the phylogeny needed to be adjusted. Your results may be approximate as a result")
		while(!DetPass(phy) && attempts<=max.attempts) {
			attempts <- attempts+1
			phy$edge.length <- phy$edge.length+0.00001*attempts
		}
	}
	if(attempts>max.attempts) {
		stop("Phylogeny could not be adjusted adequately")
	}
	return(phy)
}

GetVModified <- function(x, phy, flow, actual.params) {
	bt <- 1
	vh <- 0
	sigma.sq <- x[1]
	mu <- x[2]
	SE <- x[length(x)]
	bt.location <- which(names(actual.params)=="bt")
	if(length(bt.location)==1) {
		bt<-x[bt.location]	
	}
	vh.location <- which(names(actual.params)=="vh")
	if(length(vh.location)==1) {
		vh<-x[vh.location]	
	}
	times.original <-vcv.phylo(phy, model="Brownian") #is the initial one based on the tree alone, so just time
	V.original <- sigma.sq * times.original 
	V.modified <- V.original
	for (flow.index in sequence(dim(flow)[1])) {
		recipient.index <- which(rownames(V.modified)==flow$recipient[flow.index])
		if(length(recipient.index)!=1) {
			stop(paste("Tried to find ", flow$recipient[flow.index], " but instead found ", paste(rownames(V.modified)[recipient.index], sep=" ", collapse= " "), "; make sure the taxon names in the flow dataframe recipient match that of your tree", sep=""))
		}
		donor.index <- which(rownames(V.modified)==flow$donor[flow.index])
		if(length(donor.index)!=1) {
			stop(paste("Tried to find ", flow$donor[flow.index], " but instead found ", paste(rownames(V.modified)[donor.index], sep=" ", collapse= " "), "; make sure the taxon names in the flow dataframe donor match that of your tree", sep=""))
		}
		V.modified[recipient.index, donor.index] <- (1-flow$m[flow.index]) * V.original[recipient.index, donor.index] + (flow$m[flow.index]) * (flow$time.from.root.recipient[flow.index]) * sigma.sq #covariance is the weighted sum of the covariance from evolution along the tree plus evolution along the migration path
		V.modified[donor.index, recipient.index] <- V.modified[recipient.index, donor.index]
		#covariance managed, now to manage the variance
		V.modified[recipient.index, recipient.index] <- V.original[recipient.index, recipient.index] + vh
	}
	diag(V.modified) <- diag(V.modified)+SE
	return(V.modified)
}

GetMeansModified <- function(x, phy, flow, actual.params) {
	badval<-(0.5)*.Machine$double.xmax
	bt <- 1
	vh <- 0
	sigma.sq <- x[1]
	mu <- x[2]
	SE <- x[length(x)]
	bt.location <- which(names(actual.params)=="bt")
	if(length(bt.location)==1) {
		bt<-x[bt.location]	
	}
	vh.location <- which(names(actual.params)=="vh")
	if(length(vh.location)==1) {
		vh<-x[vh.location]	
	}
	times.original <-vcv.phylo(phy, model="Brownian") #is the initial one based on the tree alone, so just time
	V.original <- sigma.sq * times.original 

	means.original <- rep(mu, Ntip(phy))
	names(means.original) <- rownames(V.original)
	means.modified <- means.original

	means.original <- rep(mu, Ntip(phy))
	names(means.original) <- rownames(V.original)
	means.modified <- means.original
	for (flow.index in sequence(dim(flow)[1])) {
		recipient.index <- which(rownames(V.original)==flow$recipient[flow.index])
		if(length(recipient.index)!=1) {
			stop(paste("Tried to find ", flow$recipient[flow.index], " but instead found ", paste(rownames(V.original)[recipient.index], sep=" ", collapse= " "), "; make sure the taxon names in the flow dataframe recipient match that of your tree", sep=""))
		}
		donor.index <- which(rownames(V.original)==flow$donor[flow.index])
		if(length(donor.index)!=1) {
			stop(paste("Tried to find ", flow$donor[flow.index], " but instead found ", paste(rownames(V.original)[donor.index], sep=" ", collapse= " "), "; make sure the taxon names in the flow dataframe donor match that of your tree", sep=""))
		}
		means.modified[recipient.index] <- means.original[recipient.index] + log(bt)
	}
	return(means.modified)
}



#precision is the cutoff at which we think the estimates become unreliable due to ill conditioned matrix
CalculateLikelihood <- function(x, data, phy, flow, actual.params, precision=2, proportion.mix.with.diag=0) {
	badval<-(0.5)*.Machine$double.xmax
	bt <- 1
	vh <- 0
	sigma.sq <- x[1]
	mu <- x[2]
	SE <- x[length(x)]
	bt.location <- which(names(actual.params)=="bt")
	if(length(bt.location)==1) {
		bt<-x[bt.location]	
	}
	vh.location <- which(names(actual.params)=="vh")
	if(length(vh.location)==1) {
		vh<-x[vh.location]	
	}
	V.modified <- GetVModified(x, phy, flow, actual.params)
	means.modified <- GetMeansModified(x, phy, flow, actual.params)
	if(sigma.sq <0 || vh<0 || bt <= 0.0000001 || SE < 0) {
    	return(badval)
	}
	data <- data[match(names(means.modified), names(data))]
	if(length(data)!=length(means.modified)) {
		stop("Mismatch between names of taxa in data vector and on phy")	
	}
	#if(length(which(eigen(V.modified)$values<0))>0) {
	#	last.bad <- V.modified
	#	return(badval)	
	#}
	NegLogML <- (Ntip(phy)/2)*log(2*pi)+(1/2)*t(data-means.modified)%*%pseudoinverse(V.modified)%*%(data-means.modified) + (1/2)*log(abs(det(V.modified))) 
	
	if(min(V.modified)<0 || sigma.sq <0 || vh<0 || bt <= 0.0000001 || !is.finite(NegLogML) || SE<0) {
    	NegLogML<-badval 
	}
	matrix.condition <- kappa(V.modified, exact=TRUE)
	#print("condition")
	#print(kappa(V.modified, exact=TRUE))
	#print("log(condition)")
	#print(log(kappa(V.modified, exact=TRUE)))
	
	#pretty<-c(NegLogML, log(matrix.condition))
	#names(pretty) <- c("NegLogL", "log(matrix.condition")
	#print(pretty)
	#The ratio  of the largest to smallest singular value in the singular value decomposition of a matrix. The base- logarithm of  is an estimate of how many base- digits are lost in solving a linear system with that matrix. In other words, it 
	#estimates worst-case loss of precision. A system is said to be singular if the condition number is infinite, and ill-conditioned if it is too large, where "too large" means roughly  the precision of matrix entries.
	#if(rcond(V.modified) < .Machine$double.eps^2){
	if(log(matrix.condition) > precision) {
		proportions <- seq(from=1, to=0, length.out=101) 
		lnl.vector <- rep(NA, length(proportions))
		max.diff <- 0
		for(i in sequence(length(proportions))) {
			V.modified.by.proportions<-(1-proportions[i]) * V.modified + proportions[i] * diag(dim(V.modified)[1]) * diag(V.modified)
			local.lnl <- (Ntip(phy)/2)*log(2*pi)+(1/2)*t(data-means.modified)%*%pseudoinverse(V.modified.by.proportions)%*%(data-means.modified) + (1/2)*log(abs(det(V.modified.by.proportions))) 
			if(i>6) {
				very.local.lnl <- lnl.vector[(i-6):(i-1)]
				max.diff <- max(abs(very.local.lnl[-1] - very.local.lnl[-length(very.local.lnl)])) #looking locally for jumps in the likelihood
				current.diff <- abs(local.lnl - lnl.vector[i-1])
				if(current.diff > 2 * max.diff) {
					#print(paste("breaking after ", i))
					break() #the modified matrix is still poorly conditioned, so stop here	
				}	
			}
			lnl.vector[i] <- local.lnl
		}
		proportions<-proportions[which(!is.na(lnl.vector))]
		lnl.vector<-lnl.vector[which(!is.na(lnl.vector))]
		NegLogML <- predict(smooth.spline(proportions, lnl.vector), data.frame(proportions =0.000))$y
		#plot(c(0, proportions), c(NegLogML, lnl.vector), type="n")
		#points(proportions, lnl.vector, pch=20)
		#points(0, NegLogML, col="red")

		#print(paste("Did interpolation, got ", NegLogML))
		warning("VCV matrix was ill-conditioned, so used splines to estimate its likelihood")
	}
	#print("datadiff")
	#print(quantile(data-means.modified))
	#print("middle")
	#print((1/2)*t(data-means.modified)%*%pseudoinverse(V.modified)%*%(data-means.modified))
	#print("end")
	#print((1/2)*log(abs(det(V.modified))) )
	#print(x)
	#print(V.modified[1:10,1:10])
	#print(means.modified)
	if(NegLogML< (-1000)) {
		bad<-V.modified
	}
	return(NegLogML[1])
}

AdaptiveConfidenceIntervalSampling <- function(par, fn, lower=-Inf, upper=Inf, desired.delta = 2, n.points=5000, verbose=TRUE, ...) {
	starting<-fn(par, ...)
	if(length(lower) < length(par)) {
		lower<-rep(lower, length(par))
	}
	if(length(upper) < length(par)) {
		upper<-rep(upper, length(par))
	}
	min.multipliers <- rep(1, length(par))
	max.multipliers <- rep(1, length(par))
	results<-data.frame(data.frame(matrix(nrow=n.points+1, ncol=1+length(par))))
	results[1,]<-unname(c(starting, par))
	for (i in sequence(n.points)) {
		sim.points<-NA
		while(is.na(sim.points[1])) {
			sim.points<-GenerateValues(par, lower, upper, examined.max=max.multipliers*apply(results[which(results[,1]-min(results[,1], na.rm=TRUE)<=desired.delta),-1], 2, max, na.rm=TRUE), examined.min=min.multipliers*apply(results[which(results[,1]-min(results[,1], na.rm=TRUE)<=desired.delta),-1], 2, min, na.rm=TRUE))
		}
		results[i+1,] <- c(fn(sim.points, ...), sim.points)
		if (i%%20==0) {
			for (j in sequence(length(par))) {
				returned.range <- range(results[which((results[,1]-min(results[,1], na.rm=TRUE))<desired.delta), j+1], na.rm=TRUE)
				total.range <- range(results[,j+1], na.rm=TRUE)
				width.ratio <- diff(returned.range)/diff(total.range)
				if(is.na(width.ratio)) {
					width.ratio=1	
				}
				if(width.ratio > 0.5) { #we are not sampling widely enough
					min.multipliers[j] <- min.multipliers[j] * 0.9
					max.multipliers[j] <- max.multipliers[j] * 1.1 #expand the range
				} else {
					min.multipliers[j] <- 1
					max.multipliers[j] <- 1
				}
			}
		}
		if (verbose && i%%100==0) {
			print(paste(i, "of", n.points, "done"))	
		}
	}
	return(results)
}

GenerateValues <- function(par, lower, upper, max.tries=100, expand.prob=0, examined.max, examined.min) {
	pass=FALSE
	tries=0
	while(!pass && tries<=max.tries) {
		tries <- tries+1
		pass=TRUE
		new.vals <- rep(NA, length(par))
		for(i in sequence(length(par))) {
			examined.max[i]<-max(0.001, examined.max[i])
			new.vals[i]<-runif(1, max(lower[i], 0.9*examined.min[i]), min(upper[i], 1.1*examined.max[i]))
			if(new.vals[i]<lower[i]) {
				pass=FALSE
			}
			if(new.vals[i]>upper[i]) {
				pass=FALSE
			}
		}
	}
	if(tries>max.tries) {
		return(NA)
	}
	return(new.vals)
}

GetClade <- function(phy, clade.size) {
	nodes <- phy$edge[,1]
	subtrees <- lapply(nodes, extract.clade, phy=phy)
	counts <- sapply(subtrees, Ntip)
	matches<-subtrees[which(counts==clade.size)]
	if(length(matches)==0) {
		return(NA)	
	}
	lucky <- matches[sample.int(length(matches),1)][[1]]
	return(findMRCA(phy, tips=lucky$tip.label, type="node"))	
}

GetAncestor <- function(phy, node) {
	return(phy$edge[which(phy$edge[,2]==node),1])	
}


#allow.ghost allows ghost lineage: something that persists for awhile, hybridizes, goes extinct. Otherwise, hybridization events must between coeval edges with extant descendants
SimulateNetwork <- function(ntax.nonhybrid=100, ntax.hybrid=10, flow.proportion=0.5, origin.type=c("clade", "individual"), birth = 1, death = 0.5, sample.f = 0.5, tree.height = 1, allow.ghost=FALSE) {
	library(TreeSim)
	library(phytools)
	done = FALSE
	used.recipients <- c()
	available.recipient.ids <- sequence(ntax.nonhybrid + ntax.hybrid)
	flow <- data.frame()
	phy<-NA
	phy <-  sim.bd.taxa.age(n=ntax.nonhybrid+ntax.hybrid, numbsim=1, lambda=birth, mu=death, frac = sample.f, age=tree.height, mrca = TRUE)[[1]]
	if(origin.type=="clade" && ntax.hybrid==1) {
		warning("For ntax.hybrid = 1 and clade sampling, this will do individual sampling instead (which is equivalent in this case)")
		origin.type<-"individual"	
	}
	if(origin.type=="clade") {
		while(is.na(GetClade(phy, ntax.hybrid))) { #not all trees of a given size have a clade of a given size, so may need to resimulate it
			phy <-  sim.bd.taxa.age(n=ntax.nonhybrid+ntax.hybrid, numbsim=1, lambda=birth, mu=0.5, frac = sample.f, age=tree.height, mrca = TRUE)[[1]]
		}
	}
	while(!done) {
		donors <- c()
		recipients <- c()
		recipient.ids <- c()
		focal.node <- c()
		if (origin.type=="clade") {
			focal.node <- GetClade(phy, ntax.hybrid)
			if(is.na(focal.node)) {
				done=FALSE
				break()	
			}
			recipients <- phy$tip.label[getDescendants(phy, node=focal.node)]
			recipients <- recipients[!is.na(recipients)] #since we want just the tips
			recipient.ids <- which(phy$tip.label %in% recipients)
			used.recipients <- append(used.recipients, recipients)
		} else {
			focal.node<-sample(available.recipient.ids, 1, replace=FALSE)
			recipient.ids <- focal.node
			recipients <- phy$tip.label[focal.node]
			used.recipients <- append(used.recipients, recipients)
		}
		available.recipient.ids <- available.recipient.ids[!available.recipient.ids %in% recipient.ids]
		longest.from.root <- nodeheight(phy, node=focal.node)
		shortest.from.root <- nodeheight(phy, node=GetAncestor(phy, focal.node))
		all.heights <- nodeHeights(phy)
		#idea here: take a recipient clade. The gene flow must happen on its stem edge, which starts at shortest.from.root and goes up to longest.from.root. Gene flow can't go back in time
		qualifying.lower <- which(all.heights[,1]<longest.from.root) #if this is false, gene flow goes back in time
		qualifying.upper <- sequence(dim(all.heights)[1]) #in general, gene flow can go forward in time via ghost lineages
		if(!allow.ghost) {
			qualifying.upper <- which(all.heights[,2]>shortest.from.root) #if no ghost lineages, then there must be temporal overlap between the donor and recipient lineages. So the tipward end of the donor edge must be later than the rootward end of the recipient edge
		}
		qualifying.upper <- qualifying.upper[which(phy$edge[qualifying.upper,2]!=focal.node)] #let's not hybridize with ourselves
		qualifying.all <- qualifying.upper[qualifying.upper %in% qualifying.lower]
		if(length(qualifying.all)==0) {
			break()	
		}
		donor.edge <- sample(qualifying.all, 1)
		donors <- phy$tip.label[getDescendants(phy, phy$edge[donor.edge,2])]
		donors <- donors[!is.na(donors)] #getDescendants includes all descendant nodes, including internal ones. We just want the terminal taxa
		time.in <- runif(1, min=shortest.from.root, max=longest.from.root)
		time.out <- runif(1, min=all.heights[donor.edge,1], max=min(time.in, all.heights[donor.edge,2]))
		if (!allow.ghost) {
			time.in <- runif(1, min=max(shortest.from.root, all.heights[donor.edge,1]), max=min(longest.from.root, all.heights[donor.edge,2])) #if no ghost lineages, must move from the overlapping interval	
			time.out <- time.in
		}		
		pairs <- expand.grid(donors, recipients)
		for (pairs.index in sequence(dim(pairs)[1])) {
			flow <- rbind(flow, data.frame(donor=pairs[pairs.index,1], recipient=pairs[pairs.index,2], m=flow.proportion, time.from.root.donor=time.out, time.from.root.recipient=time.in, stringsAsFactors=FALSE))	
		}
		if(length(used.recipients)==ntax.hybrid) {
			done=TRUE
		}
	}
	flow$donor <- as.character(flow$donor)
	flow$recipient <- as.character(flow$recipient)
	flow$m <- as.numeric(as.character(flow$m))
	flow$time.from.root.donor <-as.numeric(as.character(flow$time.from.root.donor))
	flow$time.from.root.recipient <-as.numeric(as.character(flow$time.from.root.recipient))
	return(list(phy=phy, flow=flow))
}

PlotNetwork <- function(phy, flow, col.non="black", col.hybrid="red") {
	library(phylobase)
	phy<-reorder(phy, "pruningwise")
	phy4 <- as(phy, "phylo4")
	xxyy <- phyloXXYY(phy4)
	#plot(phy4)
	plot(x=c(min(xxyy$xx), 1.1*max(xxyy$xx)), y=range(xxyy$yy), type="n", xaxt="n", xlab="", yaxt="n", ylab="", bty="n")
	arrows(x0=xxyy$segs$v0x, x1=xxyy$segs$v1x, y0=xxyy$segs$v0y, y1=xxyy$segs$v1y, length=0)
	arrows(x0=xxyy$segs$h0x, x1=xxyy$segs$h1x, y0=xxyy$segs$h0y, y1=xxyy$segs$h1y, length=0)
	label.colors <- rep(col.non, Ntip(phy))
	for (i in sequence(Ntip(phy))) {
		if	(names(getNode(phy4, xxyy$torder))[i] %in% flow$recipient) {
			label.colors[i]<-col.hybrid	
		}
	}
	text(x=rep(max(xxyy$xx), Ntip(phy)), y=xxyy$yy[which(phylobase:::edges(phy4)[xxyy$eorder,2] %in% sequence(Ntip(phy)))], names(getNode(phy4, xxyy$torder)), col=label.colors, pos=4)
	for (i in sequence(dim(flow)[1])) {
		recipient.node <- getNode(phy4, flow$recipient[i])
		recipient.path <- c(recipient.node, ancestors(phy4, recipient.node))
		recipient.path.heights <- nodeDepth(phy4, recipient.path)
		valid.recipients <- recipient.path[which(recipient.path.heights > flow$time.from.root.recipient[i])]
		recipient <- valid.recipients[length(valid.recipients)] #do it from the earliest qualifying tipward node
		if(length(recipient.path)>1 && length(which(recipient.path.heights==flow$time.from.root.recipient[i]))>0) { #the latter condition means we're moving to an existing node
				recipient <- valid.recipients[length(valid.recipients)-1]
		}
		y1 <- xxyy$yy[which(phylobase:::edges(phy4)[xxyy$eorder,2] == recipient)]
		donor.node <- getNode(phy4, flow$donor[i])
		donor.path <- c(donor.node, ancestors(phy4, donor.node))
		donor.path.heights <- nodeDepth(phy4, donor.path)
		valid.donors <- donor.path[which(donor.path.heights > flow$time.from.root.donor[i])]
		donor <- valid.donors[length(valid.donors)] #do it from the earliest qualifying tipward node
		y0 <- xxyy$yy[which(phylobase:::edges(phy4)[xxyy$eorder,2] == donor)]
		arrows(x0=flow$time.from.root.donor[i]/max(branching.times(phy)), x1=flow$time.from.root.recipient[i]/max(branching.times(phy)), y1=y1, y0=y0, col="red") #rescale since it goes from zero to 1 in height
		#grid.arrows(x=c(flow$time.from.root[i],flow$time.from.root[i]), y=c(y0, y1))
	}
}

#		names(free.parameters) <- c("sigma.sq", "mu", "bt", "vh", "SE")
#GetMeansModified <- function(x, phy, flow, actual.params) {
#params must be named vector
SimulateTipData <- function(phy, flow, params) {
	library(mvtnorm)
	VCV.modified <- GetVModified(params, phy, flow, params)
	Means.modified <- GetMeansModified(params, phy, flow, params)
	tips <- rmvnorm(n=1, mean = Means.modified, sigma = VCV.modified)[1,]
	names(tips) <- rownames(VCV.modified)
	return(tips)
}
