newick2phylog<-function (x.tre, call = match.call()) 
{
  complete <- function(x.tre) {
    if (length(x.tre) > 1) {
      w <- ""
      for (i in 1:length(x.tre)) w <- paste(w, x.tre[i], 
                                            sep = "")
      x.tre <- w
    }
    ndroite <- nchar(gsub("[^)]", "", x.tre))
    ngauche <- nchar(gsub("[^(]", "", x.tre))
    if (ndroite != ngauche) 
      stop(paste(ngauche, "( versus", ndroite, ")"))
    if (regexpr(";", x.tre) == -1) 
      stop("';' not found")
    i <- 0
    kint <- 0
    kext <- 0
    arret <- FALSE
    if (regexpr("\\[", x.tre) != -1) {
      x.tre <- gsub("\\[[^\\[]*\\]", "", x.tre)
    }
    x.tre <- gsub(" ", "", x.tre)
    while (!arret) {
      i <- i + 1
      if (substr(x.tre, i, i) == ";") 
        arret <- TRUE
      if (substr(x.tre, i, i + 1) == "(,") {
        kext <- kext + 1
        add <- paste("Ext", kext, sep = "")
        x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, 
                                                              i + 1), sep = "")
        i <- i + 1
      }
      else if (substr(x.tre, i, i + 1) == ",,") {
        kext <- kext + 1
        add <- paste("Ext", kext, sep = "")
        x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, 
                                                              i + 1), sep = "")
        i <- i + 1
      }
      else if (substr(x.tre, i, i + 1) == ",)") {
        kext <- kext + 1
        add <- paste("Ext", kext, sep = "")
        x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, 
                                                              i + 1), sep = "")
        i <- i + 1
      }
      else if (substr(x.tre, i, i + 1) == "(:") {
        kext <- kext + 1
        add <- paste("Ext", kext, sep = "")
        x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, 
                                                              i + 1), sep = "")
        i <- i + 1
      }
      else if (substr(x.tre, i, i + 1) == ",:") {
        kext <- kext + 1
        add <- paste("Ext", kext, sep = "")
        x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, 
                                                              i + 1), sep = "")
        i <- i + 1
      }
      else if (substr(x.tre, i, i + 1) == "),") {
        kint <- kint + 1
        add <- paste("I", kint, sep = "")
        x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, 
                                                              i + 1), sep = "")
        i <- i + 1
      }
      else if (substr(x.tre, i, i + 1) == "))") {
        kint <- kint + 1
        add <- paste("I", kint, sep = "")
        x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, 
                                                              i + 1), sep = "")
        i <- i + 1
      }
      else if (substr(x.tre, i, i + 1) == "):") {
        kint <- kint + 1
        add <- paste("I", kint, sep = "")
        x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, 
                                                              i + 1), sep = "")
        i <- i + 1
      }
      else if (substr(x.tre, i, i + 1) == ");") {
        add <- "Root"
        x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, 
                                                              i + 1), sep = "")
        i <- i + 1
      }
    }
    lab.points <- strsplit(x.tre, "[(),;]")[[1]]
    lab.points <- lab.points[lab.points != ""]
    no.long <- (regexpr(":", lab.points) == -1)
    if (all(no.long)) {
      lab.points <- paste(lab.points, ":", c(rep("1", length(no.long) - 
        1), "0.0"), sep = "")
    }
    else if (no.long[length(no.long)]) {
      lab.points[length(lab.points)] <- paste(lab.points[length(lab.points)], 
                                              ":0.0", sep = "")
    }
    else if (any(no.long)) {
      
      stop("Non convenient ancdes.ancdes.array leaves or nodes with and without length")
    }
    w <- strsplit(x.tre, "[(),;]")[[1]]
    w <- w[w != ""]
    leurre <- make.names(w, unique = TRUE)
    leurre <- gsub("[.]", "_", leurre)
    for (i in 1:length(w)) {
      old <- paste(w[i])
      x.tre <- sub(old, leurre[i], x.tre, fixed = TRUE)
    }
    w <- strsplit(lab.points, ":")
    label <- function(x) {
      lab <- x[1]
      lab <- gsub("[.]", "_", lab)
      return(lab)
    }
    longueur <- function(x) {
      long <- x[2]
      return(long)
    }
    labels <- unlist(lapply(w, label))
    longueurs <- unlist(lapply(w, longueur))
    labels <- make.names(labels, TRUE)
    labels <- gsub("[.]", "_", labels)
    w <- labels
    for (i in 1:length(w)) {
      new <- w[i]
      x.tre <- sub(leurre[i], new, x.tre)
    }
    cat <- rep("", length(w))
    for (i in 1:length(w)) {
      new <- w[i]
      if (regexpr(paste("\\)", new, sep = ""), x.tre) != 
        -1) 
        cat[i] <- "int"
      else if (regexpr(paste(",", new, sep = ""), x.tre) != 
        -1) 
        cat[i] <- "ext"
      else if (regexpr(paste("\\(", new, sep = ""), x.tre) != 
        -1) 
        cat[i] <- "ext"
      else cat[i] <- "unknown"
    }
    return(list(tre = x.tre, noms = labels, poi = as.numeric(longueurs), 
                cat = cat))
  }
  res <- complete(x.tre)
  poi <- res$poi
  nam <- res$noms
  names(poi) <- nam
  cat <- res$cat
  res <- list(tre = res$tre)
  res$leaves <- poi[cat == "ext"]
  names(res$leaves) <- nam[cat == "ext"]
  res$nodes <- poi[cat == "int"]
  names(res$nodes) <- nam[cat == "int"]
  listclass <- list()
  dnext <- c(names(res$leaves), names(res$nodes))
  listpath <- as.list(dnext)
  names(listpath) <- dnext
  x.tre <- res$tre
  while (regexpr("[(]", x.tre) != -1) {
    a <- regexpr("\\([^\\(\\)]*\\)", x.tre)
    n1 <- a[1] + 1
    n2 <- n1 - 3 + attr(a, "match.length")
    chasans <- substring(x.tre, n1, n2)
    chaavec <- paste("\\(", chasans, "\\)", sep = "")
    nam <- unlist(strsplit(chasans, ","))
    w1 <- strsplit(x.tre, chaavec)[[1]][2]
    parent <- unlist(strsplit(w1, "[,\\);]"))[1]
    listclass[[parent]] <- nam
    x.tre <- gsub(chaavec, "", x.tre)
    w2 <- which(unlist(lapply(listpath, function(x) any(x[1] == 
      nam))))
    for (i in w2) {
      listpath[[i]] <- c(parent, listpath[[i]])
    }
  }
  res$parts <- listclass
  res$paths <- listpath
  dnext <- c(res$leaves, res$nodes)
  names(dnext) <- c(names(res$leaves), names(res$nodes))
  res$droot <- unlist(lapply(res$paths, function(x) sum(dnext[x])))
  res$call <- call
  class(res) <- "phylog"
  #if (!add.tools) 
  return(res)
  #return(newick2phylog.addtools(res))
}

comb<-function(n,r){
  return(factorial(n)/(factorial(n-r)*factorial(r)))
}


pair_fcn<-function(tmp){ # return pair for "tmp" sequences.
  numl=comb(length(tmp),2)   
  count=0
  posit<-array(0,c(numl))
  for(i in 1:length(tmp)){
    for(j in 1:length(tmp)){
      if(i<j){
        count=count+1
        posit[count]=paste(c(tmp[i]),c(tmp[j]),sep=",")
      }
    }
  }  
  return(posit)
}# end of pair_fcn.


pair_array<-function(tmp){#generate pair in array format. 
  pair <-pair_fcn(tmp)
  p_arr<-matrix(c(0),nrow=length(pair),ncol=2)
  for(i in 1:length(pair)){
    p_arr[i,1]=unlist(strsplit(pair[i],","))[1]
    p_arr[i,2]=unlist(strsplit(pair[i],","))[2]
  }
  return(p_arr)
}# end of function pair_array.



ord_fcn<-function(res,tipnames){ #this function gets the acenstor-descendants relationship.
  bmtp<-matrix(rev(res$nde),ncol=1) 
  rvpt <-rev((res$parts))
  rept<-array(0,c(length(rvpt),2))
  for(i in 1:length(rvpt)){
    rept[i,]=unlist(rvpt[i])
  }           
  cmb<-cbind(bmtp,rept)
  brnlen<-res$droot[(length(tipnames)+1):length(res$droot)]
  root<-matrix(cmb[1,],nrow=1)
  cmb<-cmb[-1,]
  brnlen<-brnlen[1:(length(brnlen)-1)]
  new_ord<-order(brnlen,decreasing=TRUE)
  cmb<-cmb[new_ord,]
  cmb<-rbind(root,cmb)
  return(cmb)
}# end of function ord_fcn.

getntn<-function(res){# this function gets rid of unnecessarily "_" symbol. 
  size<-length(res$parts)
  relarr<-array(0,c(size,3))
  rvpt <-(res$parts)
  rept<-array(0,c(length(rvpt),2))
  for(i in 1:length(rvpt)){
    rept[i,]=unlist(rvpt[i])
  }
  for(i in 1:size){
    relarr[i,1]<-names(res$parts)[i]
  }
  relarr[,2:3]<-rept
  temp<-matrix(0,row<-size)
  
  for(j in 2:3){
    for (i in 1: size){ 
      stmp<-unlist(strsplit(relarr[,j][i], "_" ))
      temp[i]<-stmp[1]
    }
    relarr[,j]<-temp
  }
  
  ndlen<- res$droot[!(res$droot==1)]
  nam<-names(ndlen)
  ck1<-array(0,c(length(nam)))
  count<-0
  for (ele in c(nam)){
    count<-count+1
    len <- length( unlist(strsplit(ele ,"_" )))
    if( len==2 ){ck1[count]<-1}
  }
  ndlen<-ndlen[!ck1]
  new_ord<-order(ndlen)
  relarr<-relarr[new_ord,]
  
  return(relarr)
}# end of function getntn.

getbrnlen<-function(res){#this function is used to obtain branch length.
  ndlen<- res$droot[!(res$droot==1)]
  
  nam<-names(ndlen)
  
  ck1<-array(0,c(length(nam)))
  count<-0
  for (ele in c(nam)){
    count<-count+1
    len <- length( unlist(strsplit(ele ,"_" )))
    if( len==2 ){ck1[count]<-1}
  }
  ndlen<-ndlen[!ck1]
  
  ndlen<-sort(ndlen)
  ck2<-array(0,c(length(ndlen)))
  for(i in 1:(length(ndlen)-1)){
    if(abs(ndlen[i]-ndlen[i+1])<10^(-5)){ck2[i]=1}
  }
  
  ndlen<-ndlen[!ck2]
  
  
  brnlen<-array(0,c(length(ndlen)))
  tmplen<-ndlen
  
  for(i in 1:(length(brnlen)-1)){
    brnlen[i]<-tmplen[i+1]-tmplen[i]
  }
  brnlen[length(brnlen)] <-1-tmplen[(length(tmplen))]
  return(brnlen)
}# end of function getbrnlen.

### The cov_mtx function will return the covaraince matrix. 
cov_mtx<-function(x,branchlength=branchlength,ancdes.array=ancdes.array,nleaves=nleaves,tipnames=tipnames,model.Index=model.Index){    #now it is a function of bt,  h, sigma^2 and sigma_H^2
  bt<-x[1]
  h<-x[2]
  sigma_sq<-x[3]
  sigma.H_sq<-x[4]
  
  h<-1/2
  if(model.Index==1){bt<-1}
  if(model.Index==2){sigma.H_sq<-0}
  if(model.Index==3){bt<-1;sigma.H_sq<-0}
  
  ins_fcn<-function(ist,sqc){#finds position to insert between two parents, for hybrdization only.
    ist<-as.numeric(unlist(strsplit(ist,"X"))[2])
    arr<-array(0,c(length(otmp)))
    for(i in 1:length(arr)){
      arr[i]<-as.numeric(unlist(strsplit(sqc[i],"X"))[2])  
    }
    insp<-which(arr==(ist-1))+1
    return(insp)
  }
  var_fcn<-function(){#return the variance.
    for(i in 1:length(otmp)){#use to fill other diagonal. 
      newi<-which(rownames(mtx)%in%otmp[i])              
      oldi<-which(rownames(omtx)%in%otmp[i])
      mtx[newi,newi]<-omtx[oldi,oldi]
    }#fill in old value from omtx exclude the new hyd.
    
    prn1<-tmp[which(tmp%in%ins)-1]#grab elements.
    prn2<-tmp[which(tmp%in%ins)+1]
    prn1<-which(rownames(omtx) %in% prn1)#grab position according to prn1.
    prn2<-which(rownames(omtx) %in% prn2)
    
    
    vhii<- bt^2*h^2*omtx[prn1,prn1]+bt^2*(1-h)^2*omtx[prn2,prn2]+2*bt^2*h*(1-h)*omtx[prn1,prn2] 
    
    hii<-which(!(tmp %in% otmp))#use to insert variance for hyd.
    mtx[hii,hii]<-vhii      #fill in the diagonal hyd. 
    return(mtx)
  }#formula for insertion hyd variance.
  
  
  fillspcmtx<-function(){#fill matrix due to sepciation.       
    elm<-function(){ #use to cut one row of the pair array which the speciation happens.
      ck<-c(tmp[nsi],tmp[nsj])
      for(i in 1:dim(pn_arr)[1]){
        if(sum(pn_arr[i,]==ck)==2){break}
      }
      return(i)}
    
    pn_arr<-pair_array(tmp)
    po_arr<-pair_array(otmp)
    
    #search new speciate position.
    nsi<-which(!(tmp %in% otmp))[1]
    nsj<-which(!(tmp %in% otmp))[2]
    osii<-which(!(otmp %in% tmp))
    mtx[nsi,nsj]<- omtx[osii,osii]
    #Fill in value: the covariance for 2 speciated species equal the variance of the parent.
    
    pn_arr<-pn_arr[-elm(),]#delete the ancdes.array that is already used.
    
    #The following fills covaraince components by the previous matrix.
    while(length(pn_arr[,1])>0){
      newi<-which(rownames(mtx) %in% pn_arr[1,1])
      newj<-which(rownames(mtx) %in% pn_arr[1,2])
      
      if( tmp[nsi] %in% pn_arr[1,]){
        otg<-which(!(pn_arr[1,] %in%  tmp[nsi]))
        oldi<- which( rownames(omtx) %in% otmp[osii])
        oldj<-which(rownames(omtx) %in% pn_arr[1,otg])
      }
      
      if( tmp[nsj] %in% pn_arr[1,] ){
        otg<-which(!(pn_arr[1,] %in%  tmp[nsj]))
        oldi<- which( rownames(omtx) %in% otmp[osii])
        oldj<-which(rownames(omtx) %in% pn_arr[1,otg])
      }
      
      if(!(tmp[nsi] %in% pn_arr[1,]) && !(tmp[nsj] %in% pn_arr[1,])){
        #detect common between omtx and mtx.   
        oldi<-which(rownames(omtx) %in% pn_arr[1,1])
        oldj<-which(rownames(omtx) %in% pn_arr[1,2])
      }
      mtx[newi,newj]<-omtx[oldi,oldj]
      pn_arr<-pn_arr[-1,]#delete row. 
      if(length(pn_arr)==2){pn_arr<-matrix(pn_arr,nrow=1)}
    }#end of while loop.
    
    mtx<-mtx+t(mtx)
    
    mtx[nsi,nsi]<-omtx[osii,osii]+ branchlength[length(tmp)-1]
    mtx[nsj,nsj]<-omtx[osii,osii]+ branchlength[length(tmp)-1]
    dianew<-which(tmp %in% otmp )
    diaold<-which(otmp %in% tmp )
    for(i in 1:length(dianew)){
      mtx[dianew[i],dianew[i]]<-omtx[diaold[i],diaold[i]]+branchlength[length(tmp)-1]
    }
    return(mtx)
  }#end of fillspcmtx.
  
  fillhydmtx<-function(){#fill in value into matrix due to hybridzation.   
    pn_arr<-pair_array(tmp)
    
    while(length(pn_arr[,1])>0){
      newi<-which(rownames(mtx) %in% pn_arr[1,1])
      newj<-which(rownames(mtx) %in% pn_arr[1,2])
      if (ins %in% pn_arr[1,]){#ins is the hybridized node. 
        otg<-pn_arr[1,which(!(pn_arr[1,] %in% ins ))]
        otgj<-which(rownames(omtx) %in% otg)
        #the other guy, could be the hybrdized nodes parent or others.
        #find the parent of ins.
        
        prn1<-tmp[which(tmp%in%ins)-1]#grab element.
        prn2<-tmp[which(tmp%in%ins)+1]        
        prn1<-which(rownames(omtx) %in% prn1)#grab position.
        prn2<-which(rownames(omtx) %in% prn2)
        
        mtx[newi,newj]<-bt*h*omtx[prn1,otgj] +bt*(1-h)*omtx[prn2,otgj] # cov(X, bt*hX+bt*(1-h)Y) we are going to use h=1/2.
        
        
      }else{#this is not hyd node, just read from previous mtx.
        #only need to read from rownames().
        oldi<-which(rownames(omtx) %in% pn_arr[1,1])
        oldj<-which(rownames(omtx) %in% pn_arr[1,2])
        mtx[newi,newj]<-omtx[oldi,oldj]
      }#end of else loop .
      pn_arr<-pn_arr[-1,] # delete ancdes.array array after using it.
      if(length(pn_arr)==2){pn_arr<-matrix(pn_arr,nrow=1)}
    }#end of while loop.
    return(mtx)
  }#end of fillhydmtx.
  
  #THE MAIN PROGRAM for covariance matrix.
  ckins<-FALSE # use to check the hybrdized event.
  rept<-ancdes.array[,2:3]# the descedant nodes.
  bmtp<-matrix((ancdes.array)[,1],ncol=1) #the acenstor node.
  
  loop<-2
  tmp=array(0,c(loop))
  if(loop==2){tmp=rept[1,]
              otmp<-tmp
              mtx<-diag(branchlength[1],c(length(tmp)))
              rownames(mtx)<-c(tmp)
              colnames(mtx)<-c(tmp)
              omtx<-mtx
  }#end of loop==2 
  while(loop<length(bmtp)){#loaded the acenstor-descendant ancdes.array. 
    loop<-loop+1#use loop to use the ancdes.array
    tmp=array(0,c(length(otmp)+1))#the new seq.
    mtx<-matrix(0,nrow=length(tmp),ncol=length(tmp))
    q=loop-1#index for matching the right element: will use below. 
    op<-which(otmp==bmtp[q])#index for insertion position.
    if(length(op)!=0){#op!=0 means that weve  detected speciation.
      tmp[op:(op+1)]=rept[q,] #insertion the new speciation species.
      
      if(op==1){tmp[(op+2):length(tmp)]=otmp[(op+1):(length(tmp)-1)]}
      if((op+1)==length(tmp)){tmp[1:(op-1)]=otmp[1:(op-1)] }
      if(op!=1 && (op+1)!=length(tmp)){
        tmp[(op+2):length(tmp)]=otmp[(op+1):(length(tmp)-1)] 
        tmp[1:(op-1)]=otmp[1:(op-1)]}
      
      
      rownames(mtx)<-c(tmp)
      colnames(mtx)<-c(tmp)
      mtx<- fillspcmtx()
      otmp<-tmp
      omtx<-mtx
      #above generate sequence and cov. matrix for speciation.           
      
    }else{#  op = 0 means that we have detected the hybridize event.
      ins<-(bmtp[q])#grab the insertion element, ins will be used in the fillhydmtx function.
      insp<-ins_fcn(ins,otmp)#catch the position for insertion.
      tmp[insp]<-ins #insert the hyd element.
      tmp[(insp+1):length(tmp)]=otmp[insp:(length(tmp)-1)]
      tmp[1:(insp-1)]=otmp[1:(insp-1)]
      rownames(mtx)<-c(tmp)
      colnames(mtx)<-c(tmp)
      diamtx<-var_fcn()
      mtx<- fillhydmtx()
      mtx<-mtx+t(mtx)+diamtx
      #fill in the diagonal elements.
      
      otmp<-tmp
      omtx<-mtx
      #above generate the sequnce and fill in the value into matrix for hybrdization.      
      
      ckins<-TRUE #since we did an insertion, the next step is to replace 3 elements. 
    }#end of the length(op)!=0 if-else. 
    
    if(ckins){#replace 3 elements in tmp sequence.  
      tmp<-array(0,c(length(tmp)))
      tmp[which(otmp==ins)]<- rept[loop-1,1] # replaced with hyd element.
      tmp[which(otmp == bmtp[loop])] = rept[loop,which(rept[loop,]!=ins)]
      tmp[which(otmp == bmtp[loop+1])] = rept[loop+1,which(rept[loop+1,]!=ins)]            
      #replace 3 new nodes.
      tx1<-which(otmp==ins)
      tx2<-which(otmp == bmtp[loop])
      tx3<-which(otmp == bmtp[loop+1])
      for(i in 1:length(tmp)){
        if (i != tx1 && i!=tx2 && i!=tx3){
          tmp[i]=otmp[i]
        }
      }
      
      otmp<-tmp      
      rownames(mtx)<-c(tmp)
      colnames(mtx)<-c(tmp)
      
      mtx<-mtx+diag(branchlength[length(tmp)-1],c(length(tmp)) )
      
      omtx<-mtx
      ckins<-FALSE
      loop<-loop+2          
    }#end of replace 3 elements
  }#end of while loop
  
  if(sum(tipnames%in%tmp)!=nleaves){#catches the last speciation event.
    tmp<-tipnames
    mtx<-matrix(0,nrow=length(tmp),ncol=length(tmp))
    rownames(mtx)<-c(tmp)
    colnames(mtx)<-c(tmp)
    mtx<-fillspcmtx()
  }#end of if (sum(tipnames%in%tmp)!=nleaves).
  
  
  
  if(dim(mtx)[1]==22){ #for cichlid data only
  del1<-which(rownames(mtx)=="a")
  mtx<-mtx[-del1,-del1]
  del2<-which(rownames(mtx)=="b")
  mtx<-mtx[-del2,-del2]
  del3<-which(rownames(mtx)=="c")   
  mtx<-mtx[-del3,-del3]
  
  srd<-c("X1","X2", "X3","X4","X5","X6","X7","X8", "X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19")
  rnmt<-rownames(mtx)
  newmtx<-array(0,c(dim(mtx)))
  rownames(newmtx)<-srd
  colnames(newmtx)<-srd
  for( i in 1: 19){
    for(j in 1:19){
      newmtx[i,j]= mtx[ which(rnmt==srd[i]), which(rnmt==srd[j]) ]
    }
  }  
  mtx<-newmtx
  }#end of if n==19
  
  
  
  mtx<-mtx*sigma_sq   
  
  hybrid.Index<-hybrid.node(ancdes.array,nleaves) #extra burst for hybrid
  for(i in hybrid.Index){
    mtx[i,i]<-mtx[i,i]+sigma.H_sq
  }
  
  if(model.Index==1){rm(bt)}
  if(model.Index==2){rm(sigma.H_sq)}
  if(model.Index==3){rm(bt);rm(sigma.H_sq)}
  return(mtx)
}#end of cov_mtx 


NegLogLike<-function(x,Y=Y,n=n,branchlength=branchlength,ancdes.array=ancdes.array,nleaves=nleaves,tipnames=tipnames,model.Index=model.Index){
  badval<-(0.5)*.Machine$double.xmax
  mu<-x[1]
  bt<-x[2]
  h<-x[3]
  sigma_sq<-x[4]
  sigma.H_sq<-x[5]
  
  h<-1/2
  if(model.Index==1){bt<-1}
  if(model.Index==2){sigma.H_sq<-0}
  if(model.Index==3){bt<-1;sigma.H_sq<-0}
  
  W <- cov_mtx(c(bt,h,sigma_sq,sigma.H_sq),branchlength=branchlength,ancdes.array=ancdes.array,nleaves=nleaves,tipnames=tipnames,model.Index=model.Index) #NOTE TO BCO: CHECK THAT t IS BEING USED CORRECTLY
  hybrid.Index<-hybrid.node(ancdes.array,nleaves) #extra burst for hybrid
  muone<- mu*matrix(1,nrow=n)
  for(i in hybrid.Index){
    muone[i]<-bt*muone[i] 
  }   
  
  NegLogML <- n/2*log(2*pi)+1/2*t(Y-muone)%*%pseudoinverse(W)%*%(Y-muone) + 1/2*log(abs(det(W))) 
  
  if(min(W)<0 || h<0 || h>1 || sigma_sq <0 || sigma.H_sq<0 ||  bt <= 0.0000001) {
    NegLogML<-badval 
  }
  
  if(model.Index==1){rm(bt)}
  if(model.Index==2){rm(sigma.H_sq)}
  if(model.Index==3){rm(bt);rm(sigma.H_sq)}
  
  return(NegLogML[1]) #need to put this in to get scalar output
}#end of NegLogLike.




Hessian.cov_mtx<-function(x,branchlength=branchlength,ancdes.array=ancdes.array,nleaves=nleaves,tipnames=tipnames,model.Index=model.Index){    #now it is a function of bt,  h, sigma^2 and sigma_H^2
  bt<-x[1]
  h<-x[2]
  sigma_sq<-x[3]
  sigma.H_sq<-x[4]
  
  
  
  ins_fcn<-function(ist,sqc){#finds position to insert between two parents, for hybrdization only.
    ist<-as.numeric(unlist(strsplit(ist,"X"))[2])
    arr<-array(0,c(length(otmp)))
    for(i in 1:length(arr)){
      arr[i]<-as.numeric(unlist(strsplit(sqc[i],"X"))[2])  
    }
    insp<-which(arr==(ist-1))+1
    return(insp)
  }
  var_fcn<-function(){#return the variance.
    for(i in 1:length(otmp)){#use to fill other diagonal. 
      newi<-which(rownames(mtx)%in%otmp[i])              
      oldi<-which(rownames(omtx)%in%otmp[i])
      mtx[newi,newi]<-omtx[oldi,oldi]
    }#fill in old value from omtx exclude the new hyd.
    
    prn1<-tmp[which(tmp%in%ins)-1]#grab elements.
    prn2<-tmp[which(tmp%in%ins)+1]
    prn1<-which(rownames(omtx) %in% prn1)#grab position according to prn1.
    prn2<-which(rownames(omtx) %in% prn2)
    
    vhii<- bt^2*h^2*omtx[prn1,prn1]+bt^2*(1-h)^2*omtx[prn2,prn2]+2*bt^2*h*(1-h)*omtx[prn1,prn2] 
    
    hii<-which(!(tmp %in% otmp))#use to insert variance for hyd.
    mtx[hii,hii]<-vhii      #fill in the diagonal hyd. 
    return(mtx)
  }#formula for insertion hyd variance.
  
  
  fillspcmtx<-function(){#fill matrix due to sepciation.       
    elm<-function(){ #use to cut one row of the pair array which the speciation happens.
      ck<-c(tmp[nsi],tmp[nsj])
      for(i in 1:dim(pn_arr)[1]){
        if(sum(pn_arr[i,]==ck)==2){break}
      }
      return(i)}
    
    pn_arr<-pair_array(tmp)
    po_arr<-pair_array(otmp)
    
    #search new speciate position.
    nsi<-which(!(tmp %in% otmp))[1]
    nsj<-which(!(tmp %in% otmp))[2]
    osii<-which(!(otmp %in% tmp))
    mtx[nsi,nsj]<- omtx[osii,osii]
    #Fill in value: the covariance for 2 speciated species equal the variance of the parent.
    
    pn_arr<-pn_arr[-elm(),]#delete the ancdes.array that is already used.
    
    #The following fills covaraince components by the previous matrix.
    while(length(pn_arr[,1])>0){
      newi<-which(rownames(mtx) %in% pn_arr[1,1])
      newj<-which(rownames(mtx) %in% pn_arr[1,2])
      
      if( tmp[nsi] %in% pn_arr[1,]){
        otg<-which(!(pn_arr[1,] %in%  tmp[nsi]))
        oldi<- which( rownames(omtx) %in% otmp[osii])
        oldj<-which(rownames(omtx) %in% pn_arr[1,otg])
      }
      
      if( tmp[nsj] %in% pn_arr[1,] ){
        otg<-which(!(pn_arr[1,] %in%  tmp[nsj]))
        oldi<- which( rownames(omtx) %in% otmp[osii])
        oldj<-which(rownames(omtx) %in% pn_arr[1,otg])
      }
      
      if(!(tmp[nsi] %in% pn_arr[1,]) && !(tmp[nsj] %in% pn_arr[1,])){
        #detect common between omtx and mtx.   
        oldi<-which(rownames(omtx) %in% pn_arr[1,1])
        oldj<-which(rownames(omtx) %in% pn_arr[1,2])
      }
      mtx[newi,newj]<-omtx[oldi,oldj]
      pn_arr<-pn_arr[-1,]#delete row. 
      if(length(pn_arr)==2){pn_arr<-matrix(pn_arr,nrow=1)}
    }#end of while loop.
    
    mtx<-mtx+t(mtx)
    
    mtx[nsi,nsi]<-omtx[osii,osii]+ branchlength[length(tmp)-1]
    mtx[nsj,nsj]<-omtx[osii,osii]+ branchlength[length(tmp)-1]
    dianew<-which(tmp %in% otmp )
    diaold<-which(otmp %in% tmp )
    for(i in 1:length(dianew)){
      mtx[dianew[i],dianew[i]]<-omtx[diaold[i],diaold[i]]+branchlength[length(tmp)-1]
    }
    return(mtx)
  }#end of fillspcmtx.
  
  fillhydmtx<-function(){#fill in value into matrix due to hybridzation.   
    pn_arr<-pair_array(tmp)
    
    while(length(pn_arr[,1])>0){
      newi<-which(rownames(mtx) %in% pn_arr[1,1])
      newj<-which(rownames(mtx) %in% pn_arr[1,2])
      if (ins %in% pn_arr[1,]){#ins is the hybridized node. 
        otg<-pn_arr[1,which(!(pn_arr[1,] %in% ins ))]
        otgj<-which(rownames(omtx) %in% otg)
        #the other guy, could be the hybrdized nodes parent or others.
        
        #find the parent of ins.
        prn1<-tmp[which(tmp%in%ins)-1]#grab element.
        prn2<-tmp[which(tmp%in%ins)+1]        
        prn1<-which(rownames(omtx) %in% prn1)#grab position.
        prn2<-which(rownames(omtx) %in% prn2)
        
        ###########################################        
        mtx[newi,newj]<-bt*h*omtx[prn1,otgj] +bt*(1-h)*omtx[prn2,otgj] # cov(X, bt*hX+bt*(1-h)Y) we are going to use h=1/2.
        ############################################ 
        
      }else{#this is not hyd node, just read from previous mtx.
        #only need to read from rownames().
        oldi<-which(rownames(omtx) %in% pn_arr[1,1])
        oldj<-which(rownames(omtx) %in% pn_arr[1,2])
        mtx[newi,newj]<-omtx[oldi,oldj]
      }#end of else loop .
      pn_arr<-pn_arr[-1,] # delete ancdes.array array after using it.
      if(length(pn_arr)==2){pn_arr<-matrix(pn_arr,nrow=1)}
    }#end of while loop.
    return(mtx)
  }#end of fillhydmtx.
  
  #THE MAIN PROGRAM for covariance matrix.
  ckins<-FALSE # use to check the hybrdized event.
  rept<-ancdes.array[,2:3]# the descedant nodes.
  bmtp<-matrix((ancdes.array)[,1],ncol=1) #the acenstor node.
  
  loop<-2
  tmp=array(0,c(loop))
  if(loop==2){tmp=rept[1,]
              otmp<-tmp
              mtx<-diag(branchlength[1],c(length(tmp)))
              rownames(mtx)<-c(tmp)
              colnames(mtx)<-c(tmp)
              omtx<-mtx
  }#end of loop==2 
  while(loop<length(bmtp)){#loaded the acenstor-descendant ancdes.array. 
    loop<-loop+1#use loop to use the ancdes.array
    tmp=array(0,c(length(otmp)+1))#the new seq.
    mtx<-matrix(0,nrow=length(tmp),ncol=length(tmp))
    q=loop-1#index for matching the right element: will use below. 
    op<-which(otmp==bmtp[q])#index for insertion position.
    if(length(op)!=0){#op!=0 means that weve  detected speciation.
      tmp[op:(op+1)]=rept[q,] #insertion the new speciation species.
      
      if(op==1){tmp[(op+2):length(tmp)]=otmp[(op+1):(length(tmp)-1)]}
      if((op+1)==length(tmp)){tmp[1:(op-1)]=otmp[1:(op-1)] }
      if(op!=1 && (op+1)!=length(tmp)){
        tmp[(op+2):length(tmp)]=otmp[(op+1):(length(tmp)-1)] 
        tmp[1:(op-1)]=otmp[1:(op-1)]}
      
      
      rownames(mtx)<-c(tmp)
      colnames(mtx)<-c(tmp)
      mtx<- fillspcmtx()
      otmp<-tmp
      omtx<-mtx
      
      #above generate sequence and cov. matrix for speciation.           
      
    }else{#  op = 0 means that we have detected the hybridize event.
      ins<-(bmtp[q])#grab the insertion element, ins will be used in the fillhydmtx function.
      insp<-ins_fcn(ins,otmp)#catch the position for insertion.
      tmp[insp]<-ins #insert the hyd element.
      tmp[(insp+1):length(tmp)]=otmp[insp:(length(tmp)-1)]
      tmp[1:(insp-1)]=otmp[1:(insp-1)]
      rownames(mtx)<-c(tmp)
      colnames(mtx)<-c(tmp)
      diamtx<-var_fcn()
      mtx<- fillhydmtx()
      mtx<-mtx+t(mtx)+diamtx
      #fill in the diagonal elements.
      
      otmp<-tmp
      omtx<-mtx
      #above generate the sequnce and fill in the value into matrix for hybrdization.      
      
      ckins<-TRUE #since we did an insertion, the next step is to replace 3 elements. 
    }#end of the length(op)!=0 if-else. 
    
    if(ckins){#replace 3 elements in tmp sequence.  
      tmp<-array(0,c(length(tmp)))
      tmp[which(otmp==ins)]<- rept[loop-1,1] # replaced with hyd element.
      tmp[which(otmp == bmtp[loop])] = rept[loop,which(rept[loop,]!=ins)]
      tmp[which(otmp == bmtp[loop+1])] = rept[loop+1,which(rept[loop+1,]!=ins)]            
      #replace 3 new nodes.
      tx1<-which(otmp==ins)
      tx2<-which(otmp == bmtp[loop])
      tx3<-which(otmp == bmtp[loop+1])
      for(i in 1:length(tmp)){
        if (i != tx1 && i!=tx2 && i!=tx3){
          tmp[i]=otmp[i]
        }
      }
      
      otmp<-tmp      
      rownames(mtx)<-c(tmp)
      colnames(mtx)<-c(tmp)
      
      mtx<-mtx+diag(branchlength[length(tmp)-1],c(length(tmp)) )
      
      omtx<-mtx
      ckins<-FALSE
      loop<-loop+2          
    }#end of replace 3 elements
  }#end of while loop
  
  if(sum(tipnames%in%tmp)!=nleaves){#catches the last speciation event.
    tmp<-tipnames
    mtx<-matrix(0,nrow=length(tmp),ncol=length(tmp))
    rownames(mtx)<-c(tmp)
    colnames(mtx)<-c(tmp)
    mtx<-fillspcmtx()
  }#end of if (sum(tipnames%in%tmp)!=nleaves).
  
  if(dim(mtx)[1]==22){ #for cichlid data only
  del1<-which(rownames(mtx)=="a")
  mtx<-mtx[-del1,-del1]
  del2<-which(rownames(mtx)=="b")
  mtx<-mtx[-del2,-del2]
  del3<-which(rownames(mtx)=="c")   
  mtx<-mtx[-del3,-del3]
  
  srd<-c("X1","X2", "X3","X4","X5","X6","X7","X8", "X9","X10","X11","X12","X13","X14","X15","X16","X17","X18","X19")
  rnmt<-rownames(mtx)
  newmtx<-array(0,c(dim(mtx)))
  rownames(newmtx)<-srd
  colnames(newmtx)<-srd
  for( i in 1: 19){
    for(j in 1:19){
      newmtx[i,j]= mtx[ which(rnmt==srd[i]), which(rnmt==srd[j]) ]
    }
  }  
  mtx<-newmtx
  }#end of if n==19
  
  
  
  mtx<-mtx*sigma_sq   
  
  hybrid.Index<-hybrid.node(ancdes.array,nleaves) #extra burst for hybrid
  for(i in hybrid.Index){
    mtx[i,i]<-mtx[i,i]+sigma.H_sq
    #print(mtx[i,i])
  }
  
  return(mtx)
}#end of cov_mtx 


hybrid.node<-function(ancdes.array,nleaves){
  hyd.sigma_h<-c()
  for (i in 1:dim(ancdes.array)[1]){
    if(ancdes.array[i,2]==ancdes.array[i,3]){
      hyd.idx<-as.numeric(unlist(strsplit(ancdes.array[i,2],"X"))[2])
      if ( hyd.idx<=nleaves){						
        hyd.sigma_h<-c(hyd.sigma_h,hyd.idx )}else{
          
          hyd.speciation<-which(ancdes.array[,1]==ancdes.array[i,2])
          candidate.hyd.des<-as.numeric(unlist(strsplit(ancdes.array[hyd.speciation,],"X")[2:3]))[c(2,4)]
          for(hyd.des in candidate.hyd.des){	
            if(hyd.des<=nleaves){
              hyd.sigma_h<-c(hyd.sigma_h,hyd.des	)}}
          
        }
      
    }
  }
  return(hyd.sigma_h )
}


Hessian.NegLogLike<-function(x,Y=Y,n=n,branchlength=branchlength,ancdes.array=ancdes.array,nleaves=nleaves,tipnames=tipnames,model.Index=model.Index){
  badval<-(0.5)*.Machine$double.xmax
  mu<-x[1]
  bt<-x[2]
  h<-x[3]
  sigma_sq<-x[4]
  sigma.H_sq<-x[5]
  
  W <- Hessian.cov_mtx(c(bt,h,sigma_sq,sigma.H_sq),branchlength=branchlength,ancdes.array=ancdes.array,nleaves=nleaves,tipnames=tipnames,model.Index=model.Index) #NOTE TO BCO: CHECK THAT t IS BEING USED CORRECTLY
  hybrid.Index<-hybrid.node(ancdes.array,nleaves) #extra burst for hybrid
  muone<- mu*matrix(1,nrow=n)
  for(i in hybrid.Index){
    muone[i]<-bt*muone[i] 
  }   
  
  NegLogML <- n/2*log(2*pi)+1/2*t(Y-muone)%*%pseudoinverse(W)%*%(Y-muone) + 1/2*log(abs(det(W))) 
  
  return(NegLogML[1]) #need to put this in to get scalar output
}#end of NegLogLike.



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

#Para.Var.Matrix<-
#function(k){
#if(k==1){return(cov.set.bt.1.s.H.0)}
#if(k==2){return(cov.set.bt.1)}
#if(k==3){return(cov.set.s.H.0)}
#if(k==4){return(cov.set.free)}
#}

se.ave.weight.para<-
  function(para.vect,var.para.vect,weight.vect){
    num.para<-length(para.vect)
    sum.para<-0
    if(num.para>=1){
      for(i in 1:num.para){
        sum.para<-sum.para+weight.vect[i]*sqrt(var.para.vect[i] + ( para.vect[i]- mean (para.vect) )^2)  
      }
      return(sum.para)
    }else{return(0)}
  }



####################################################################
###################### MAIN PROGRAM ################################
####################################################################

BMhyd<-function(data, network){
  Y<-data
  x.tre<-network
  res<-newick2phylog(x.tre)
  ancdes.array<-getntn(res)
  branchlength<-getbrnlen(res)
  tipnames<-sort(names(res$droot[which(res$droot==1)]))
  nleaves<-length(tipnames)
  if(sum(is.na(Y))!=0){
  n<-nleaves-sum(is.na(Y))
  Y<-Y[-which(is.na(Y))]}else{n<-nleaves}
  output.array<-array(0,c(4,7))
  rownames(output.array)<-c("bt=1","v.H=0","bt=1;v.H=0","free")
  opt.method<-c("Nelder-Mead")
  Hessian.mtx<-array(0,c(4,5,5))
  para.cov.matrix<-array(0,c(4,5,5))
  p0 = c(mean(Y),1,1/2,var(Y),var(Y))#starting point
  cat("Begin optimx optimization routine -- 
    Starting value: (mu,beta,h,sigma^2,v.H)=(",round(p0,digits=2),")","\n\n")
    convergence.record<-array(0,dim(output.array)[1])
  for (model.Index in 1:dim(output.array)[1]){
    if(model.Index==1){k<-3;bt<-1}
    if(model.Index==2){k<-3;sigma.H_sq<-0}
    if(model.Index==3){k<-2;bt<-1;sigma.H_sq<-0}
    if(model.Index==4){k<-4}
    MLE.ALL<-optim(p0,NegLogLike,method="Nelder-Mead", Y=Y,n=n,branchlength=branchlength,ancdes.array=ancdes.array,nleaves=nleaves,tipnames=tipnames,model.Index=model.Index )
    print(MLE.ALL)
    convergence.record[model.Index]<- MLE.ALL$convergence
    if(MLE.ALL$convergence==0){cat("The MLE estimations converge","\n\n")}
    output.array[model.Index,1:5]<-MLE.ALL$par
    output.array[model.Index,6]<-MLE.ALL$value
    output.array[model.Index,7]<-AICc(n,k,MLE.ALL$value)
    if(model.Index==1){rm(k);rm(bt)}
    if(model.Index==2){rm(k);rm(sigma.H_sq)}
    if(model.Index==3){rm(k);rm(bt);rm(sigma.H_sq)}
    
  }#end of for model.Index
  
  obj<-NULL
  output.array<-cbind(output.array, matrix(output.array[,7]-min(output.array[,7]),ncol=1))
  output.array[,3]<-0.5
  output.array[1,2]<-1;output.array[2,5]<-0
  output.array[3,2]<-1;output.array[3,5]<-0;
  
  output.array<-cbind(output.array, matrix(AkaikeWeight(output.array[,8]),ncol=1) )
  output.array<-cbind(output.array, matrix(cumsum(output.array[,9]),ncol=1))
  
  cat("Finished MLE Estimation","\n\n")
  cat("Calculating standard error for MLEs","\n\n")
  for (model.Index in 1:dim(output.array)[1]){
    Hessian.mtx[model.Index,,]<-hessian(Hessian.NegLogLike, c(output.array[model.Index,1:5]) , method="Richardson",    Y=Y,n=n,branchlength=branchlength,ancdes.array=ancdes.array,nleaves=nleaves,tipnames=tipnames,model.Index=model.Index)
    para.cov.matrix[model.Index,,]<-solve(Hessian.mtx[model.Index,,])
  }
  
  mu.para.vect<-output.array[,1]
  mu.var.para.vect<-c(abs(para.cov.matrix[,1,1]))
  mu.weight.vect<-output.array[,9]
  mu.weight<-sum(mu.para.vect*mu.weight.vect)
  mu.se<-se.ave.weight.para(mu.para.vect,mu.var.para.vect,mu.weight.vect)
  
  beta.para.vect<-output.array[,2]
  beta.var.para.vect<-c(abs(para.cov.matrix[,2,2]))
  beta.weight.vect<-output.array[,9]
  beta.weight<-sum(beta.para.vect*beta.weight.vect)
  beta.se<-se.ave.weight.para(beta.para.vect,beta.var.para.vect,beta.weight.vect)
  
  h.var.para.vect<-c(abs(para.cov.matrix[,3,3]))
  
  sigma_sq.para.vect<-output.array[,4]
  sigma_sq.var.para.vect<-c(abs(para.cov.matrix[,4,4]))
  sigma_sq.weight.vect<-output.array[,9]
  sigma_sq.weight<-sum(sigma_sq.para.vect*sigma_sq.weight.vect)
  sigma_sq.se<-se.ave.weight.para(sigma_sq.para.vect,sigma_sq.var.para.vect,sigma_sq.weight.vect)
  
  
  sigma.H_sq.para.vect<-output.array[,5]
  sigma.H_sq.var.para.vect<-c(abs(para.cov.matrix[,5,5]) )
  sigma.H_sq.weight.vect<-output.array[,9]
  sigma.H_sq.weight<-sum(sigma.H_sq.para.vect*sigma.H_sq.weight.vect)
  sigma.H_sq.se<-se.ave.weight.para(sigma.H_sq.para.vect,sigma.H_sq.var.para.vect,sigma.H_sq.weight.vect)
  
  
  
  output.array<-output.array[order(output.array[,8]),]
  output.array<-cbind(output.array, matrix(mu.var.para.vect,ncol=1),
                      matrix(beta.var.para.vect,ncol=1),
                      matrix(h.var.para.vect,ncol=1),
                      matrix(sigma_sq.var.para.vect,ncol=1),
                      matrix(sigma.H_sq.var.para.vect,ncol=1)
  )
  
  
  col.names<-c("mu" ,"bt","h","sigma_sq","v.H","negloglike","AICc","Delta.AICc","w","cum(w)", "var(mu)","var(bt)","var(h)","var(sigma_sq)","var(v.H)")
  colnames(output.array)<-col.names
  
  summary.weight<-array(0,c(2,15))
  rownames(summary.weight)<-c("weighted average","s.e")
  summary.weight[1,1]<-mu.weight;summary.weight[2,1]<-mu.se;
  summary.weight[1,2]<-beta.weight;summary.weight[2,2]<-beta.se;
  summary.weight[1,4]<-sigma_sq.weight;summary.weight[2,4]<-sigma_sq.se;
  summary.weight[1,5]<-sigma.H_sq.weight;summary.weight[2,5]<-sigma.H_sq.se;
  
  summary.weight[,6:15]<-NA
  
  output.array<-rbind(output.array,summary.weight)
  cat("Finished.","\n\n")
  cat("Fix parameter h = 0.5","\n\n")
  #print(output.array)
  #print(summary.weight)
  obj$param.est<-round(output.array[1:4,c(1:2,4:5)],2)
  obj$param.se<-round(output.array[1:4,c(11:12,14:15)],2)
  obj$loglik<-round(-output.array[1:4,6],2)
  obj$AICc<-round(output.array[1:4,7],2)
  obj$Akaike.weight<-round(output.array[1:4,9],2)
  obj$wgt.param.est.se<-round(output.array[5:6,c(1:2,4:5)],2)
  names(convergence.record)<- c("bt=1","v.H=0","bt=1;v.H=0","free")
  obj$convergence<- convergence.record
  print("0: estimation is convergent")
  #library(ade4)
  #plot.phylog(newick2phylog(x.tre, FALSE))
  #return(output.array)
  return(obj)
}
