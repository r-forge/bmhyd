library(corpcor)
library(optimx)
library(numDeriv)

newick2phylog<-function (x.tre, add.tools = FALSE, call = match.call()) 
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
    if (!add.tools) 
        return(res)
    return(newick2phylog.addtools(res))
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
    


ord_fcn<-function(){ #this function gets the acenstor-descendants relationship.
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
  
  
  #if(FALSE){  
  #The following delete the unnecessarily taxa in Fish ancdes.array set.
srd<- c()

for(i in 1:nleaves){
   srd<-c(srd,paste("X",i,sep=""))
   }
#print(srd)
   rnmt<-rownames(mtx)
#  print(rnmt)
   newmtx<-array(0,c(dim(mtx)))
   rownames(newmtx)<-srd
   colnames(newmtx)<-srd
   for( i in 1: nleaves){
     for(j in 1:nleaves){
       newmtx[i,j]= mtx[ which(rnmt==srd[i]), which(rnmt==srd[j]) ]
       }
     }  
  mtx<-newmtx
#print(mtx)
  #       }#end of if FALSE
  
    
    mtx<-mtx*sigma_sq   

    hybrid.Index<-hybrid.node(ancdes.array) #extra burst for hybrid
#print(ancdes.array)	
#	print(hybrid.Index)
   #print(mtx)
   #print("no problem here")
    for(i in hybrid.Index){
          # print(i)
 	   mtx[i,i]<-mtx[i,i]+sigma.H_sq
	   #print(mtx[i,i])
     }

    #print(mtx)
   if(model.Index==1){rm(bt)}
   if(model.Index==2){rm(sigma.H_sq)}
   if(model.Index==3){rm(bt);rm(sigma.H_sq)}
#print(mtx)
    return(mtx)
    }#end of cov_mtx 


NegLogLike<-function(x,Y=Y,n=n,branchlength=branchlength,ancdes.array=ancdes.array,nleaves=nleaves,tipnames=tipnames,model.Index=model.Index){
  badval<-(0.5)*.Machine$double.xmax
    mu<-x[1]
    bt<-x[2]
    h<-x[3]
    sigma_sq<-x[4]
    sigma.H_sq<-x[5]
    
     #FOR DEBUGGING
    # mu<-0
    #bt<-1
      h<-1/2
    #sigma_sq<-1
    #sigma.H_sq<-0 
    if(model.Index==1){bt<-1}
    if(model.Index==2){sigma.H_sq<-0}
    if(model.Index==3){bt<-1;sigma.H_sq<-0}

    #print(h)
    W <- cov_mtx(c(bt,h,sigma_sq,sigma.H_sq),branchlength=branchlength,ancdes.array=ancdes.array,nleaves=nleaves,tipnames=tipnames,model.Index=model.Index) #NOTE TO BCO: CHECK THAT t IS BEING USED CORRECTLY
   #which(W==)
     
   #print(W)

    hybrid.Index<-hybrid.node(ancdes.array) #extra burst for hybrid
    #print(hybrid.Index)
    muone<- mu*matrix(1,nrow=n)
    for(i in hybrid.Index){
 	   muone[i]<-bt*muone[i] 
	   }   
    #debugging things by Brian
    #print("Doing value")
    #print(x)
    #print(cbind(Y,muone))
    #print(W)
#    print(solve(W))
   #print(det(W))
    NegLogML <- n/2*log(2*pi)+1/2*t(Y-muone)%*%pseudoinverse(W)%*%(Y-muone) + 1/2*log(abs(det(W))) 
    #print(logML)   
  if(min(W)<0 || h<0 || h>1 || sigma_sq <0 || sigma.H_sq<0 ||  bt <= 0.0000001) {
      NegLogML<-badval 
    }
    #print(NegLogML[1])
    if(model.Index==1){rm(bt)}
    if(model.Index==2){rm(sigma.H_sq)}
    if(model.Index==3){rm(bt);rm(sigma.H_sq)}
   
 return(NegLogML[1]) #need to put this in to get scalar output
    }#end of NegLogLike.
    



### The cov_mtx function will return the covaraince matrix. 
Hessian.cov_mtx<-function(x,branchlength=branchlength,ancdes.array=ancdes.array,nleaves=nleaves,tipnames=tipnames,model.Index=model.Index){    #now it is a function of bt,  h, sigma^2 and sigma_H^2
    bt<-x[1]
    h<-x[2]
    sigma_sq<-x[3]
    sigma.H_sq<-x[4]
    
    
    #bt<-1
    #h<-1/2
    #sigma_sq<-1
    #sigma.H_sq<-0 
    #if(model.Index==1){bt<-1}
    #if(model.Index==2){sigma.H_sq<-0}
    #if(model.Index==3){bt<-1;sigma.H_sq<-0}

    
    
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
        
       #####################
        vhii<- bt^2*h^2*omtx[prn1,prn1]+bt^2*(1-h)^2*omtx[prn2,prn2]+2*bt^2*h*(1-h)*omtx[prn1,prn2] 
        #print(vhii)
        #fill in value into matrix with formula var(R)=  
        ######################        
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
        #print(tmp)
        mtx<-diag(branchlength[1],c(length(tmp)))
        rownames(mtx)<-c(tmp)
        colnames(mtx)<-c(tmp)
        omtx<-mtx
        #print(mtx)
        #print(eigen(mtx)$values)
        }#end of loop==2 
    while(loop<length(bmtp)){#loaded the acenstor-descendant ancdes.array. 
        loop<-loop+1#use loop to use the ancdes.array
        #print(loop)
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
                
                #print(tmp)
                #print(eigen(mtx)$values)
                #print(mtx)
                #above generate sequence and cov. matrix for speciation.           

           }else{#  op = 0 means that we have detected the hybridize event.
                ins<-(bmtp[q])#grab the insertion element, ins will be used in the fillhydmtx function.
                insp<-ins_fcn(ins,otmp)#catch the position for insertion.
                tmp[insp]<-ins #insert the hyd element.
                tmp[(insp+1):length(tmp)]=otmp[insp:(length(tmp)-1)]
                tmp[1:(insp-1)]=otmp[1:(insp-1)]
                #print(tmp)
                rownames(mtx)<-c(tmp)
                colnames(mtx)<-c(tmp)
                diamtx<-var_fcn()
                mtx<- fillhydmtx()
                mtx<-mtx+t(mtx)+diamtx
                #fill in the diagonal elements.
               
                otmp<-tmp
                omtx<-mtx
                #print(tmp)
                #print(eigen(mtx)$values)
                #print(mtx)
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
             #print(tmp) 
             #print(eigen(mtx)$values)  
             #print(mtx)
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
            #print(tmp)
            #print(eigen(mtx)$values)  
            #print(mtx)
            }#end of if (sum(tipnames%in%tmp)!=nleaves).
  
  
  #if(FALSE){  
  #The following delete the unnecessarily taxa in Fish ancdes.array set.
   

srd<- c()

for(i in 1:nleaves){
   srd<-c(srd,paste("X",i,sep=""))
   }
#print(srd)
   rnmt<-rownames(mtx)
   newmtx<-array(0,c(dim(mtx)))
   rownames(newmtx)<-srd
   colnames(newmtx)<-srd
   for( i in 1: nleaves){
     for(j in 1:nleaves){
       newmtx[i,j]= mtx[ which(rnmt==srd[i]), which(rnmt==srd[j]) ]
       }
     }  
  mtx<-newmtx
     #       }#end of if FALSE
  
    
    mtx<-mtx*sigma_sq   

    hybrid.Index<-hybrid.node(ancdes.array) #extra burst for hybrid
#print(hybrid.Index)
    for(i in hybrid.Index){
 	   mtx[i,i]<-mtx[i,i]+sigma.H_sq
	   #print(mtx[i,i])
     }

    #print(mtx)
   #if(model.Index==1){rm(bt)}
   #if(model.Index==2){rm(sigma.H_sq)}
   #if(model.Index==3){rm(bt);rm(sigma.H_sq)}
   
   return(mtx)
    }#end of cov_mtx 


hybrid.node<-function(ancdes.array){
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
    
     #FOR DEBUGGING
    # mu<-0
    #bt<-1
    #h<-1/2
    #sigma_sq<-1
    #sigma.H_sq<-0 
    #if(model.Index==1){bt<-1}
    #if(model.Index==2){sigma.H_sq<-0}
    #if(model.Index==3){bt<-1;sigma.H_sq<-0}


    #print(x)
    W <- Hessian.cov_mtx(c(bt,h,sigma_sq,sigma.H_sq),branchlength=branchlength,ancdes.array=ancdes.array,nleaves=nleaves,tipnames=tipnames,model.Index=model.Index) #NOTE TO BCO: CHECK THAT t IS BEING USED CORRECTLY
   #which(W==)
     
   #print(W)

    hybrid.Index<-hybrid.node(ancdes.array) #extra burst for hybrid
    #print(hybrid.Index)
    muone<- mu*matrix(1,nrow=n)
    for(i in hybrid.Index){
 	   muone[i]<-bt*muone[i] 
	   }   
    #debugging things by Brian
    #print("Doing value")
    #print(x)
    #print(cbind(Y,muone))
    #print(W)
#    print(solve(W))
   #print(det(W))
    NegLogML <- n/2*log(2*pi)+1/2*t(Y-muone)%*%pseudoinverse(W)%*%(Y-muone) + 1/2*log(abs(det(W))) 
    #print(logML)   
  #if(min(W)<0 || h<0 || h>1 || sigma_sq <0 || sigma.H_sq<0 ||  bt < 1) {
    #  NegLogML<-badval 
    #}
    #print(NegLogML[1])
    #if(model.Index==1){rm(bt)}
    #if(model.Index==2){rm(sigma.H_sq)}
    #if(model.Index==3){rm(bt);rm(sigma.H_sq)}
   
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

Para.Var.Matrix<-
function(k){
if(k==1){return(cov.set.bt.1.s.H.0)}
if(k==2){return(cov.set.bt.1)}
if(k==3){return(cov.set.s.H.0)}
if(k==4){return(cov.set.free)}
}

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

bmhyd<-function(Y,x.tre){
res<<-newick2phylog(x.tre)
#print(res)
ancdes.array<<-getntn(res)
#print(ancdes.array)
branchlength<<-getbrnlen(res)
#print(branchlength)
tipnames<<-sort(names(res$droot[which(res$droot==1)]))
#print(tipnames)
nleaves<<-length(tipnames)
n<<-nleaves 
#print(n)
output.array<-array(0,c(4,7))
rownames(output.array)<-c("bt=1","v.H=0","bt=1;v.H=0","free")
opt.method<-c("Nelder-Mead")
#print(output.array)
Hessian.mtx<-array(0,c(4,5,5))
para.cov.matrix<-array(0,c(4,5,5))
p0 = c(mean(Y),1,1/2,sd(Y),sd(Y))#starting point
cat("Begin optimx optimization routine -- 
    Starting value: (mu,beta,h,sigma^2,v.H)=(",round(p0,digit=2),")","\n\n")
    convergence.record<-array(0,dim(output.array)[1])
for (model.Index in 1:dim(output.array)[1]){
    if(model.Index==1){k<-3;bt<-1}
    if(model.Index==2){k<-3;sigma.H_sq<-0}
    if(model.Index==3){k<-2;bt<-1;sigma.H_sq<-0}
    if(model.Index==4){k<-4}
    #print(model.Index)
    #print(p0)
    #print(NegLogLike(p0,Y=Y,n=n,branchlength=branchlength,ancdes.array=ancdes.array,nleaves=nleaves,tipnames=tipnames,model.Index=model.Index ))
    #print("Go to optimx")
    #MLE.ALL<-optimx(p0,NegLogLike,method="Nelder-Mead", Y=Y,n=n,branchlength=branchlength,ancdes.array=ancdes.array,nleaves=nleaves,tipnames=tipnames,model.Index=model.Index )
    MLE.ALL<-optim(p0,NegLogLike,method="Nelder-Mead", Y=Y,n=n,branchlength=branchlength,ancdes.array=ancdes.array,nleaves=nleaves,tipnames=tipnames,model.Index=model.Index )
   #  MLE.ALL<-optimx(p0,NegLogLike,all.methods=TRUE, Y=Y,n=n,branchlength=branchlength,ancdes.array=ancdes.array,nleaves=nleaves,tipnames=tipnames,model.Index=model.Index )

    
    print(MLE.ALL)
    convergence.record[model.Index]<- MLE.ALL$convergence
    if(MLE.ALL$convergence==0){cat("The MLE estimation converge")}
     #for(j in 1:5){
    #          output.array[model.Index,j]<- as.numeric(unlist(MLE.ALL[1,1])[j])
    #          }#end of j
    output.array[model.Index,1:5]<-MLE.ALL$par
    #output.array[model.Index,6]<-as.numeric(unlist(MLE.ALL[1,2])[1])
    output.array[model.Index,6]<-MLE.ALL$value
    #print(c(n,k, output.array[model.Index,6]))
    output.array[model.Index,7]<-AICc(n,k,MLE.ALL$value)
               #print( output.array[model.Index,7])


    if(model.Index==1){rm(k);rm(bt)}
    if(model.Index==2){rm(k);rm(sigma.H_sq)}
    if(model.Index==3){rm(k);rm(bt);rm(sigma.H_sq)}

    }#end of for model.Index

obj<-NULL
#print(output.array)
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


col.names<-c("mu" ,"bt","h","sigma_sq","v.H","negloglike","AICc","Delta.AICc","w","cum(w)", "se(mu)","se(bt)","se(h)","se(sigma_sq)","se(v.H)")
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
obj$loglik<-round(output.array[1:4,6],2)
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

	
	
	
#5 taxa 1 HYD event EX 1
#x.tre<-c("(((1:0.5,(2:0.5)7:0)6:0.1,(7:0,3:0.5)8:0.1)10:0.4,(4:0.3,5:0.3)9:0.7)11:0;")
#Y<-matrix(c(1.2,3.5,1.6,2.7,4.3),ncol=1)
	
	
#5 taxa  2 HYD events nongall  Example in Tex
#x.tre=c("((1:0.3,(2:0.3)10:0)9:0.7,((10:0,(3:0.1,(4:0.1)7:0)6:0.2)11:0.3,(7:0,5:0.1)8:0.5)12:0.4)13:0; ")
#Y<-rnorm(5)
	
#6 taxa 2 HYD 1 spc:  Algorithm example
#x.tre<-c( "((1:0.65,((2:0.3,3:0.3)7:0.35)12:0)11:0.35,((12:0,(4:0.5,(5:0.5)9:0)8:0.15)13:0.25,(9:0,6:0.5)10:0.4)14:0.1)15:0;")
#Y<-rnorm(6)
	
	
#8taxa 2HYD event Cardona et al 208  Paper ancdes.array EX 4
#x.tre=c("((1:0.85,((2:0.35,(3:0.15,(4:0.15)10:0)9:0.2)12:0.25,(((10:0,5:0.15)11:0.25,6:0.4)13:0.2)15:0)14:0.25)17:0.15,((15:0,7:0.6)16:0.15,8:0.75)18:0.25)19:0;")
#Y<-rnorm(8)
	
	
# 9 taxa ancdes.array

#x.tre=c("(((1:0.35,2:0.35)10:0.5,(3:0.7,(4:0.45,(5:0.1,(6:0.1)13:0)12:0.35)11:0.25)16:0.15)18:0.15,((13:0,7:0.1)14:0.55,(8:0.2,9:0.2)15:0.45)17:0.35)19:0;")
#Y<-matrix(c( 16.00,12.00,16.00,23.00,18.00,18.00,25.00,20.00,14.00  ),ncol=1)

#result<-bmhyd(Y,x.tre)
#print(result)



if(FALSE){
#-------------bt vs v.H



NegLogLike.contour<-function(x,Y=Y,n=n,branchlength=branchlength,ancdes.array=ancdes.array,nleaves=nleaves,tipnames=tipnames,fixed.param.value=fixed.param.value,assigned.param.Index=assigned.param.Index){
  badval<-(0.5)*.Machine$double.xmax
    mu<-x[1]
    bt<-x[2]
    h<-x[3]
    sigma_sq<-x[4]
    sigma.H_sq<-x[5]
    
    if(assigned.param.Index==1){
    bt<-fixed.param.value[1]
    sigma.H_sq<-fixed.param.value[2]}
  
   if(assigned.param.Index==2){
    h<-fixed.param.value[1]
    sigma.H_sq<-fixed.param.value[2]}
  
   if(assigned.param.Index==3){
    h<-fixed.param.value[1]
    bt<-fixed.param.value[2]}
  
    #print(c(bt,h,sigma_sq,sigma.H_sq))
  
    #print(bt)
  
    #FOR DEBUGGING
    # mu<-0
    #bt<-1
    # h<-1/2
    #sigma_sq<-1
    #sigma.H_sq<-0 
    #if(model.Index==1){bt<-1}
    #if(model.Index==2){sigma.H_sq<-0}
    #if(model.Index==3){bt<-1;sigma.H_sq<-0}

    #print(h)
    W <- cov_mtx.contour(c(bt,h,sigma_sq,sigma.H_sq),branchlength=branchlength,ancdes.array=ancdes.array,nleaves=nleaves,tipnames=tipnames,fixed.param.value=fixed.param.value,assigned.param.Index=assigned.param.Index) #NOTE TO BCO: CHECK THAT t IS BEING USED CORRECTLY
   #which(W==)
     
   #print(W)

    hybrid.Index<-hybrid.node(ancdes.array) #extra burst for hybrid
    #print(hybrid.Index)
    muone<- mu*matrix(1,nrow=n)
    for(i in hybrid.Index){
      muone[i]<-bt*muone[i] 
     }   
    #debugging things by Brian
    #print("Doing value")
    #print(x)
    #print(cbind(Y,muone))
    #print(W)
#    print(solve(W))
   #print(det(W))
    NegLogML <- n/2*log(2*pi)+1/2*t(Y-muone)%*%pseudoinverse(W)%*%(Y-muone) + 1/2*log(abs(det(W))) 
    #print(logML)   
  if(min(W)<0 || h<0 || h>1 || sigma_sq <0 || sigma.H_sq<0 ||  bt <= 0.0000001) {
      NegLogML<-badval 
    }
    #print(NegLogML[1])
    #if(model.Index==1){rm(bt)}
    #if(model.Index==2){rm(sigma.H_sq)}
    #if(model.Index==3){rm(bt);rm(sigma.H_sq)}
   
 return(NegLogML[1]) #need to put this in to get scalar output
    }#end of NegLogLike.





cov_mtx.contour<-function(x,branchlength=branchlength,ancdes.array=ancdes.array,nleaves=nleaves,tipnames=tipnames,fixed.param.value=fixed.param.value,assigned.param.Index=assigned.param.Index){    #now it is a function of bt,  h, sigma^2 and sigma_H^2
    bt<-x[1]
    h<-x[2]
    sigma_sq<-x[3]
    sigma.H_sq<-x[4]
        if(assigned.param.Index==1){
    bt<-fixed.param.value[1]
    sigma.H_sq<-fixed.param.value[2]}
  
   if(assigned.param.Index==2){
    h<-fixed.param.value[1]
    sigma.H_sq<-fixed.param.value[2]}
  
   if(assigned.param.Index==3){
    h<-fixed.param.value[1]
    bt<-fixed.param.value[2]}
    #print(c(bt,h,sigma_sq,sigma.H_sq))
    #bt<-1
     #h<-1/2
    #sigma_sq<-1
    #sigma.H_sq<-0 
   # if(model.Index==1){bt<-1}
  #  if(model.Index==2){sigma.H_sq<-0}
  #  if(model.Index==3){bt<-1;sigma.H_sq<-0}

    #print(h)
    
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
        
       #####################
        vhii<- bt^2*h^2*omtx[prn1,prn1]+bt^2*(1-h)^2*omtx[prn2,prn2]+2*bt^2*h*(1-h)*omtx[prn1,prn2] 
        #print(vhii)
        #fill in value into matrix with formula var(R)=  
        ######################        
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
    #print(rept)    


    loop<-2
    tmp=array(0,c(loop))
    if(loop==2){tmp=rept[1,]
        otmp<-tmp
        #print(tmp)
        mtx<-diag(branchlength[1],c(length(tmp)))
        rownames(mtx)<-c(tmp)
        colnames(mtx)<-c(tmp)
        omtx<-mtx
        #print(mtx)
        #print(eigen(mtx)$values)
        }#end of loop==2 
    while(loop<length(bmtp)){#loaded the acenstor-descendant ancdes.array. 
        loop<-loop+1#use loop to use the ancdes.array
        #print(loop)
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
                
#print(tmp)
                #print(eigen(mtx)$values)
                #print(mtx)
                #above generate sequence and cov. matrix for speciation.           

           }else{#  op = 0 means that we have detected the hybridize event.
                ins<-(bmtp[q])#grab the insertion element, ins will be used in the fillhydmtx function.
                insp<-ins_fcn(ins,otmp)#catch the position for insertion.
                tmp[insp]<-ins #insert the hyd element.
                tmp[(insp+1):length(tmp)]=otmp[insp:(length(tmp)-1)]
                tmp[1:(insp-1)]=otmp[1:(insp-1)]
                #print(tmp)
                rownames(mtx)<-c(tmp)
                colnames(mtx)<-c(tmp)
                diamtx<-var_fcn()
                mtx<- fillhydmtx()
                mtx<-mtx+t(mtx)+diamtx
                #fill in the diagonal elements.
               
                otmp<-tmp
                omtx<-mtx
#print(tmp)
                #print(eigen(mtx)$values)
                #print(mtx)
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
             #print(tmp) 
             #print(eigen(mtx)$values)  
             #print(mtx)
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
            #print(tmp)
            #print(eigen(mtx)$values)  
            #print(mtx)
            }#end of if (sum(tipnames%in%tmp)!=nleaves).
  
  
  #if(FALSE){  
  #The following delete the unnecessarily taxa in Fish ancdes.array set.
srd<- c()

for(i in 1:nleaves){
   srd<-c(srd,paste("X",i,sep=""))
   }
#print(srd)
   rnmt<-rownames(mtx)
#  print(rnmt)
   newmtx<-array(0,c(dim(mtx)))
   rownames(newmtx)<-srd
   colnames(newmtx)<-srd
   for( i in 1: nleaves){
     for(j in 1:nleaves){
       newmtx[i,j]= mtx[ which(rnmt==srd[i]), which(rnmt==srd[j]) ]
       }
     }  
  mtx<-newmtx
#print(mtx)
  #       }#end of if FALSE
  
    
    mtx<-mtx*sigma_sq   

    hybrid.Index<-hybrid.node(ancdes.array) #extra burst for hybrid
#print(ancdes.array)  
#  print(hybrid.Index)
   #print(mtx)
   #print("no problem here")
    for(i in hybrid.Index){
          # print(i)
 	   mtx[i,i]<-mtx[i,i]+sigma.H_sq
	   #print(mtx[i,i])
     }

    #print(mtx)
  # if(model.Index==1){rm(bt)}
  # if(model.Index==2){rm(sigma.H_sq)}
  # if(model.Index==3){rm(bt);rm(sigma.H_sq)}
#print(mtx)
    return(mtx)
    }#end of cov_mtx 

#par(mfrow=c(3,1))


# assigned.param==1)bt, sigma.H_sq
# assigned.param==2)h,  sigma.H_sq
# assigned.param==3)h,  bt

assigned.param.Index<-1

sigma.H_sq.array<-seq(0,4.8,4.8/10)
bt.array<-seq(1.71-2*0.42,1.71+2*0.42,(2*0.42)/10)  

NegLoglike.bt.v.H.array<-array(0,c(length(bt.array),length(sigma.H_sq.array)))
NegLoglike.bt.v.H.array.mle<-array(0,c(length(bt.array),length(sigma.H_sq.array),5))

for(bt.arrayIndex in 1:length(bt.array)){
    bt<-bt.array[bt.arrayIndex]
   for(v.H.arrayIndex in 1:length(sigma.H_sq.array)){
     sigma.H_sq<-sigma.H_sq.array[v.H.arrayIndex]
     fixed.param.value<-c(bt,sigma.H_sq)
     print(c(bt,sigma.H_sq))
     sim.result<-optim(p0,NegLogLike.contour,method="Nelder-Mead", Y=Y,n=n,branchlength=branchlength,ancdes.array=ancdes.array,nleaves=nleaves,tipnames=tipnames,fixed.param.value=fixed.param.value,assigned.param.Index=assigned.param.Index)
     print(sim.result$value)
     NegLoglike.bt.v.H.array.mle[bt.arrayIndex,v.H.arrayIndex,]<-sim.result$par
     NegLoglike.bt.v.H.array[bt.arrayIndex,v.H.arrayIndex]<- sim.result$value
     #MLE.ALL<-optimx(p0,NegLogLike,all.methods=TRUE, Y=Y,n=n,branchlength=branchlength,ancdes.array=ancdes.array,nleaves=nleaves,tipnames=tipnames,model.Index=model.Index ) 
    }
  }

NegLoglike.bt.v.H.array<-NegLoglike.bt.v.H.array-min(NegLoglike.bt.v.H.array)
#pdf("contour.grey.bt.v.H.pdf",width=10,height=10)

#contour(bt.array,sigma.H_sq.array,NegLogLike.array,xlab=expression(beta),ylab=expression(v[H]))

# plain
#contour(x,y,z) 
# adjusting levels
for(rowIndex in 1:dim(NegLoglike.bt.v.H.array)[1]){
  for(colIndex in 1:dim(NegLoglike.bt.v.H.array)[2]){
    if( NegLoglike.bt.v.H.array[rowIndex,colIndex]==min(NegLoglike.bt.v.H.array)){
      min.Index<-c(rowIndex,colIndex )
      print( c(rowIndex,colIndex ))
      }  
    }
  }


mylevels <- round(seq(0,max(NegLoglike.bt.v.H.array), max(NegLoglike.bt.v.H.array)/20),4)                               

contour(bt.array,sigma.H_sq.array,NegLoglike.bt.v.H.array,levels=mylevels,xlab=expression(beta),ylab=expression(v[H]))

points(bt.array[min.Index[1]],sigma.H_sq.array[min.Index[2]],pch=8)




points(bt.array[min.Index[1]],sigma.H_sq.array[min.Index[2]],pch=8)


# assigned.param==1)bt, sigma.H_sq
# assigned.param==2)h,  sigma.H_sq
# assigned.param==3)h,  bt

p0<-c(mean(Y),1,0.5,sd(Y),0)

sigma.H_sq.array<-seq(0,2*4.8,4.8/10)
h.array<-seq(0.01,1,0.1)

assigned.param.Index<-2
NegLoglike.h.v.H.array<-array(0,c(length(h.array),length(sigma.H_sq.array)))
NegLoglike.h.v.H.array.mle<-array(0,c(length(h.array),length(sigma.H_sq.array),5))

for(h.arrayIndex in 1:length(h.array)){
    h<-h.array[h.arrayIndex]
   for(v.H.arrayIndex in 1:length(sigma.H_sq.array)){
     sigma.H_sq<-sigma.H_sq.array[v.H.arrayIndex]
     fixed.param.value<-c(h,sigma.H_sq)
     print(c(h,sigma.H_sq))
     sim.result<-optim(p0,NegLogLike.contour,method="Nelder-Mead", Y=Y,n=n,branchlength=branchlength,ancdes.array=ancdes.array,nleaves=nleaves,tipnames=tipnames,fixed.param.value=fixed.param.value,assigned.param.Index=assigned.param.Index)
     print(sim.result$value)
     NegLoglike.h.v.H.array.mle[h.arrayIndex,v.H.arrayIndex,]<-sim.result$par
     NegLoglike.h.v.H.array[h.arrayIndex,v.H.arrayIndex]<- sim.result$value
     #MLE.ALL<-optimx(p0,NegLogLike,all.methods=TRUE, Y=Y,n=n,branchlength=branchlength,ancdes.array=ancdes.array,nleaves=nleaves,tipnames=tipnames,model.Index=model.Index ) 
    }
  }

NegLoglike.h.v.H.array<-NegLoglike.h.v.H.array-min(NegLoglike.h.v.H.array)
#pdf("contour.grey.bt.v.H.pdf",width=10,height=10)

#contour(bt.array,sigma.H_sq.array,NegLogLike.array,xlab=expression(beta),ylab=expression(v[H]))

# plain
#contour(x,y,z) 
# adjusting levels
for(rowIndex in 1:dim(NegLoglike.h.v.H.array)[1]){
  for(colIndex in 1:dim(NegLoglike.h.v.H.array)[2]){
    if( NegLoglike.h.v.H.array[rowIndex,colIndex]==min(NegLoglike.h.v.H.array)){
      min.Index<-c(rowIndex,colIndex )
      print( c(rowIndex,colIndex ))
      }  
    }
  }


mylevels <- round(seq(0,max(NegLoglike.h.v.H.array), max(NegLoglike.h.v.H.array)/20),4)                               

contour(h.array,sigma.H_sq.array,NegLoglike.h.v.H.array,levels=mylevels,xlab=expression(h),ylab=expression(v[H]))

points(h.array[min.Index[1]],sigma.H_sq.array[min.Index[2]],pch=8)





# assigned.param==1)bt, sigma.H_sq
# assigned.param==2)h,  sigma.H_sq
# assigned.param==3)h,  bt

p0<-c(mean(Y),1,0.5,sd(Y),0)

sigma.H_sq.array<-seq(0,10,1)
bt.array<-seq(0.01,5,0.5)  
h.array<-seq(0.01,1,0.1)

assigned.param.Index<-3
NegLoglike.h.bt.array<-array(0,c(length(h.array),length(bt.array)))
NegLoglike.h.bt.array.mle<-array(0,c(length(h.array),length(bt.array),5))

for(h.arrayIndex in 1:length(h.array)){
    h<-h.array[h.arrayIndex]
   for(bt.arrayIndex in 1:length(bt.array)){
     bt<-bt.array[bt.arrayIndex]
     fixed.param.value<-c(h,bt)
     print(c(h,bt))
     sim.result<-optim(p0,NegLogLike.contour,method="Nelder-Mead", Y=Y,n=n,branchlength=branchlength,ancdes.array=ancdes.array,nleaves=nleaves,tipnames=tipnames,fixed.param.value=fixed.param.value,assigned.param.Index=assigned.param.Index)
     print(sim.result$value)
     NegLoglike.h.bt.array.mle[h.arrayIndex,bt.arrayIndex,]<-sim.result$par
     NegLoglike.h.bt.array[h.arrayIndex,bt.arrayIndex]<- sim.result$value
     #MLE.ALL<-optimx(p0,NegLogLike,all.methods=TRUE, Y=Y,n=n,branchlength=branchlength,ancdes.array=ancdes.array,nleaves=nleaves,tipnames=tipnames,model.Index=model.Index ) 
    }
  }

NegLoglike.h.bt.array<-NegLoglike.h.bt.array-min(NegLoglike.h.bt.array)
#pdf("contour.grey.bt.v.H.pdf",width=10,height=10)

#contour(bt.array,sigma.H_sq.array,NegLogLike.array,xlab=expression(beta),ylab=expression(v[H]))

# plain
#contour(x,y,z) 
# adjusting levels
for(rowIndex in 1:dim(NegLoglike.h.bt.array)[1]){
  for(colIndex in 1:dim(NegLoglike.h.bt.array)[2]){
    if( NegLoglike.h.bt.array[rowIndex,colIndex]==min(NegLoglike.h.bt.array)){
      min.Index<-c(rowIndex,colIndex )
      print( c(rowIndex,colIndex ))
      }  
    }
  }


mylevels <- round(seq(0,max(NegLoglike.h.bt.array), max(NegLoglike.h.bt.array)/20),4)                               

contour(h.array,bt.array,NegLoglike.h.bt.array,levels=mylevels,xlab=expression(h),ylab=expression(beta))

points(h.array[min.Index[1]],bt.array[min.Index[2]],pch=8)







#simulate point with radius of se
#sim.ests.set<-array(0,c(dim(result$param.est),100))
#sim.loglik.set<-array(0,c(dim(result$param.est)[1],100))
#est.p.array<-array(0,c(5,dim(result$param.est)[1],100))
#for(model.Index in 1:dim(result$param.est)[1]){
#  for(neigh.Index in 1:100){
#    est.p<-result$param.est +runif(result$param.se)
    #print(est.p)
#    if(model.Index==1){est.p<-c(est.p[1:2],0.5,est.p[3:4])} #bt=1;v.H=0     
#    if(model.Index==2){est.p<-c(est.p[1:2],0.5,est.p[3:4])} #v.H=0       
#    if(model.Index==3){est.p<-c(est.p[1:2],0.5,est.p[3:4])} #bt=1       
#    if(model.Index==4){est.p<-c(est.p[1:2],0.5,est.p[3:4])} #free
    #print(est.p)
#    est.p.array[,model.Index,neigh.Index]<-est.p
#    sim.loglik.set[model.Index,neigh.Index]<- NegLogLike( est.p ,Y=Y,n=n,branchlength=branchlength,ancdes.array=ancdes.array,nleaves=nleaves,tipnames=tipnames,model.Index=model.Index)
#    }
#  }


}#end of if FALSE





