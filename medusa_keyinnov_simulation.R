
require(phylo)
require(diversitree)
require(geiger)

dropextinct.phylo=function(phy, tol=1e-7){
	phy <- reorder(phy)
    n <- length(phy$tip.label)
    n.node <- phy$Nnode
    xx <- numeric(n + n.node)
    for (i in 1:nrow(phy$edge)) xx[phy$edge[i, 2]] <- xx[phy$edge[i, 1]] + phy$edge.length[i]
	
	paths=xx[1:n]
	depth=max(paths)
	tt=abs(depth-paths)>tol
	if(any(tt)){
		phy=drop.tip(phy, phy$tip.label[tt])
	}
	phy
}

richnesses_and_ages.phylo=function(phy){
	if(!any("clade.tree"%in%class(phy))) stop("supply 'phy' as a cladetree")
	if(!is.ultrametric(phy)) {
		phy=ultrametricize(phy)
	}
	rr=sapply(phy$clades, length)
	ww=phy$edge[,2]<=Ntip(phy)
	aa=phy$edge.length[ww]
	names(aa)=phy$tip.label[phy$edge[ww,2]]
	dat=merge(data.frame(rr), data.frame(aa), by=0)
	rownames(dat)=dat$Row.names
	dat=dat[,-1]
	names(dat)=c("richness","age")
	dat
}

cladeages.phylo=function(phy, nodes, type=c("crown", "stem")){
# nodes: defines clades
	require(phylo)
	type=match.arg(type, choices=c("crown", "stem"))
	tt=heights.phylo(phy)
	N=Ntip(phy)
	y=sapply(nodes, function(x) {
			 if(type=="stem" || x<=N){
			 return(tt[x,"start"])
			 } else {
			 return(tt[x,"end"])
			 }
			 })
	names(y)=nodes
	y
}

.ages_and_richnesses.phylo=function(phy, cldphy, nodes, ...){
	nm=sapply(cldphy$clades, function(x) hash.tip(x, phy$tip.label))
	rr=sapply(cldphy$clades, length)
	names(rr)=nm
	phy=hashes.phylo(phy, phy$tip.label)
	ages=cladeages.phylo(phy, nodes, ...)
	names(ages)=phy$hash[nodes]
	dat=merge(data.frame(rr), data.frame(ages), by=0)
	rownames(dat)=dat$Row.names
	dat=dat[,-1]
	names(dat)=c("richness","age")
	dat
}



tocladetree=function (phy, clades) {
	if (!identical(class(phy), "phylo")) 
	stop("tree must be a plain 'phylo' tree")
	if (!all(names(clades) %in% phy$tip.label)) 
	stop("Unknown clade representatives")
	if (!all(sapply(clades, is.character))) 
	stop("'clades' must be a list of character vectors")
	if (any(duplicated(unlist(clades)))) 
	stop("Duplicated species names")
	phy$clades <- clades
	class(phy) <- c("clade.tree", class(phy))
	phy
}

collapse.phylo=function(phy, nodes){
	desctips=lapply(nodes, function(nd) {
					if(nd<=Ntip(phy)) {
						return(phy$tip.label[nd]) 
					} else {
						tmp=get.descendants.of.node(nd, phy, tips=TRUE)
						return(phy$tip.label[tmp])
					}
				})
	droptips=lapply(desctips, function(x) if(length(x)==1) return(NULL) else return(x[2:length(x)]))
	
	prn=drop.tip(phy, unlist(droptips))
	names(droptips)=prn$tip.label
	prn$tip.label=paste("CLADE",1:Ntip(prn),sep="_")
	
	names(desctips)=prn$tip.label
	cldphy=tocladetree(prn, desctips)
	cldphy
}

pull_nodes=function(nd, phy, p=0.5){
	nds=c()
	if(nd<=Ntip(phy)) {
		nds=c(nds, nd)
	} else {
		dd=get.desc.of.node(nd, phy)
		for(i in 1:length(dd)){
			if(runif(1)<p){
				nds=c(nds, dd[i])
				dd=dd[-i]
			}
		}
		if(length(dd)){
			nds=c(nds, sapply(dd, function(x) pull_nodes(nd=x, phy=phy, p=p)))
		}
	}
	x=unlist(nds)
	return(x[!is.na(x)])
}

split_node=function(nd, phy, p=0.5){
	if(nd>Ntip(phy) & runif(1)<p) return(get.desc.of.node(nd, phy)) else return(NULL)
}

rcladetree.phylo=function(phy, clades, p=0.5, type=c("crown","stem")){
	# p: probability of collapsing a node into a clade (starting from root)
	# clades: upper limit on number of clades to collapse tree into
	type=match.arg(type, c("crown","stem"))
	require(phylo)
	N=Ntip(phy)
	nodes=N+1
	while(length(nodes)<clades){
		if(all(nodes<=N)) break()
		tmp=nodes[nodes>N]
		for(i in 1:length(nodes)){
			xx=split_node(nodes[i], phy, p)
			if(!is.null(xx)){
				nodes=c(nodes[-i], xx)
				if(length(nodes)==clades) break()
			}
		}
	}
	
	cldphy=collapse.phylo(phy, nodes)
	dat=.ages_and_richnesses.phylo(phy, cldphy, nodes, type)
	return(list(phy=cldphy, nodes=nodes, dat=dat))
}

rcladetree.time=function(phy, time=0.5, relative=TRUE, type=c("crown","stem")){
	type=match.arg(type, c("crown","stem"))

	require(phylo)
	root=Ntip(phy)+1
	phy <- reorder(phy)
    n <- length(phy$tip.label)
    n.node <- phy$Nnode
    xx <- numeric(n + n.node)
    for (i in 1:nrow(phy$edge)) xx[phy$edge[i, 2]] <- xx[phy$edge[i, 1]] + phy$edge.length[i]
	
	ee=phy$edge.length[match(1:length(xx), phy$edge[,2])]
	ee[is.na(ee)]=0
	tt=cbind(xx, xx-ee)
	depth=max(tt)
	tt=depth-tt
	if(relative==TRUE) tt=tt/depth
	colnames(tt)=c("to","from")
	yy = cbind(tt, inrange = (abs(sign(time - tt[, "from"]) + sign(time - tt[, "to"]))) != 2)
	
	nodes=which(yy[,"inrange"]==1)
	
	tocladetree=function (phy, clades) {
		if (!identical(class(phy), "phylo")) 
		stop("tree must be a plain 'phylo' tree")
		if (!all(names(clades) %in% phy$tip.label)) 
		stop("Unknown clade representatives")
		if (!all(sapply(clades, is.character))) 
		stop("'clades' must be a list of character vectors")
		if (any(duplicated(unlist(clades)))) 
		stop("Duplicated species names")
		phy$clades <- clades
		class(phy) <- c("clade.tree", class(phy))
		phy
	}
	
	desctips=lapply(nodes, function(nd) {
					if(nd<=Ntip(phy)) {
					return(phy$tip.label[nd]) 
					} else {
					tmp=get.descendants.of.node(nd, phy, tips=TRUE)
					return(phy$tip.label[tmp])
					}
					})
	droptips=lapply(desctips, function(x) if(length(x)==1) return(NULL) else return(x[2:length(x)]))
	
	prn=drop.tip(phy, unlist(droptips))
	names(droptips)=prn$tip.label
	prn$tip.label=paste("CLADE",1:Ntip(prn),sep="_")
	
	names(desctips)=prn$tip.label
	cldphy=tocladetree(prn, desctips)
	dat=.ages_and_richnesses.phylo(phy, cldphy, nodes, type)
	
	return(list(phy=cldphy,nodes=nodes,dat=dat))
}

`set.seed.clock` <-
function(print=F){
	date = date()
 	seed1 = as.numeric(strsplit(substring(date,12,19),":")[[1]])%*%c(1,100,10000)
 	seed <- runif(1, min=0, max=50) * seed1
 	set.seed(seed)
 	if(print) cat("Seed = ", seed, "\n");
 	seed[1,1]
}

`bd_keyinnov` <-
function (b, d, change, effect=2, timeStop=0, taxaStop=0, seed=0, print.seed=FALSE, return.all.extinct=FALSE, keepAppearances=FALSE){
	
#	b=.1; d=0.01; change=0.2; effect=2; timeStop=0; taxaStop=2000; seed=0; print.seed=FALSE; return.all.extinct=FALSE; keepAppearances=FALSE
	if(seed==0) seed=set.seed.clock(print=print.seed);
	
	if(timeStop==0 & taxaStop==0)
	stop("Must have stopping criterion\n");
	
	
	
	while(1) {
#		cat("a")
		edge <- rbind(c(1, 2), c(1, 3)) # this is a starting edge matrix
		edge.length <- rep(NA, 2)
		stem.depth <- numeric(2)
		sRates<-rep(1, 2)
		eRates<-rep(1, 2)
		alive<-rep(TRUE, 2)
		t <- 0 #time at any point in the tree
		next.node<-4
		if(keepAppearances) {
			app<-cbind(c(2, 3), c(0, 0), c(NA, NA))
		}
		branchChanges<-rbind(c(1, 0, 1, 1, 1), c(1, 0, 2, 1, 1))
		colnames(branchChanges)<-c("Lineage", "Time", "Type", "Beginning", "End")
		
############
		repeat{
#			cat(sum(alive), "\n")			
			if(taxaStop) if(sum(alive)>=taxaStop) break;
			if(sum(alive)==0) break;
			rateVector<-c(sRates[alive]*b, eRates[alive]*d, rep(change, sum(alive)))
			dt<-rexp(1, sum(rateVector));
	 	 	t<-t+dt;
	  		if(timeStop) if(t>=timeStop) {
				t<-timeStop;
				break;
	  		}
	  		probVector<-cumsum(rateVector/sum(rateVector))
			
	  		r<-runif(1)
	  		event<-min(which(probVector>r))
			
	  		
         	if(event<=sum(alive)) {#this creates a bifucation in the tree
#				cat("bifurc\n")         		
	        	random_lineage <- event
				e<-matrix(edge[alive,], ncol=2)
				parent<-e[random_lineage,2]
				
				## getting parent rate
				ps<-sRates[alive][random_lineage]
				pe<-sRates[alive][random_lineage]
				sRates<-c(sRates, ps, ps)
				eRates<-c(eRates, pe, pe)
				
				alive[alive][random_lineage]<-FALSE
				edge<-rbind(edge, c(parent, next.node), c(parent, next.node+1))
				if(keepAppearances) {
					app<-rbind(app, c(next.node+1, t, NA))
					app[app[,1]==parent,1]<-next.node
				}
				next.node<-next.node+2
				
				
				alive<-c(alive, TRUE, TRUE)
				stem.depth<-c(stem.depth, t, t)
				x<-which(edge[,2]==parent)
				edge.length[x]<-t-stem.depth[x]
				edge.length<-c(edge.length, NA, NA)
				
				
            } else if(event<=sum(alive)*2) {#This terminates one of the current lineages on the tree
#				cat("death\n")         		
            	
            	random_lineage <- event-sum(alive)
		    	edge.length[alive][random_lineage]<-t-stem.depth[alive][random_lineage];
          	    if(keepAppearances) {
          	    	killed<-edge[alive,2][random_lineage]
					app[app[,1]==killed,3]<-t
				}
				alive[alive][random_lineage]<-FALSE
				
            } else {
#				cat("change\n")         		
            	
            	random_lineage <- event-2*sum(alive)
            	rr<-runif(1)
            	if(rr<1) {  ## FORCE EFFECT ON SPECIATION ONLY
#					cat("c1\n")            		
            		et<-1
            		or<-sRates[alive][random_lineage]
            		while(1) {
            			nr<-effect*or
            			break
#						cat("c2\n")            		
            			
					}
            		sRates[alive][random_lineage]<-nr
				} else {
					et<-2
					or<-eRates[alive][random_lineage]
					while(1) {
						nr<-sample(rateCats)[1]
            			if(nr!=or) break
					}
					eRates[alive][random_lineage]<-sample(rateCats)[1]
					nr<-eRates[alive][random_lineage]
				}
				if(or!=nr)
				branchChanges<-rbind(branchChanges, c(edge[alive,2][random_lineage], t, et, or, nr))
            }
#            cat(sum(alive), "\n")
      	}
		
		if(return.all.extinct==T | sum(alive)>1) break;
	}
	edge.length[alive]<-t-stem.depth[alive]
#	cat("b")	
	
	n<--1;
	for(i in 1:max(edge)) {
		if(any(edge[,1]==i)) {
			edge[which(edge[,1]==i), 1]<-n
			edge[which(edge[,2]==i), 2]<-n
			if(keepAppearances) app[which(app[,1]==i), 1]<-n
			n<-n-1
		}
	}
	
	edge[edge[,2]>0,2]<-1:sum(edge>0)
	
	
	tips<-!(edge[,2] %in% edge[,1])
	tip.label<-edge[tips,2]
	
	if(keepAppearances) {
		mm<-match(app[,1], edge[,2])
		app[,1]<-edge[mm,2]
	}
	
#tip.label<-1:sum(edge>0)
	
	mode(tip.label) <- "character"
	obj <- list(edge = edge, edge.length = edge.length, tip.label=tip.label)
	class(obj) <- "phylo"
	obj<-old2new.phylo(obj)
	obj<-read.tree(text=write.tree(obj))
#	cat("c")    
    if(!keepAppearances) {
		app=NULL
	}
	
	branchChanges=as.data.frame(branchChanges, stringsAsFactors=FALSE)
	branchChanges$Type=c("b","d")[branchChanges$Type]
	
	drp=dropextinct.phylo(obj)
	tt=list(drp=drp, obj=obj)
	class(tt)="multiPhylo"
	hh=hashes.phylo(tt, drp$tip.label)
	branchChanges$hash=hh$obj$hash[hh$obj$edge[branchChanges$Lineage,2]]
	branchChanges$ID=match(branchChanges$hash, hh$drp$hash)
	return(list(phy=hh$drp, app=app, t=t, changes=branchChanges[-c(1:2),]))	
}

`bd_taxonreplacement` <-
function (b, d, change, effect=2, timeStop=0, taxaStop=0, seed=0, print.seed=FALSE, return.all.extinct=FALSE, keepAppearances=FALSE){
	
#	b=.1; d=0.01; change=0.2; effect=2; timeStop=0; taxaStop=2000; seed=0; print.seed=FALSE; return.all.extinct=FALSE; keepAppearances=FALSE
	if(seed==0) seed=set.seed.clock(print=print.seed);
	
	if(timeStop==0 & taxaStop==0)
	stop("Must have stopping criterion\n");
	
	
	
	while(1) {
#		cat("a")
		edge <- rbind(c(1, 2), c(1, 3)) # this is a starting edge matrix
		edge.length <- rep(NA, 2)
		stem.depth <- numeric(2)
		sRates<-rep(1, 2)
		eRates<-rep(1, 2)
		alive<-rep(TRUE, 2)
		t <- 0 #time at any point in the tree
		next.node<-4
		if(keepAppearances) {
			app<-cbind(c(2, 3), c(0, 0), c(NA, NA))
		}
		branchChanges<-rbind(c(1, 0, 1, 1, 1), c(1, 0, 2, 1, 1))
		colnames(branchChanges)<-c("Lineage", "Time", "Type", "Beginning", "End")
		
############
		repeat{
#			cat(sum(alive), "\n")			
			if(taxaStop) if(sum(alive)>=taxaStop) break;
			if(sum(alive)==0) break;
			rateVector<-c(sRates[alive]*b, eRates[alive]*d, rep(change, sum(alive)))
			dt<-rexp(1, sum(rateVector));
	 	 	t<-t+dt;
	  		if(timeStop) if(t>=timeStop) {
				t<-timeStop;
				break;
	  		}
	  		probVector<-cumsum(rateVector/sum(rateVector))
			
	  		r<-runif(1)
	  		event<-min(which(probVector>r))
			
	  		
         	if(event<=sum(alive)) {#this creates a bifucation in the tree
#				cat("bifurc\n")         		
	        	random_lineage <- event
				e<-matrix(edge[alive,], ncol=2)
				parent<-e[random_lineage,2]
				
## getting parent rate
				ps<-sRates[alive][random_lineage]
				pe<-sRates[alive][random_lineage]
				sRates<-c(sRates, ps, ps)
				eRates<-c(eRates, pe, pe)
				
				alive[alive][random_lineage]<-FALSE
				edge<-rbind(edge, c(parent, next.node), c(parent, next.node+1))
				if(keepAppearances) {
					app<-rbind(app, c(next.node+1, t, NA))
					app[app[,1]==parent,1]<-next.node
				}
				next.node<-next.node+2
				
				
				alive<-c(alive, TRUE, TRUE)
				stem.depth<-c(stem.depth, t, t)
				x<-which(edge[,2]==parent)
				edge.length[x]<-t-stem.depth[x]
				edge.length<-c(edge.length, NA, NA)
				
				
            } else if(event<=sum(alive)*2) {#This terminates one of the current lineages on the tree
#				cat("death\n")         		
            	
            	random_lineage <- event-sum(alive)
		    	edge.length[alive][random_lineage]<-t-stem.depth[alive][random_lineage];
          	    if(keepAppearances) {
          	    	killed<-edge[alive,2][random_lineage]
					app[app[,1]==killed,3]<-t
				}
				alive[alive][random_lineage]<-FALSE
				
            } else {
#				cat("change\n")  # increase b in one clade and d in another       		
            	
            	random_lineage <- event-2*sum(alive)
				while(1) {
					r<-runif(1)
					nother<-min(which(probVector>r))
					if(nother>sum(alive)*2){
						another_lineage <- nother-2*sum(alive)
						if(random_lineage!=another_lineage) break
					}
				}
#				cat("c1\n")            		
				et<-1
				or<-sRates[alive][random_lineage]
				nr<-effect*or
				sRates[alive][random_lineage]<-nr
				branchChanges<-rbind(branchChanges, c(edge[alive,2][random_lineage], t, et, or, nr))

				et<-2
				or<-eRates[alive][another_lineage]
				nr<-effect*or
				eRates[alive][another_lineage]<-nr
				branchChanges<-rbind(branchChanges, c(edge[alive,2][another_lineage], t, et, or, nr))
            }
#            cat(sum(alive), "\n")
      	}
		
		if(return.all.extinct==T | sum(alive)>1) break;
	}
	edge.length[alive]<-t-stem.depth[alive]
#	cat("b")	
	
	n<--1;
	for(i in 1:max(edge)) {
		if(any(edge[,1]==i)) {
			edge[which(edge[,1]==i), 1]<-n
			edge[which(edge[,2]==i), 2]<-n
			if(keepAppearances) app[which(app[,1]==i), 1]<-n
			n<-n-1
		}
	}
	
	edge[edge[,2]>0,2]<-1:sum(edge>0)
	
	
	tips<-!(edge[,2] %in% edge[,1])
	tip.label<-edge[tips,2]
	
	if(keepAppearances) {
		mm<-match(app[,1], edge[,2])
		app[,1]<-edge[mm,2]
	}
	
	#tip.label<-1:sum(edge>0)
	
	mode(tip.label) <- "character"
	obj <- list(edge = edge, edge.length = edge.length, tip.label=tip.label)
	class(obj) <- "phylo"
	obj<-old2new.phylo(obj)
	obj<-read.tree(text=write.tree(obj))
#	cat("c")    
    if(!keepAppearances) {
		app=NULL
	}
	
	branchChanges=as.data.frame(branchChanges, stringsAsFactors=FALSE)
	branchChanges$Type=c("b","d")[branchChanges$Type]
	
	drp=dropextinct.phylo(obj)
	tt=list(drp=drp, obj=obj)
	class(tt)="multiPhylo"
	hh=hashes.phylo(tt, drp$tip.label)
	branchChanges$hash=hh$obj$hash[hh$obj$edge[branchChanges$Lineage,2]]
	branchChanges$ID=match(branchChanges$hash, hh$drp$hash)
	return(list(phy=hh$drp, app=app, t=t, changes=branchChanges[-c(1:2),]))	
}

diversities.phylo=function(phy, clades=NULL, time=NULL){
	if(!is.null(clades)) tre=rcladetree.phylo(phy, clades=clades, shifted=c())
	if(!is.null(time))  tre=rcladetree.time(phy, time=time, relative=TRUE)
	dat=richnesses_and_ages.phylo(tre)
	lm=lm(log(dat$richness)~dat$age)
	return(list(N=as.integer(Ntip(tre)), r=coef(lm)[[2]], p=coefficients(summary(lm))[2,4], dat=dat, phy=tre))
}

logscale_density=function(x, ...){
	if(any(x==0)) warning(paste(length(x[x==0])/length(x), "of 'x' is zero."))
	
	
	cvrt=function(x){
		m=max(c(0,abs(log(x[x<1 & x>0]))),na.rm=TRUE)
		f=log(abs(x))+m # scale negative values (between 0 and -1) appropriately
		f[x<0]=-f[x<0]
		oo=order(f)
		return(cbind(x=f[oo], order=oo, orig=x[oo]))
	}
	d=density(x)
	dr=cvrt(d$x)
	plot(dr[,"x"], d$y[dr[,"order"]], type="l", xaxt="n", ...)
	r=pretty(dr[,"orig"])
	nr=sapply(r, function(x) {
			  y=abs(x-dr[,"orig"])
			  dr[which(y==min(y)),"x"]
			  })
	axis(1, at=nr, labels=r)
}


## SIMULATIONS ##
build=FALSE
process=TRUE
require(multicore)
dirn="medusa_AGEDIVERSITY"
if(!file.exists(dirn)) dir.create(dirn)
b=1
d=c(0, 0.25, 0.5, 0.95)
effects=c(0,10)
rates=c(0, 0.25)
reps=1:1000
taxa=1000
parms=expand.grid(rep=reps, b=b, d=d, effect=effects)
parms=cbind(parms, rate=rates[match(parms[,"effect"], effects)])

sets=which(parms[,"rep"]==1)

if(build){
	for(i in 1:length(sets)){
		start=sets[i]
		finish=min(c(sets[i+1]-1, nrow(parms)), na.rm=TRUE)
		filename=paste(dirn, paste(paste("set",i,sep="_"),"rda",sep="."), sep="/")
		res=mclapply(start:finish, function(j){
					 pp=data.frame(parms[j,])
					 sim=bd_keyinnov(b=pp$b, d=pp$d, change=pp$rate, effect=pp$effect, taxaStop=taxa)
					 out=rcladetree.phylo(phy=sim$phy, clades=15, type="crown", p=0.5)
					 f=log(out$dat[,"richness"])~out$dat[,"age"]
					 lm=lm(f)
					 slope=coef(lm)[[2]]
					 
					 sft=unique(sim$changes$ID)
					 pct=sum(out$nodes%in%sft)/length(out$nodes)
					 out$shift_pct=pct
					 out$shifts=length(sft)
					 out$slope=slope
					 out$orig=sim$phy
					 out$changes=sim$changes
					 out$parm=pp
					 return(out)
		})
		save(res, file=filename)
		rm(res)
	}
	
}

if(process){
	out=matrix(NA, nrow=nrow(parms), ncol=4)
	colnames(out)=c("slope", "shifts", "shift_pct", "p_slope")

	files=dir(dirn,full.names=TRUE)
	for(i in 1:length(files)){
		set=as.numeric(basename.noext(gsub("set_","",files[i])))
		start=sets[set]
		finish=min(c(sets[set+1]-1, nrow(parms)), na.rm=TRUE)
		cur=get(load(files[i]))
		for(j in 1:length(cur)){
			if(!inherits(sub<-cur[[j]],"try-error")){
				tmp=.ages_and_richnesses.phylo(sub$orig, sub$phy, sub$nodes, type="crown")
				lm=lm(log(tmp[,"richness"])~tmp[,"age"])
				pv=coefficients(summary(lm))[2,4]
				out[start+(j-1),]=c(sub$slope, sub$shifts, sub$shift_pct, pv)
			}
		}
	}
	
	dat=data.frame(cbind(parms, out))
	save(dat, file="agediversity_progressiveevolution.rda")
	
	layout(matrix(1:(round(length(files)/2)*2), ncol=2))
	for(i in 1:length(files)){
		start=sets[i]
		finish=min(c(sets[i+1]-1, nrow(parms)), na.rm=TRUE)
		bb=unique(dat$b[start:finish])
		dd=unique(dat$d[start:finish])
		rr=unique(dat$rate[start:finish])
		mm=paste(paste(bb,"(b)"), paste(dd, "(d)"), paste(rr, "(rt)"),  sep="; ")
		if(rr==0) {
			xlim=c(-2,2) 
			if(dd==0.95) brk=7 else brk=20
		} else {
			xlim=c(-20,20)
			if(dd==0) brk=100 else brk=200
		}
	
		
		hist(dat$slope[start:finish], breaks=brk, xlim=xlim, main=mm, xlab="slope (age~diversity)", ylab="frequency")
	}
}
