## WORKING
# add a series of fossils to a tree; waiting times drawn according to uniform prior
add_fossils.phylo=function(phy, fossils.df){
	
# fossils.df: data.frame of defining taxa, fossil_age, and name...
#	   t1 t2  time      name
#	1 t42 t53   10 fossil_a
#	2 t33 t34   10 fossil_b
#	3 t31 t55   10 fossil_c
#	4 t48 t38   10 fossil_d
#	5 t54 t55   10 fossil_e
	
	## FIXME: allow for placement of fossil on tip
	for(i in 1:nrow(fossils.df)){
		get.btimes(phy)->times
		tt=match(fossils.df[i,1:2], phy$tip.label)
		if(!length(unique(tt)->aa)==1) aa=getMRCA(phy, tt)
		
	## allow for a segmented stem at base of crown group
		nn=sort(segmented.fossil.stem(aa, phy))
		edge=cumsum(sapply(nn, function(nd) times$edge[which(times$des==nd)]))
		min_anc=times$anc[times$des==min(nn)]
		
	## define min and max for possible split times 
		tip_age=fossils.df$time[i]
		max_age=times$time[times$des==min_anc]
		min_age=max(tip_age, times$time[times$des==max(nn)])
		interval=max_age-min_age
		
	## SPLIT TIME
		split_interval=runif(1,0,interval)
		
		
		split_edge=nn[min(which(split_interval<edge))]
		split_edge_anc=times$anc[times$des==split_edge]
		split_age=max_age-split_interval
		
	## ADD TIP to TREE
		phy=add.fossiltip(phy=phy, stem=split_edge, tip_age=tip_age, split_age=split_age, label=fossils.df$name[i])
	}
	return(phy)
}


## WORKING
# add a single tip to tree given the stem lineage, age of the tip (mya), age at the split (mya), and a taxon label 
# avoids ape:::bind.tree()
add.fossiltip=function(phy, stem, tip_age, split_age, label=""){
	
# stem: numeric edge label
# tip_age: mya for tip to be added
# split_age: mya for where tip connects to stem
# label: character to attach to fossil tip
	
	edge=phy$edge
	lengths=phy$edge.length
	labels=phy$tip.label
	
	times=get.btimes(phy)
	ii=which(times$des==stem)
	stem_min=times$time[ii]
	stem_max=stem_min+times$edge[ii]
	
	if(tip_age>stem_max) stop("fossil 'tip_age' appears inconsistent with supplied 'stem'")
	if(!withinrange(split_age, stem_min, stem_max)) stop("fossil 'split_age' appears inconsistent with timing of supplied 'stem'")
	
# update edge matrix
	ii=which(edge[,2]==stem)
	n=Ntip(phy)
	nn=Nnode(phy)
	if(stem<=n) stop("Cannot yet handle fossil attachment to tips.")
	stem=stem+1
	edge[edge>n]=edge[edge>n]+1
	edge[edge>stem]=edge[edge>stem]+1
	edge[edge[,1]==stem,1]=stem+1
	edge=rbind(edge, c(stem, n+1))
	insert_edge=c(stem, stem+1)	
	
# update edge.length array
	lengths[ii]=stem_max-split_age
	lengths=c(lengths, split_age-tip_age)
	insert_length=split_age-stem_min
	
# update tip.label array
	labels=c(labels, label)
	
# update fossil array
	if(!any(names(phy)=="fossil")) phy$fossil=rep(FALSE,n)
	fossil=c(phy$fossil, TRUE)
	
## finalize phylo object ##
	slot=ii+1
	edge=rbind(edge[1:(slot-1),], insert_edge, edge[slot:nrow(edge),])
	dimnames(edge)=NULL
	lengths=c(lengths[1:(slot-1)], insert_length, lengths[slot:length(lengths)])
	
	new=phy
	new$edge=edge
	new$edge.length=lengths
	new$tip.label=labels
	new$Nnode=nn+1
	new$fossil=fossil
	attr(new, "order")=NULL
	return(reorder(new))
}

## WORKING
# slide internal node attaching a tip according to uniform and given lengths to bounding nodes [max <------------- nd --> min]
# intended to work in MCMC -- movement rootward or tipward is equiprobable
# FIXME: negative brlens in SOME cases ??
push.node=function(phy, tip.nd, rootward.nd, tipward.nd, rootward_interval, tipward_interval){
	
# .nd: a numeric edge ID
# _interval: edge.length for associated edge
	
	rootward.nd_x=phy$edge[,2]==rootward.nd
	tipward.nd_x=phy$edge[,2]==tipward.nd
	tip.nd_x=phy$edge[,2]==tip.nd
	
	if(tipward_interval>0) mv=sample(c(0,1),1) else mv=0
	
	if(mv==0){	## MOVE ROOTWARD
		slice=runif(1,0,rootward_interval)
		phy$edge.length[rootward.nd_x]=phy$edge.length[rootward.nd_x]-slice
		phy$edge.length[tipward.nd_x]=phy$edge.length[tipward.nd_x]+slice
		phy$edge.length[tip.nd_x]=phy$edge.length[tip.nd_x]+slice
	} else {	## MOVE TIPWARD
		slice=runif(1,0,tipward_interval)
		phy$edge.length[rootward.nd_x]=phy$edge.length[rootward.nd_x]+slice
		phy$edge.length[tipward.nd_x]=phy$edge.length[tipward.nd_x]-slice
		phy$edge.length[tip.nd_x]=phy$edge.length[tip.nd_x]-slice
	}
	print(mv)
	
	return(phy)
}

## WORKING
# slide internal node attaching a tip according to uniform [max <-- nd --> min]
push.fossil=function(fossil="", phy){
	
# fossil: character tip.label
	
	nd=which(phy$tip.label==fossil)
	a=get.ancestor.of.node(nd,phy)
	dd=get.desc.of.node(a,phy)
	rootward_interval=phy$edge.length[phy$edge[,2]==a] 
	
# if multifurcation, do not push tipward of ancestral node
	if(length(dd)>2) {
		tipward_interval=0
	} else {
		d=dd[dd!=nd]
		tipward_interval=phy$edge.length[phy$edge[,2]==d] 
	}
	
	return(push.node(phy, nd, a, d, rootward_interval, tipward_interval))
}

## WORKING 
# working rootward from INTERNAL node, find stem(s) associated with a crown group (that may have been segmented by several fossil lineages)
segmented.fossil.stem=function(node, phy){
	
	if(!any(names(phy)=="fossil")) return(node)
	n=Ntip(phy)
	nodes<-nd<-node
	while(1){
		a=get.ancestor.of.node(nd, phy)
		dd=get.desc.of.node(a, phy)
		ff=sapply(dd, function(d) {
				  if(d<=n) {
				  return(phy$fossil[d])
				  } else {
				  return(FALSE)
				  }
				  })
		if(any(ff)) {
			nodes=c(nodes,a)
			nd=a
		} else {
			break()
		}
	}
	return(sort(nodes))
}


## WORKING
# starting from FOSSIL tip, grabs all segments of a formerly continuous stem (broken up by fossils), moving tipward and rootward
# could be useful for NNI or SPR swapping of fossil lineages along a stem
collect.stems=function(fossil="", phy){
# fossil: character tip label
	nd=which(phy$tip.label==fossil)
	a=get.ancestor.of.node(nd,phy)
	dd=get.desc.of.node(a,phy)
	n=Ntip(phy)
	nn=dd[dd!=nd]
	while(1){
		dd=get.desc.of.node(nn, phy)
		ff=sapply(dd, function(d) {
				  if(d<=n) {
				  return(phy$fossil[d])
				  } else {
				  return(FALSE)
				  }
				  })
		if(any(ff) & !all(ff)) nn=dd[!ff] else break()
	}
	stems=segmented.fossil.stem(nn, phy)
	return(stems)
}



