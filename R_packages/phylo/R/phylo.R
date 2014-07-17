compare.phylo=function(phy){
# returns most tipward set of subtrees for each tree in 'phy' that contributes to topological inconsistencies between trees
	if(class(phy)=="phylo") stop("Multiple trees must be supplied")
	class(phy)="multiPhylo"
	hphy=hashes.phylo(phy)
	tmp=table(c(sapply(hphy, function(x) x$hash)))
	common_keys=names(tmp[tmp==length(phy)])
	nds=lapply(hphy, function(x) which(!x$hash%in%common_keys))
	subtrees=lapply(1:length(hphy), function(idx) {
					curtree=hphy[[idx]]
					curnds=nds[[idx]]
					if(!length(curnds)) return(NULL)
					
					# find which subtrees contribute to multiple inconsistencies
					contributes_to=lapply(curnds, function(x) {y=get.ancestors.of.node(x, curtree); return(y[y%in%curnds])})
					drop=sapply(1:length(contributes_to), function(idx) {
								cur=contributes_to[[idx]]
								if(length(cur)>0){
								any(sapply(contributes_to[-idx], function(x) all(cur%in%x)))
								} else {
								return(FALSE)
								}
								})
					problem_nodes=curnds[-c(which(drop))]
					subtrees=lapply(problem_nodes, function(x) ladderize(extract.clade(curtree, x),right=FALSE))
					names(subtrees)=problem_nodes
					return(subtrees)
					})
	names(subtrees)=names(phy)
	return(subtrees)
}

.match.subtrees=function(subtrees, tips){
# subtrees: from compare.phylo
# attempts to find trees consistent in 'subtrees' with 'tips'
	if(class(tips)=="phylo") tips=tips$tip.label
	res=lapply(subtrees, function(cur) {
		   z=sapply(cur, function(y) {
					Ny=Ntip(y)
					Nt=length(tips)
					return((sum(y$tip.label%in%tips)*2)/(Nt+Ny))
		   })
		   return(z[z!=0])
	})
	names(res)=names(subtrees)
	res
}


nodelookup.phylo=function(phy, label){
	N=Ntip(phy)
	ww=which(phy$node.label==label)
	if(length(ww)) {
		return(N+ww)
	} else {
		warning(paste("Encountered no node.label for ",label,sep=""))
		return(NULL)
		
	}
}


# utility for sinking nodes by node.label or by ID
nodesink.phylo = function (phy, node=NULL, label=NULL) {
# label: a character string matching tip.label or node.label
# node: numeric id matching in phy$edge[,2]
	
    if ((is.null(node) & is.null(label)) | (is.numeric(node) & all(is.character(label)))) 
	stop("Please supply either 'node' or 'label'")
    if (class(phy) != "phylo") 
	stop("The tree does not appear to be a valid phylo object")
	
	# find 'nodes' of each 'label'
	if(is.null(node) & all(is.character(label))){
		labs=c(phy$tip.label, phy$node.label)
		if(!all(label%in%labs)){
			stop(paste(paste(label[which(!label%in%labs)], collapse=", "), "could not be matched to any branches in 'phy'"))				
		} 
		node=unique(unlist(lapply(label, function(x){
								  ww=which(labs==x)
								  if(length(ww)>1){
								  warning(paste(x, "is encountered for multiple branches in 'phy'"))
								  }
								  ww
								  })))			
	}
	
	ind = node
    if (any(ind == Ntip(phy) + 1)) {
        ind = ind[-which(ind == Ntip(phy) + 1)]
        warning("Root node will be left unsunk")
    }
    
	if (any(ind <= Ntip(phy))) {
     	tips.drop=ind[ww<-which(ind <= Ntip(phy))]
        ind = ind[-ww]
        warning(paste(paste(tips.drop, collapse=", "), "encountered in tips of 'phy'"))
    } else {
    	tips.drop=NULL
    }
	
	if(!all(ind%in%phy$edge[,2])) {
		stop(paste(paste(ind[which(!ind%in%phy$edge[,2])], collapse=", "), "could not be matched to any branches in 'phy'"))		
	}
	
    n <- length(ind)
    
    if (n) {
        ind.tmp = match(ind, phy$edge[, 2])
        ind = ind.tmp[!is.na(ind.tmp)]
		orig.edge = phy$edge
		orig.phy = phy
		ntips = Ntip(phy)
		
		## primary function for modifying 'edge' matrix
		reedge <- function(ancestor, des.to.drop) {
			wh <- which(phy$edge[, 1] == des.to.drop)
			dd <- which(orig.edge[, 2] == des.to.drop)
			dropped.branch <- phy$edge.length[dd]
			d.d <- c(get.desc.of.node(des.to.drop, orig.phy))
			if (length(d.d)) 
			phy$edge.length[match(d.d, orig.edge[, 2])] <<- phy$edge.length[match(d.d, 
																				  orig.edge[, 2])] + dropped.branch
			for (k in wh) {
				if (phy$edge[k, 2] %in% node.to.drop) {
					reedge(ancestor, phy$edge[k, 2])
				}
				else {
					phy$edge[k, 1] <<- ancestor
				}
			}
		}
		node.to.drop <- phy$edge[ind, 2]
		anc <- phy$edge[ind, 1]
		for (i in 1:n) {
			if (anc[i] %in% node.to.drop) 
			next
			reedge(anc[i], node.to.drop[i])
		}
		phy$edge <- phy$edge[-ind, ]
		phy$edge.length <- phy$edge.length[-ind]
		phy$Nnode <- phy$Nnode - n
		sel <- phy$edge > min(node.to.drop)
		for (i in which(sel)) phy$edge[i] <- phy$edge[i] - sum(node.to.drop < 
															   phy$edge[i])
		if (!is.null(phy$node.label)) 
		phy$node.label <- phy$node.label[-(node.to.drop - length(phy$tip.label))]
		
    }
	if(!is.null(tips.drop)) phy=drop.tip(phy, tips.drop)
    phy
}


#general phylogenetic utility for returning the first (usually, unless a polytomy exists) two descendants of the supplied node 
#author: JM EASTMAN 2010

get.desc.of.node <-
function(node, phy) {
	return(phy$edge[which(phy$edge[,1]==node),2])
}


get.epochs<-function(phy) 
# ordered: tip to root
{
	ancestors=lapply(1:Ntip(phy), function(x) c(x,rev(sort(get.ancestors.of.node(x,phy)))))
	phylo=get.times(phy)
	tt=phylo$time[order(phylo$des)]
	e=lapply(ancestors, function(x) tt[x])
	names(e)=phy$tip.label
	return(e)
}

get.times<-function(phy) 
{
	df=data.frame(phy$edge)
	names(df)=c("anc","des")
	anc=lapply(df$des,get.ancestors.of.node,phy)
	order.edges=c(0,phy$edge.length)[order(c(Ntip(phy)+1,df$des))]
	df$time=0
	for(n in 1:length(anc)){
		df$time[n]=sum(order.edges[anc[[n]]])+order.edges[df$des[n]]
	}
	
	df=rbind(c(0,Ntip(phy)+1,0),df)
	
	return(df)
}

#general phylogenetic utility for determining whether a node is the root of the phylogeny
#author: JM EASTMAN 2011

is.root <-
function(node,phy) {
	if(node==Ntip(phy)+1) return(TRUE) else return(FALSE)
}



phy.timeintersect=function(phy, time=135, relative=FALSE){
	if (is.null(phy$edge.length)) 
	stop("The tree does not appear to have branch lengths")
    if (class(phy) != "phylo") 
	stop("The tree does not appear to be a valid phylo object")
    if (!is.ultrametric(phy))
	stop("'phy' must be ultrametric")
    orig=phy
    phy$node.label=NULL
    bb <- branching.times(phy)
    mm <- match(phy$edge[,1],names(bb))
    nt <- bb[mm]
    xx=as.data.frame(round(cbind(phy$edge, from=nt, to=nt-phy$edge.length),4))
    names(xx)=c("anc", "des", "from", "to")
    if(any(xx<0)) stop("FIXME: branching times are incorrect")
    tmp=apply(xx, 1, function(x) withinrange(time, x[["from"]], x[["to"]]))
    if(any(tmp)){
    	nd=xx$des[which(tmp)]
    	nd=nd[nd>Ntip(phy)]
    }
    N=Ntip(phy)
    res=lapply(nd, function(x) extract.clade(phy, x))
    nm=orig$node.label[nd-N]
	class(res)="multiPhylo"
    names(res)=ifelse(nm=="", paste("node_", nd, sep=""), nm)
    
## plotting ##
#   if(plot){
#  	plot(phy, show.tip=FALSE)
#	nodelabels(node=nd, frame="n", pch=21, bg=1)
#   }
	
    return(list(phy=res, nd=nd))
}

get.btimes <- function (phy) 
{
    df = data.frame(phy$edge)
    names(df) = c("anc", "des")
    anc = lapply(df$des, get.ancestors.of.node, phy)
    order.edges = c(0, phy$edge.length)[order(c(Ntip(phy) + 1, df$des))]
    df$time = 0
    for (n in 1:length(anc)) {
        df$time[n] = sum(order.edges[anc[[n]]]) + order.edges[df$des[n]]
    }
    df = rbind(c(0, Ntip(phy) + 1, 0), df)
	df$edge=c(0,phy$edge.length)
	df$time=max(df$time)-df$time
    return(df)
}


collect_times=function(phy){
	if(!is.ultrametric(phy)) stop("Must supply an ultrametric 'phy'.")
	dat=get.times(phy)
	dat$time=max(dat$time)-dat$time
	return(dat)
}

get.multifurcations=function(phy){
	nn=phy$edge[phy$edge[,2]>Ntip(phy),2]
	mm=sapply(nn, function(n) if(length(get.desc.of.node(n, phy))>2) return(TRUE) else return(FALSE))
	names(mm)=nn
	mm
}

relabel.phylo=function(taxonomy, phy, label="family"){
	tmp=lapply(phy$tip.label, function(o) unique(taxonomy[which(taxonomy==o, arr=TRUE)[,1],label]))
	tips=sapply(tmp, function(t) if(length(t)==1) return(t) else return("null"))
	phy$tip.label=tips
	phy=drop.tip(phy, phy$tip.label[tips=="null"])
	return(phy)
}

polytomy.phylo=function(tips, age=1){
	N=length(tips)+1
	edge=cbind(N,1:(N-1))
	length=rep(age, N-1)
	phy=list(edge=edge, edge.length=length, tip.label=tips, Nnode=1)
	class(phy)="phylo"
	return(phy)
}

subset.phylo=function(phy){
# phy: has tip.labels that are not unique
# returns phylogeny with unique tip labels
# if a taxon (indicated by multiple tip labels of same string) is non-monophyletic, all tips of that taxon are removed with warning
# likely used where an exemplar phylogeny is needed and phy$tip.label is 'tricked' into a non-unique vector of names
	tt=table(phy$tip.label)
	if(!any(tt>1)) {
		return(phy)
	} else {
		todo=names(tt[tt>1])
		cat("checking monophyly of groups...\n")
		mon=sapply(todo, function(t) {
				   cat(paste("\n\t",t, sep=""))
				   exemplar.monophyly(t, phy)
		})
		if(any(!mon)) phy=drop.tip(phy, names(!mon))
		return(unique.phylo(phy))
	}
}

unique.phylo=function(phy){
	# phy is assumed to have tips that are redundant (and whose exemplars are monophyletic)
	# prunes tree to leave one member of each unique label
	
	mon=table(phy$tip.label)
	if(all(mon==1)) return(phy)
	mon=mon[mon>1]
	for(i in 1:length(mon)){
		m=names(mon[i])
		drop=which(phy$tip.label==m)
		drop=drop[2:length(drop)]
		phy$tip.label[drop]="null"
	}
	phy=drop.tip(phy, phy$tip.label[phy$tip.label=="null"])
	return(phy)
}

exemplar.monophyly=function(tip, phy) {
# tip: occurs multiply in phy$tip.label
	nn=which(phy$tip.label==tip)
	if(length(nn)==1) return(TRUE) 
	anc=getMRCA(phy, nn)
	dd=unique(phy$tip.label[get.descendants.of.node(anc, phy, tips=TRUE)])
	if(length(dd)==1) return(TRUE) else return(FALSE)
}

nodelabel.phylo=function(phy, taxonomy, rank=NULL){
# all phy$tip.label must be in taxonomy
# rank: a column in 'taxonomy' showing where phy$tip.label are expected; if NULL, rownames are assumed to be rank
# taxonomy: exclusivity highest on left, lowest on right (species, genus, family, etc., as columns)
# columns in 'taxonomy' should ONLY be taxonomic ranks
	
## FIXME: add multicore support
## FIXME: warn on missing labels (due to non-monophyly)

	if(is.null(rank)){
		taxonomy=cbind(rownames=rownames(taxonomy),taxonomy)
		rank="rownames"
	}
	
	taxonomy=as.data.frame(as.matrix(taxonomy),stringsAsFactors=FALSE)
	
	oo=order(apply(taxonomy, 2, function(x) length(unique(x))),decreasing=TRUE)
	if(!all(oo==c(1:ncol(taxonomy)))){
		warning("Assuming 'taxonomy' is not from most to least exclusive")
		taxonomy=taxonomy[,oo]
	}
	
	if(!all(xx<-phy$tip.label%in%taxonomy[,rank])) {
		warning(paste("taxa not found in 'taxonomy':\n\t", paste(phy$tip.label[!xx], collapse="\n\t"), sep=""))
	}
	taxonomy[taxonomy==""]=NA
	options(expressions=50000)
	
	unmatched=phy$tip.label[!xx]
	idx=match(phy$tip.label, taxonomy[,rank])
	tax=taxonomy[idx,]
	
	labels=unique(unlist(tax[,-which(names(tax)==rank)]))
	labels=labels[!is.na(labels)]
	
	cat("resolving descendants of groups...\n\t")
	dat=tax[,-which(names(tax)==rank)]
	hashes_labels=character(length(labels))
	zz=tax[,rank]
	index=0
	tips=phy$tip.label[xx]
	for(i in 1:ncol(dat)){
		uu=unique(dat[,i])
		uu=uu[!is.na(uu)]
		for(j in uu){
			index=index+1
			cur=zz[which(dat[,i]==j)]
			if(length(cur)>1){
				hashes_labels[index]=hash.tip(cur, tips)
			} 		
		}
	}
	names(hashes_labels)=labels
	if(any(hashes_labels=="")) hashes_labels=hashes_labels[-which(hashes_labels=="")]
	tmp=table(hashes_labels)
	if(any(tmp[tmp>1]->redund)){
		for(r in names(redund)){
			rdx=which(hashes_labels==r)
			warning(paste("redundant labels encountered:\n\t", paste(names(hashes_labels[rdx]), collapse="\n\t"), sep=""))
			hashes_labels=hashes_labels[-rdx[!rdx==max(rdx)]]
		}
	}
	
	cat("\n\nresolving descendants for splits in tree...\n\t")
	tmp=hashes.phylo(phy, phy$tip.label[xx])
	hashes_tree=tmp$hash
	phy$node.label=rep("",max(phy$edge))
	mm=match(hashes_tree, hashes_labels)
	nodelabels=ifelse(is.na(mm), "", names(hashes_labels[mm]))
	nodelabels[is.na(nodelabels)]=""
	nodelabels=nodelabels[(Ntip(phy)+1):max(phy$edge)]
	tmp=table(nodelabels[nodelabels!=""])
	if(any(tmp[tmp>1]->redund)){
		for(r in names(redund)){
			rdx=which(nodelabels==r)
			nodelabels[rdx[rdx!=min(rdx)]]=""
		}
	}
	
	phy$node.label=nodelabels
	phy
}


root.phylo=function(phy, outgroup, taxonomy=NULL){
## GENERAL FUNCTION FOR ROOTING (based on outgroup)
# taxonomy: classification data.frame with 'species' minimally as a column

	sys=taxonomy
	if(!is.null(sys)) {
		rows=unique(unlist(sapply(outgroup, function(o) which(sys==o, arr=TRUE)[,1])))
		outgroup=sys$species[rows]
		outgroup=outgroup[outgroup%in%phy$tip.label]
	} else {
		if(!all(outgroup%in%phy$tip.label)) stop("Some 'outgroup' appear missing from 'phy'.")
	}
	
	tips=match(outgroup, phy$tip.label)
	node=getMRCA(phy,tips)
	if(node==Ntip(phy)+1 & !is.null(sys)){
		node=getMRCA(phy, (1:Ntip(phy))[-tips])
	}
	rooted=root(phy, node=node, resolve.root=TRUE)
	rooted
}



polytomy_sampler=function(phy){
	
# stores node labels from original phylogeny
	phy=reorder(phy)
	if(is.null(phy$node.label)) phy$node.label=sort(unique(phy$edge[,1]))
	
	sampler=function(node, phy, N){
		dd=get.desc.of.node(node, phy)
		if(length(dd)>2) {
			tips=dd[dd<=N]
			if(length(tips)>1){
				drop=tips[sample(1:length(tips), length(tips)-1)]
				return(list(drop=drop, keep=tips[!tips%in%drop]))
			}
		}
	}
	
	fill_mapping=function(mapping, todo, phy){
		keep=sapply(todo, "[[", "keep")
		for(i in 1:length(keep)){
			if(!is.null(keep[[i]])){
				taxon=phy$tip.label[keep[[i]]]
				tips.idx=todo[[i]]$drop
				mapping[[taxon]]=unique(c(mapping[[taxon]], phy$tip.label[tips.idx], unlist(lapply(phy$tip.label[tips.idx], function(x) mapping[[x]]))))
			}
		}
		return(mapping)
	}
	
	mapping=vector("list", Ntip(phy))
	names(mapping)=phy$tip.label
	
	while(1){
		N=Ntip(phy)
		todo=lapply((N+1):max(phy$edge), function(n) sampler(n , phy, N))
		tips=unlist(sapply(todo, "[[", "drop"))
		if(length(tips)){
			mapping=fill_mapping(mapping, todo, phy)
			phy=drop.tip(phy, tips)
		} else {
			break()
		}
	}
	mapping=mapping[names(mapping)%in%phy$tip.label]
	return(list(phy=phy, mapping=mapping))
}


is.monophyletic.for.group=function(node, tips, phy) {
	if(node>Ntip(phy)) {
		if(all(phy$tip.label[get.descendants.of.node(node, phy, tips=TRUE)]%in%tips))return(TRUE) else return(FALSE)
	} else {
		if(phy$tip.label[node]%in%tips) return(TRUE) else return(FALSE)
	}
}


#general phylogenetic utility for tabulating non-overlapping lineages (by node reference), descended from a node and given a taxon set (tips)
#counts clusters of monophyly (singletons count for 1 as well as nested monophyletic groupings; a paraphyletic group (A(B(C,1))) will count for 3; a paraphyletic group (A(1(C,B))) will count for 2)
#author: JM EASTMAN 2011

proc.polyphyly <-
function (node, phy, tips) 
{
    node <- as.numeric(node)
    n <- Ntip(phy)
	
	if(is.monophyletic.for.group(node, tips, phy)) return(1)
	
	l <- c()
    d <- get.desc.of.node(node, phy)
    for (j in d) {
		l <- c(l, proc.polyphyly(j, phy, tips))
    }
    return(sum(l))
}


#general phylogenetic utility for tabulating non-overlapping lineages (by node reference), descended from a node and given a taxon set (tips)
#computes proportion of edges that belong to a group [... strict=TRUE: minimal path lengths for the focal group; strict=FALSE: exclude tips only
#author: JM EASTMAN 2011

proc.monophyly=function(tips, phy) {
	yy=phy$edge[,2]; yy[]=0
	
	if(is.monophyletic(phy, tips, reroot=FALSE)) {
		val=1
		nn=collect.nodes(tips, phy, strict=TRUE)
		yy[match(nn,phy$edge[,2])]=1
	} else {
		N=Ntip(phy)
		nn=lapply(c(TRUE, FALSE), function(x) {y=collect.nodes(tips, phy, strict=x)})
		ee=sapply(nn, function(x) {return(sum(phy$edge.length[match(x[x!=N+1], phy$edge[,2])]))})
		yy[match(nn[[1]],phy$edge[,2])]=1
		val=ee[1]/ee[2]
		nn=nn[[1]]
	}
	return(list(val=val, bin=yy, nodes=nn))
}

drop.zerotips=function(phy, threshold=0, all=FALSE){
# all=FALSE: retains one tip of a pair of sister species that each have 0-length tips
	
	t=phy
	N=Ntip(t)
	tipedges=t$edge.length[t$edge[,2]<=N]
	names(tipedges)=t$tip.label[t$edge[t$edge[,2]<=N,2]]
	w=which(tipedges<=threshold)
	if(length(w)){
		nn=match(names(w), t$tip.label)
		if(!all){
			anc=unique(sapply(nn, get.ancestor.of.node, t))
			drop=unlist(sapply(anc, function(x) {dd=get.desc.of.node(x, t); dd=dd[dd<=N]; dd[sample(1:length(dd), length(dd)-1)]}))
		} else {
			drop=nn
		}
		tre=drop.tip(t, drop)
	} else {
		tre=t
	}
	
	return(tre)
}

ultrametricize=function(phy, tol=1e-8, trim=c("min","max","mean"), depth=NULL){
	
	phy <- reorder(phy)
    n <- length(phy$tip.label)
    n.node <- phy$Nnode
    xx <- numeric(n + n.node)
    for (i in 1:dim(phy$edge)[1]) xx[phy$edge[i, 2]] <- xx[phy$edge[i, 1]] + phy$edge.length[i]
	
	paths=xx[1:n]
	trim=switch(match.arg(trim),
				min = min(paths),
				max = max(paths),
				mean = mean(paths))
	if(!is.null(depth)) trim=depth
	if(diff(range(paths))<=tol) {
		ww=which(phy$edge[,2]<=n)
		phy$edge.length[ww]=phy$edge.length[ww]+(trim-paths[phy$edge[ww,2]])
		return(phy)
	} else {
		stop("Difference in path lengths is larger than supplied tolerance")
	}
}

#general phylogenetic utility for finding the edge labels (from phy$edge[,2]) associated with the most exclusive monophyletic grouping containing at least supplied 'tips' where 'tips' are elements of the phy$tip.label vector
#author: JM EASTMAN 2010

collect.nodes <-
function(tips,phy,strict=FALSE){
	if(length(tips)==1) stop("Must supply at least two tips.")
	tt=phy$tip.label
	tt=lapply(tips, function(x) which(phy$tip.label==x))
	nodes=lapply(tt, function(x) get.ancestors.of.node(x, phy))
	all.nodes=nodes[[1]]
	for(i in 2:length(nodes)) {
		all.nodes=intersect(all.nodes,nodes[[i]])
	}
	base.node=max(all.nodes)
	if(strict) {
		u=unlist(nodes)
		return(c(sort(unique(u[u>=base.node])), unlist(tt)))
	} else {
		return(c(base.node,get.descendants.of.node(base.node,phy)))
	}
}

#general phylogenetic utility for finding the phy$tip.labels associated with the most exclusive monophyletic grouping containing at least supplied 'tips' where 'tips' are elements of the phy$tip.label vector
#author: JM EASTMAN 2010

find.relatives <-
function(tips,phy){
	nn=collect.nodes(tips, phy, strict=FALSE)
	return(phy$tip.label[nn[nn<=Ntip(phy)]])
}

#general phylogenetic utility for returning first ancestor (as a numeric referent) of the supplied node 
#author: JM EASTMAN 2010

get.ancestor.of.node <-
function(node, phy) {
	return(phy$edge[which(phy$edge[,2]==node),1])
}

#general phylogenetic utility for returning all ancestors of nodes in a phylogeny
#author: JM EASTMAN 2011

get.ancestors <-
function(phy) {
	r=lapply(1:Ntip(phy), function(x) c(x,rev(sort(get.ancestors.of.node(x,phy)))))
	r
}



#general phylogenetic utility for returning all ancestors (listed as given in phy$edge[,2]) of a node 
#author: JM EASTMAN 2010

get.ancestors.of.node <-
function(node, phy) {
	a=c()
	if(node==(Ntip(phy)+1->root)) return(NULL)
	f=get.ancestor.of.node(node, phy)
	a=c(a,f)
	if(f>root) a=c(a, get.ancestors.of.node(f, phy))
	return(a)
}

#general phylogenetic utility for returning all descendants (listed as given in phy$edge[,2]) of a node (excluding the supplied node) 
#author: JM EASTMAN 2011

get.descendants.of.node <- 
function(node, phy, tips=FALSE){
	n=Ntip(phy)
	all=ifelse(tips, FALSE, TRUE)
	out <- .Call("get_descendants", tree=list(
											  NODE = as.integer(node),
											  ROOT = as.integer(n+1),
											  ALL = as.integer(all),
											  ENDOFCLADE = as.integer(dim(phy$edge)[1]),
											  ANC = as.integer(phy$edge[,1]),
											  DES = as.integer(phy$edge[,2])),
				 PACKAGE = "phylo")
	res=out$TIPS
	if(!length(res)) res=NULL
	return(res)
} 


# alternate version of "prop_part" if only tips are needed [FIXME: put into C++]
.compile_descendants=function(phy){
	pphy=reorder(phy,"cladewise")
	ee=pphy$edge
	ee=ee[nrow(ee):1,]
	pphy$edge=ee
	desc=as.list(1:(nrow(ee)+1))
	N=Ntip(phy)
	dat=ee[ee[,2]>N,]
	if(check.multicore()) f=mclapply else f=lapply
	descendants=f(dat[,2], function(x) get.desc.of.node(x, pphy))
	for(i in 1:nrow(dat)){
		dd=dat[i,2]
		desc[[dd]]=unlist(desc[descendants[[i]]])
	}
	desc[[N+1]]=1:N
	return(desc)
}

## FUNCTIONS
## grabs most exclusive tips from clade definitions (and whose tips can then be reconciled with a taxonomic database -- a lookup table)
# allows for recursion in clade definitions (clade defined in part using another clade definition)
# returns trees representing each clade definition
# 'nested': an important variable -- this is the subset of clades that are defined (recursively) within 'clades'
phy_cladedefinition=function(clades, unplaced=TRUE){
	
	cc=unique(unlist(c(clades)))
	tt=cc%in%names(clades)
	nested=cc[tt]
	
	## FUNCTION for testing whether a clade is 'nested' within 'clades'	
	is.nestedclade=function(clade, nested){
		if(clade%in%nested) return(TRUE) else return(FALSE)
	}
	
	
	## FUNCTION for grabbing the set of nested groups within 'clades'
	fetch_nestedclades=function(clade, clades, nested){
		
		if(is.nestedclade(clade, nested)){
			desc=clades[[clade]]
			new=as.list(desc)
			names(new)=desc
			tt=sapply(new, is.nestedclade, nested)
			for(i in 1:length(new)){
				if(tt[i]) new[[i]]=fetch_nestedclades(new[[i]], clades, nested) else new[[i]]=NA
			}
		} else {
			new=NA
		}
		return(new)
	}
	
	## FUNCTION for finding nestedness of clades within clades (and their path)
	paths_through_clades=function(clades, nested){
		nn=names(clades)
		res=lapply(nn, function(x) {
				   dd=clades[[x]]
				   y=lapply(dd, fetch_nestedclades, clades, nested)
				   names(y)[1:length(dd)]=dd
				   y
				   })
		names(res)=nn
		res
	}
	
	cladepaths=paths_through_clades(clades, nested)
	if(unplaced){
		for(i in 1:length(cladepaths)){
			cur=cladepaths[[i]]
			if(!is.null(unplc<-attributes(clades[[i]])$unplaced)){
				unplaced_taxa=lapply(unplc, function(x) return(NA))
				names(unplaced_taxa)=unplc
				cladepaths[[i]]$unplaced=unplaced_taxa
			}
		}
	}
	
	unplaced_phy=function(phy, cladepath){
		if(any(names(cladepath)=="unplaced")){
			tips=names(cladepath$unplaced)
			y=polytomy.phylo(tips)
			y$edge.length=NULL
			new=bind.tree(phy, y)
			new=drop.tip(new, "unplaced")
			return(new)
		} else {
			return(phy)
		}
	}
	
	## FUNCTION for writing a Newick tree from cladepath
	tree_cladepath=function(cladepath, nested){
		require(ape)
		
		# write newick string
		print_group=function(cladepath) {
			xx=sapply(names(cladepath), is.nestedclade, nested)
			middle=character(length(xx))
			for(i in 1:length(xx)){
				if(xx[i]) {
					new=cladepath[[names(xx)[i]]]
					middle[i]=print_group(new)
				} else {
					middle[i]=names(cladepath)[i]
				}
			}
			paste("(", paste(middle, collapse=", "), ")", sep="")
		}
		
		tmp=paste(print_group(cladepath),";",sep="")
		phy=read.tree(text=tmp)
		return(unplaced_phy(phy, cladepath))
	}
	
	phy=lapply(cladepaths, tree_cladepath, nested)
	
	## USE node.depth() to find maximum node depth in tree
	
	return(phy)
}


build_lookup.phylo.DEFUNCT=function(phy, taxonomy, clades=NULL){
## taxonomy expected to have first column at same level as tip labels in phy
## first row in taxonomy is most exclusive
## clade_defs are phylogenetic trees of a clade representation
	
	if(!any(taxonomy[,1]%in%phy$tip.label) & !is.null(rownames(taxonomy))) {
		taxonomy=as.data.frame(as.matrix(cbind(species=rownames(taxonomy), taxonomy)),stringsAsFactors=FALSE)
	}
	
	related_tips=function(tips, phy){
		tips=tips[tips%in%phy$tip.label]
		if(length(tips)<2) stop("related_tips(): 'tips' to be found in 'phy' are too few")
		nd=unique(sapply(tips, function(x) match(x, phy$tip.label)))
		if(length(nd)==1) return(nd)
		anc=getMRCA(phy, nd)
		dd=get.descendants.of.node(anc, phy, tips=TRUE)
		return(phy$tip.label[dd])
	}
	
	tips_in_group=function(phy, taxonomy, clade_def){
		tmp=sapply(clade_def$tip.label, function(grp){
				ww=which(taxonomy==grp, arr.ind=TRUE)[,"row"]
				if(length(ww)>0){
				   return(unique(taxonomy[ww,1]))
				} else {
				   return(c())
				}  
		})
		
	## most exclusive 'spp' (recognized from tree)
		spanning_taxa=unique(unlist(tmp))
		recovered_taxa=unique(related_tips(spanning_taxa, phy))
		
		return(recovered_taxa)
	}
	
	require(phylo)
	orig_phy=phy
	
	## check ordering of taxonomy
	oo=order(apply(taxonomy, 2, function(x) length(unique(x))),decreasing=TRUE)
	if(!all(oo==c(1:ncol(taxonomy)))){
		warning("Assuming 'taxonomy' is not from most to least exclusive")
		taxonomy=taxonomy[,oo]
	}
	original_taxonomy=taxonomy
	
	## check rank of tree
	phy$node.label=NULL
	gtips=sapply(phy$tip.label, function(x) unlist(strsplit(gsub(" ", "_", x), "_"))[1])
	check=sapply(list(phy$tip.label, gtips), function(x) sum(x%in%taxonomy[,1]))
	if(max(check)==check[2]) {
		genuslevel_tree=TRUE
		phy$tip.label=unname(gtips)
	} else {
		genuslevel_tree=FALSE
	}
	tips=phy$tip.label
	
	## prune taxonomy down to members in the phylogeny
	taxonomy=taxonomy[taxonomy[,1]%in%tips,]
	matching_tips=taxonomy[,1]
	
	## BUILD taxonomic mapping to each row in 'phy'
	mm=match(tips, matching_tips)
	tt=taxonomy[mm,]
	colnames(tt)=names(taxonomy)
	
	if(!is.null(clades)){
		clade_defs=phy_cladedefinition(clades)
		cat("resolving clades...\n\t")
		res=lapply(1:length(clade_defs), function(idx) {
				   def=clade_defs[[idx]]
				   cat(paste(names(clade_defs)[idx], "\n\t", sep=""))
				   tt=try(tips_in_group(phy, taxonomy, def),silent=TRUE)
				   if(inherits(tt, "try-error")) return(NA) else return(tt)
				   }) 
		names(res)=names(clade_defs)
		
		zz=sapply(res, function(x) all(is.na(x)))
		if(any(zz)){
			warning(paste("taxa not represented in 'tips':\n\t", paste(names(res)[zz], collapse="\n\t"), sep=""))
		}
		res=res[!zz]
		clade_defs=clade_defs[!zz]
		
	## BUILD clade-level mapping to each row in 'phy'
		gg=matrix("", nrow=Ntip(phy), ncol=length(res))
		for(i in 1:length(res)){
			nm=names(clade_defs)[i]
			cur_tips=res[[i]]
			zz=tips%in%cur_tips
			gg[zz,i]=nm
		}
		colnames(gg)=names(clade_defs)
		phy_mapping=cbind(tt, gg, stringsAsFactors=FALSE)
	} else {
		phy_mapping=tt
	}
	
	rownames(phy_mapping)=orig_phy$tip.label
	return(as.data.frame(as.matrix(phy_mapping), stringsAsFactors=FALSE)) 	
}

build_lookup.phylo=function(phy, taxonomy=NULL, clades=NULL){
## taxonomy expected to have first column at same level as tip labels in phy
## first row in taxonomy is most exclusive
## clade_defs are phylogenetic trees of a clade representation
	
	if(!is.null(taxonomy)){
		if(!any(taxonomy[,1]%in%phy$tip.label) & !is.null(rownames(taxonomy))) {
			taxonomy=as.data.frame(as.matrix(cbind(species=rownames(taxonomy), taxonomy)),stringsAsFactors=FALSE)
		}
	}
	
	
	related_tips=function(tips, phy){
		tips=tips[tips%in%phy$tip.label]
		if(length(tips)<2) stop("related_tips(): 'tips' to be found in 'phy' are too few")
		nd=unique(sapply(tips, function(x) match(x, phy$tip.label)))
		if(length(nd)==1) return(nd)
		anc=getMRCA(phy, nd)
		dd=get.descendants.of.node(anc, phy, tips=TRUE)
		return(phy$tip.label[dd])
	}
	
	tips_in_group=function(phy, taxonomy=NULL, clade_def){
		
		lengths=sapply(clade_def$tip.label, function(grp){
						if(!is.null(taxonomy)){
						   ww=which(taxonomy==grp, arr.ind=TRUE)[,"row"]
						   if(length(ww)>0){
						   return(TRUE)
						   } else {
						   return(FALSE)
						   } 
						} else {
							return(length(phy$tip.label[which(phy$tip.label%in%grp)])>0)
						}
		})
		
		if(sum(lengths)<2) return(NA)
		
		## FIXME: use 'lengths' to determine if 'clade_def' is satisfied (at least two members needed)
		
		
		tmp=sapply(clade_def$tip.label, function(grp){
				   if(!is.null(taxonomy)){
					ww=which(taxonomy==grp, arr.ind=TRUE)[,"row"]
					if(length(ww)>0){
						return(unique(taxonomy[ww,1]))
					} else {
						return(c())
					} 
				   } else {
					return(phy$tip.label[which(phy$tip.label%in%grp)])
				   }
		})
		
		## most exclusive 'spp' (recognized from tree)
		spanning_taxa=unique(unlist(tmp))
		recovered_taxa=unique(related_tips(spanning_taxa, phy))
		
		return(recovered_taxa)
	}
	
	require(phylo)
	orig_phy=phy
	
	if(!is.null(taxonomy)){
		## check ordering of taxonomy
		oo=order(apply(taxonomy, 2, function(x) length(unique(x))),decreasing=TRUE)
		if(!all(oo==c(1:ncol(taxonomy)))){
			warning("Assuming 'taxonomy' is not from most to least exclusive")
			taxonomy=taxonomy[,oo]
		}
		original_taxonomy=taxonomy
		
		## check rank of tree
		phy$node.label=NULL
		gtips=sapply(phy$tip.label, function(x) unlist(strsplit(gsub(" ", "_", x), "_"))[1])
		check=sapply(list(phy$tip.label, gtips), function(x) sum(x%in%taxonomy[,1]))
		if(max(check)==check[2]) {
			genuslevel_tree=TRUE
			phy$tip.label=unname(gtips)
		} else {
			genuslevel_tree=FALSE
		}
		tips=phy$tip.label
		
		## prune taxonomy down to members in the phylogeny
		taxonomy=taxonomy[taxonomy[,1]%in%tips,]
		matching_tips=taxonomy[,1]
		
		## BUILD taxonomic mapping to each row in 'phy'
		mm=match(tips, matching_tips)
		tt=taxonomy[mm,]
		colnames(tt)=names(taxonomy)
	} else {
		tt=c()
	}
	
	if(!is.null(clades)){
		tips=phy$tip.label
		clade_defs=phy_cladedefinition(clades)
		cat("resolving clades...\n\t")
		res=lapply(1:length(clade_defs), function(idx) {
				   def=clade_defs[[idx]]
				   cat(paste(names(clade_defs)[idx], "\n\t", sep=""))
				   tt=try(tips_in_group(phy, taxonomy, def),silent=TRUE)
				   if(inherits(tt, "try-error")) return(NA) else return(tt)
				   }) 
		names(res)=names(clade_defs)
		
		zz=sapply(res, function(x) all(is.na(x)))
		if(any(zz)){
			warning(paste("taxa not represented in 'tips':\n\t", paste(names(res)[zz], collapse="\n\t"), sep=""))
		}
		res=res[!zz]
		clade_defs=clade_defs[!zz]
		
## BUILD clade-level mapping to each row in 'phy'
		gg=matrix("", nrow=Ntip(phy), ncol=length(res))
		for(i in 1:length(res)){
			nm=names(clade_defs)[i]
			cur_tips=res[[i]]
			zz=tips%in%cur_tips
			gg[zz,i]=nm
		}
		colnames(gg)=names(clade_defs)
		tt=as.matrix(cbind(tt, gg))
	} 
	
	phy_mapping=tt
	rownames(phy_mapping)=orig_phy$tip.label
	return(as.data.frame(as.matrix(phy_mapping), stringsAsFactors=FALSE)) 	
}



## FUNCTIONS ##

rank_prune.phylo=function(phy, taxonomy, rank="family"){
## rank (e.g., 'family') and 'genus' must be in columns of 'taxonomy'
	
	if(!rank%in%colnames(taxonomy)){
		stop(paste(sQuote(rank), " does not appear as a column name in 'taxonomy'", sep=""))
	}
	
	xx=match(phy$tip.label, rownames(taxonomy))
	
	new=as.matrix(cbind(tip=phy$tip.label, rank=taxonomy[xx,rank]))
	drop=apply(new, 1, function(x) if( any(is.na(x)) | any(x=="")) return(TRUE) else return(FALSE))
	if(any(drop)){
		warning(paste("Information for some tips is missing from 'taxonomy'; offending tips will be pruned:\n\t", paste(phy$tip.label[drop], collapse="\n\t"), sep=""))
		phy=drop.tip(phy, phy$tip.label[drop])
		new=new[!drop,]
	}
	
	tips=phy$tip.label
	hphy=hashes.phylo(phy, tips=tips)
	tax=as.data.frame(new, stringsAsFactors=FALSE)
	stax=split(tax$tip,tax$rank)
	rank_hashes=sapply(stax, function(ss) hash.tip(ss, tips=tips))
	
	pruned=hphy
	pruned$tip.label=tax$rank

	if(!all(zz<-rank_hashes%in%hphy$hash)){
		warning(paste(paste("non-monophyletic at level of ",rank,sep=""),":\n\t", paste(sort(nonmon<-names(rank_hashes)[!zz]), collapse="\n\t"), sep=""))
		pruned=drop.tip(pruned, nonmon)
	}
		
	rank_phy=unique.phylo(pruned)
	
	return(rank_phy)	
}

rank_prune.phylo.DEFUNCT=function(phy, taxonomy, rank="family"){
## rank (e.g., 'family') and 'genus' must be in columns of 'taxonomy'
	
	if(!rank%in%colnames(taxonomy)){
		stop(paste(sQuote(rank), " does not appear as a column name in 'taxonomy'", sep=""))
	}
	
	taxonomy$rownames=rownames(taxonomy)

	xx=sapply(phy$tip.label, function(x) {
			  y=which(taxonomy==x, arr.ind=TRUE)[,"row"]
			  if(!length(y)==1) {
				tmp=split_lab(x)[1]
				yy=which(taxonomy==tmp, arr.ind=TRUE)[,"row"]
				if(!length(yy)==1){
					return(NA)
				} else {
					return(yy)
				}
			  } else {
				return(y)
			  }
	})
	
	new=as.matrix(cbind(tip=phy$tip.label, rank=taxonomy[xx,rank]))
	drop=apply(new, 1, function(x) if( any(is.na(x)) | any(x=="")) return(TRUE) else return(FALSE))
	if(any(drop)){
		warning(paste("Information for some tips is missing from 'taxonomy'; offending tips will be pruned:\n\t", paste(phy$tip.label[drop], collapse="\n\t"), sep=""))
		phy=drop.tip(phy, phy$tip.label[drop])
		new=new[!drop,]
	}
	
	mm=match(phy$tip.label, new[,"tip"])
	pruned=phy
	pruned$tip.label=new[mm,"rank"]
	
	rank_phy=subset.phylo(pruned)
	
	retained_tips=rank_phy$tip.label
	if(!all(zz<-pruned$tip.label%in%retained_tips)){
		warning(paste(paste("non-monophyletic at level of ",rank,sep=""),":\n\t", paste(unique(pruned$tip.label[!zz]), collapse="\n\t"), sep=""))
	}
	
	return(rank_phy)	
}


build_constraint.phylo=function(phy, taxonomy, rank="family"){
## phy: a 'rank' level phylogeny (tips of 'phy' should be matchable to taxonomy[,rank])
## taxonomy: a mapping from genus, family, order (columns in that order); rownames are tips to be added to constraint tree
##		-- 'taxonomy' MUST absolutely be ordered from higher exclusivity to lower (e.g., genus to order)
## rank: rank at which groups are assumed to be monophyletic (currently for 'family' only)
## returns a nodelabeled constraint tree based on 'phy' and 'rank'-level constraints
	
	require(phylo)
	
	tax=taxonomy
	nn=which(names(tax)==rank)
	tax=tax[,nn:ncol(tax)]
	tips=rownames(tax)
	
	# PRUNE 'rank'-level tree if some families unmatched in 'tips'
	if(any(zz<-!phy$tip.label%in%tax[,rank])){
		warning(paste("taxa not represented in 'tips':\n\t", paste(phy$tip.label[zz], collapse="\n\t"), sep=""))
		fam=drop.tip(phy, phy$tip.label[zz])
	}
	
	## FIX NODE LABELS (drop.tip() not trustworthy for keeping node.labels straight)
	apgtax=get(data(spermatophyta_APGIII_taxonomy))
	apgcld=get(data(spermatophyta_APGIII_clades))
	apglkp=build_lookup.phylo(fam, apgtax[,(match(rank, names(apgtax))):ncol(apgtax)], apgcld)
	fam=nodelabel.phylo(fam, apglkp)

	tmp=unique(c(fam$tip.label, fam$node.label))
	all_labels=tmp[!is.na(tmp)&tmp!=""]
	exclude=apply(tax, 1, function(x) !any(x%in%all_labels))
	missing_tips=rownames(tax)[exclude]
	warning(paste("tips missing data in 'taxonomy' and excluded from constraint tree:\n\t", paste(missing_tips, collapse="\n\t"), sep=""))
	tax=tax[!exclude,]
	
	# find tips that have data but whose data not at 'rank' level (plug in deeper in tree)
	at_rank=tax[,rank]%in%fam$tip.label
	deeper_tips=tax[!at_rank,]
	if(nrow(deeper_tips)){
		ww=apply(deeper_tips, 1, function(x) x[[min(which(x%in%all_labels))]])
		for(i in 1:length(ww)){
			nn=nodelookup.phylo(fam, ww[i])
			tmp=polytomy.phylo(names(ww[i]))
			fam=bind.tree(fam, tmp, where=nn)			
		}
		fam=compute.brlen(fam, method="Grafen")
		warning(paste("tips missing data at 'rank' level in 'taxonomy' but included in constraint tree:\n\t", paste(rownames(deeper_tips), collapse="\n\t"), sep=""))
	}
	fam$node.label=NULL
	tax=tax[at_rank,]
	if(any(is.na(tax[,rank]))) stop("Corrupted data encountered when checking taxonomy[,rank]")
	
	## CREATE 'rank'-level subtrees
	ss=split(tax, tax[,rank])
	mm=min(fam$edge.length)
	subtrees=lapply(1:length(ss), function(idx) {phy=polytomy.phylo(rownames(ss[[idx]]), mm/2); phy$root.edge=0; return(list(subtree=phy, age=mm/2, tip=names(ss)[idx]))})
	names(subtrees)=names(ss)
	
	## PASTE in 'rank'-level subtrees
	contree=glomogram.phylo(fam, subtrees)
	contree
}

.build_constraint.phylo.DEFUNCT=function(phy, tips, taxonomy, rank="family"){
## phy: a genus-level phylogeny
## tips: species to be added to constraint tree
##		-- labels in tips are expected to be subsettable by 'Genus_species'
## taxonomy: a mapping from genus, family, order (columns in that order)
## returns a nodelabeled constraint tree based on 'phy' and 'rank'-level constraints
## rank: rank at which groups are assumed to be monophyletic
	
	require(phylo)
	
	if(all(names(taxonomy)!=c("genus","family","order"))) stop("Supply 'taxonomy' with 'genus','family','order' as columns")
	
	tax=taxonomy
	seq=tips
	
	# PRUNE tree to 'rank' level
	fam=rank_prune.phylo(phy, tax, rank=rank)
	
	## create taxonomy for tips in seq
	seq_tax=matrix(NA, nrow=length(seq), ncol=ncol(tax))
	for(i in 1:length(seq)){
		nm=seq[i]
		gen=unlist(strsplit(nm, "_"))[1]
		if(any(tt<-tax$genus%in%gen)){
			ww=which(tt)
			if(length(ww)==1){
				seq_tax[i,]=unlist(tax[ww,])
			}
		}
	}
	seq_tax[is.na(seq_tax)]=""
	rownames(seq_tax)=seq
	colnames(seq_tax)=names(tax)
	seq_tax=data.frame(seq_tax, stringsAsFactors=FALSE)
	
	# PRUNE 'rank'-level tree if some families unmatched in 'tips'
	if(any(zz<-!fam$tip.label%in%seq_tax[,rank])){
		warning(paste("taxa not represented in 'tips':\n\t", paste(fam$tip.label[zz], collapse="\n\t"), sep=""))
		fam=drop.tip(fam, fam$tip.label[zz])
	}
	
	## CREATE 'rank'-level subtrees
	tmp=unique(seq_tax[,rank])
	unique_families=tmp[!tmp==""]
	tree_families=fam$tip.label[fam$tip.label%in%unique_families]
	
	mm=min(fam$edge.length)
	
	subtrees=list()
	subtreenames=c()
	t=0
	
	## CREATE SUBTREES for each 'rank' exemplar
	for(f in tree_families){
		sub=which(seq_tax[,rank]==f)
		tips=rownames(seq_tax)[sub]
		if(length(tips)==1){
			fam$tip.label[fam$tip.label==f]=tips
		} else {
			t=t+1
			stree=polytomy.phylo(tips, age=mm/2)
			stree$root.edge=0
			subtrees[[t]]=list(subtree=stree, age=mm/2, tip=f)
			subtreenames[t]=f
		}
	}
	
	names(subtrees)=subtreenames
	contree=glomogram.phylo(fam, subtrees)
	
	## NODE-LABEL constraint tree
	seq_tax=cbind(species=rownames(seq_tax), seq_tax)
	seq_tax[seq_tax==""]=NA
	return(list(contree=contree, taxonomy=seq_tax))
}

