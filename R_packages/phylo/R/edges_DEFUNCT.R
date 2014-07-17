## DEFUNCT: get.hashes() see hashes.phylo()
get.hashes <- function (phy, tips=NULL){
## GENERAL FUNCTION: find hash tag for every edge in tree (using phy$tip.label or predefined set of 'tips')
#	returns list of hash tags from node 1:(Nnode(phy)+Ntip(phy))
# relatively fast (faster than hash_proppart and hash_edges)
	
	if(is.null(tips)) {
		tips=phy$tip.label
		tips.idx=1:length(tips)
	} else {
		tips.idx=match(tips, phy$tip.label)
	}
	storage.mode(phy$Nnode) <- "integer"
	descendants <- .Call("prop_part", list(phy), 1, TRUE, PACKAGE = "ape")
	hash_internal <- sapply(descendants, function(desc) get.hash(desc, tips.idx))
	hash_tips <- sapply(1:Ntip(phy), function(tip) get.hash(tip, tips.idx))
	hashes=c(hash_tips, hash_internal)
	attr(hashes, "tips")=tips
	return(hashes)
}

store.unique_edges=function(phy, clades, classification){
## FIXME
	
# GENERIC FUNCTION for compiling all unique edges and storing name attributes given 'clades' and 'classification' across set of trees
	
# collect phylogenies from auteur output
	class(phy)="multiPhylo"
	
# collect all unique edges
	u_edges=unique_edges(phy)
	u_edges=name.attribute.for.edges(u_edges, classification, clades)
	
# store unique edge identifiers for every node in every tree
	for(i in 1:length(phy)){
		phy[[i]]=store_edges(phy[[i]], u_edges)
	}
	return(list(phy=phy, edges=u_edges))
	
}


## DEFUNCT: replace with which()
species.in.taxon=function(classification.df, contains=c("","")){
## FIXME
# contains: name is a taxon, c() is the list of contained taxa (must be in classification.df)
# classification.df: must have 'species' as a column
	set=contains
	dat=classification.df
	spp=lapply(set, function(x) {
			   cols=apply(dat, 2, function(t) any(suppressWarnings(grep(x, t, ignore.case=TRUE))))
			   if(any(cols)) {
			   return(dat$species[suppressWarnings(grep(x, dat[,which(cols)], ignore.case=TRUE))])
			   } else {
			   return(list())
			   }
			   })
	res=unique(unlist(spp))
	return(res)
}

name.attribute.for.edges=function(edges, classification.df=NULL, clades=NULL){
## FIXME: 
	species=colnames(edges)
	if(is.null(attr(edges, "name"))) attr(edges, "name")=rep("",nrow(edges))
	dat=classification.df
	ss=which(names(dat)=="species")
	if(!length(ss)) stop("'classification.df' must minimally have a 'species' column'")
	dat=dat[dat$species%in%species,]
	higher_taxa=names(dat)[-ss]
	missing=c()
	
## FROM CLASSIFICATION.DF
	if(!is.null(classification.df)) {
		for(i in 1:length(higher_taxa)){
			ii=which(names(dat)==higher_taxa[i])
			taxa=unique(dat[,higher_taxa[i]])
			taxa=taxa[taxa!=""]
			for(taxon in taxa){
				spp=dat$species[grep(taxon, dat[,ii])]
				ww=which_edge(spp[spp%in%species], edges)
				if(!length(ww)) missing=c(missing, taxon)
				attr(edges, "name")[ww]=taxon
				print(taxon)
			}
		}
		
## FROM CLADES
		if(!is.null(clades)){
			for(i in 1:length(clades)){
				taxon=names(clades)[i]
				clade=clades[[taxon]]
				spp=species.in.taxon(dat, clade)
				ww=which_edge(spp[spp%in%species], edges)
				if(!length(ww)) missing=c(missing, taxon)
				print(taxon)
				attr(edges, "name")[ww]=taxon
			}	
		}
	}
	
## FROM TIPS
	tips=which(rowSums(edges)==1)
	attr(edges, "name")[tips]=sapply(tips, function(x) species[edges[x,]==1])
	if(length(missing)) warning(paste("the following appear non-monophyletic: \n\t", paste(missing, collapse="\n\t"), sep=""))
	return(edges)
}


## DEFUNCT: see hash.tip
get.hash <- function(tips, tips_reference){
	x=as.integer(tips_reference%in%tips)
	return(digest(x))
}

## DEFUNCT: unique_edges() see hashes.phylo()
unique_edges=function(phy){
	
# phy: most usefully a set of (perfectly overlapping in species labels) trees
# returns all unique edges in tree set, stores 'hash' attribute for fast lookups
# similar to what ape:::prop.part() accomplishes (for single trees)
	
###	SET function below ('hash_proppart' will use ape:::prop.part; 'hash_edges' will use phylo:::get.edges)
#	f=hash_proppart
	f=hash_edges
	
	if(class(phy)!="multiPhylo") {
		if(class(phy)=="phylo") {
			return(f(phy))
		} else {
			stop("check that 'phy' is 'multiPhylo' or 'phylo' object")
		}
	} 
	
## ENSURE THAT ALL TREES HAVE IDENTICAL TIP SET
	checktrees=function(multiphylo){
		tips=multiphylo[[1]]$tip.label
		ck=sapply(2:length(multiphylo), function(x) if(length(intersect(tips, multiphylo[[x]]$tip.label))==length(tips)) return(TRUE) else return(FALSE))
		if(all(ck)) return(TRUE) else return(FALSE)
	}
	
	ck=checktrees(phy)
	
	if(!ck) warning("mismatch in extents of sampled species found across trees")
	
	tips=phy[[1]]$tip.label
	for(i in 2:length(phy)){
		tips=intersect(tips, phy[[i]]$tip.label)
	}
	
	edges=f(phy[[1]], tips=tips)
	tips=colnames(edges)
	
	tickerFreq=ceiling(seq(1,length(phy),length=30))
	
	cat("|",rep(" ",9),toupper("phylogenies complete"),rep(" ",9),"|","\n")
	cat(". ")
	
	if(length(phy)>1){
		for(i in 2:length(phy)){
			cur_edges=f(phy[[i]],tips=tips)
			edges=uniquify_edges(edges, cur_edges)
			if(i%in%tickerFreq){
				for(i in 1:sum(i==tickerFreq)) cat(". ")
			}
		}
		cat("\n")
	}
	
## ensure no duplicates
	hash=attributes(edges)$hash
	if(any(bb<-duplicated(hash))){
		ww=which(bb)
		hash=hash[-ww]
		edges=edges[-ww,]
		attr(edges,"hash")=hash
	}
	return(edges)
}

## DEFUNCT: store_edges() see hashes.phylo()
store_edges=function(phy,edges){
	
# phy: a phylo object
# edges: a matrix for all edges, with the attribute ('hash'), unique keys for each edge
# records edge ID in phy$label (uniquely if edges is built from multiple trees, of which 'phy' is a set)
	
	nn=max(phy$edge)
	if(is.null(phy$label)) phy$label=numeric(nn) else return(phy)
	for(i in 1:nn){
		phy$label[i]=find_edge(i,phy,edges)
	}
	return(phy)
}

## DEFUNCT: hash_edges()  see hashes.phylo()
hash_edges=function(phy, tips=NULL){
	
# phy: a phylo object
# tips: the order of species desired (order is important in hashing)
# returns an edges matrix with 'hash' attribute
# same as hash_proppart() but faster for trees fewer than 4000 tips
	
	edges=get.edges(phy)
	if(!is.null(tips)) {
		mm=match(tips,colnames(edges))
		edges=edges[,mm]
	}
	hash=apply(edges, 1, function(x) digest(as.numeric(x)))
	attr(edges, "hash")=hash
	return(edges)
}

hash_proppart=function(phy, tips=NULL){
	
# phy: a phylo object
# tips: the order of species desired (order is important in hashing)
# returns an edges matrix with 'hash' attribute
	
	N=Ntip(phy)
	nn=nrow(phy$edge)
	pp=prop.part(phy)
	s=numeric(Ntip(phy))
	edges=matrix(0,N+length(pp),N)
	diag(edges[1:N,1:N])=1
	for(i in 1:length(pp)){
		stmp=s
		stmp[pp[[i]]]=1
		edges[N+i,]=stmp
	}
	if(!is.null(tips)){
		mm=match(tips, attributes(pp)$labels)
		edges=edges[,mm]
	} else {
		tips=phy$tip.label
	}
	dimnames(edges) = list(NULL, tips)
	hash=apply(edges, 1, function(x) digest(as.numeric(x)))
	attr(edges, "hash")=hash
	return(edges)
}

## DEFUNCT: uniquify_edges  see hashes.phylo()
uniquify_edges=function(x,y){
	
# x: PRIMARY edges matrix with 'hash' attribute
# y: SECONDARY edges matrix with 'hash' attribute
# adds edges in y that are missing in x to x
	
	if(!all(colnames(x)==colnames(y))) stop("columns appear to be mismatched between 'x' and 'y'")
	xh=attr(x,"hash")
	yh=attr(y,"hash")
	add=!yh%in%xh
	if(any(add)){
		hashes=c(xh,yh[add])
		edges=rbind(x,y[add,])
		attr(edges, "hash")=hashes
	} else {
		edges=x
	}
	return(edges)
}

## DEFUNCT: find_edge()  see hashes.phylo()
find_edge=function(node, phy, edges){
	root=Ntip(phy)+1
	if(node<root) {
		tt=phy$tip.label[node]
	} else {
		tt=phy$tip.label[get.descendants.of.node(node, phy, tips=TRUE)]
	}
	tt=tt[tt%in%colnames(edges)]
	
	res=which_edge(tt, edges)
	if(!length(res)) res=NA
	return(res)
	
}

## DEFUNCT: which_edge()  see hashes.phylo()
which_edge=function(tips, edges){
	
# tips: a binary representation of an edge (a named binary vector determining which species are subtended... e.g.
#
#		 Monodelphis_americana  Monodelphis_adusta	Monodelphis_kunsi 
#		 0                      1					1 
#	
#		... is the edge subtending Monodelphis_adusta and Monodelphis_kunsi 
#
# edges: a matrix for all edges, with the attribute ('hash'), unique keys for each edge
	
	if(is.null(attr(edges, "hash"))) stop("'edges' must have 'hash' attribute")
	
	species=colnames(edges)
	if(!all(tips%in%species)) stop("tips must all occur in 'colnames(edges)'")
	tmp=numeric(length(species))
	tmp[match(tips, species)]=1
	names(tmp)=species
	tips=tmp
	key=digest(as.numeric(tips))
	
	return(which(attr(edges,"hash")==key))
}

get.edges <- 
function(phy){
	n=Ntip(phy)
	nr=dim(phy$edge)[1]
	mn=max(phy$edge)
	out <- .Call("binary_edges", tree=list(
										   NTIP = as.integer(n),
										   ROOT = as.integer(n+1),
										   ENDOFCLADE = as.integer(nr),
										   ANC = as.integer(phy$edge[,1]),
										   DES = as.integer(phy$edge[,2]),
										   MAXNODE = as.integer(mn),
										   EDGES = as.integer(array(matrix(0, mn, n)))),
				 PACKAGE = "phylo")
	out=matrix(out$EDGES,nrow=mn,byrow=TRUE)
	dimnames(out)=list(NULL, phy$tip.label)
	return(out)
} 