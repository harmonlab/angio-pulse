## PROMINENCE OF HASHES (rather than edges)
#	-- rarely are the edges themselves (the bit strings denoting which species are subtended) are used
#	-- rely (and speed up) get.hashes()

## NEEDS
#	-- identify the full set of hashes (edges) across a distribution of trees (use a particular ordering and a set of 'tips' that are consistent among trees)
#	-- identify hash key for a set of tips (or single edge in tree)
#	-- deal with edges that are non-unique within a tree (e.g., 'tips' is incomplete and thus multiple branches subtend the same subset of 'tips')
#			- deal with edges that are empty (e.g., no species in 'tips' subtended)
#	-- deal with cases where levels of sampling differ among trees (e.g., one tree is species level, another is family level but we have a 'lookup' for species-to-family)
# uniquify redundant hash keys
.uniquify_hashes = function(phy) {
	
	hash=phy$hash
	if(is.null(hash)) stop("Must supply 'phy' with 'hash' object.")
	
	tt = table(hash[!is.na(hash)])
	if (any(tt > 1)) { # find each tipmost hash to store as unique
		N=Ntip(phy)
		subset = names(tt[tt > 1])
		for (s in subset) {
			idx = which((hash == s) == TRUE)
			if (s %in% hash[1:N]) {
				hash[idx[idx != min(idx)]] = NA
			}
			else {
				hash[idx[idx != max(idx)]] = NA
			}
		}
		phy$hash=hash
	}
	return(phy)
}


hashes.phylo<-
function(phy, tips=NULL){
	## FIXME: distinguish between an empty edge (one that subtends node of the tips in 'tips' and a redundant edge (one that subtends same set of tips in 'tips' as another edge)? 
	# 
	# tips: an ordering (and set) from which to determine unique hashes
	
	if( !(class(phy)%in%c("phylo")->is.single) & !(class(phy)%in%"multiPhylo"->is.multi) ) stop("Supply 'phy' as a 'phylo' or 'multiPhylo' object.")
	
	# returns 'phy' with 'hashtips' attribute and 'hash' object
	## USE edges.phylo() if edges matrix is needed
	
	## REPLACES: 
	#		get.hashes() (returns 'hashes' with 'tips' attribute)
	#		store_edges() (returns 'phy' with 'label' object 
	
	## near REDUNDANCIES: 
	#		edges.phylo()
	
	if(is.single) trees=list(phy) else trees=phy
	
	if(is.null(tips)) {
		tips_tmp=trees[[1]]$tip.label
		if(is.multi){
			for(i in 2:length(trees)){
				tips_tmp=intersect(tips_tmp, trees[[i]]$tip.label)
			}
		} 
		tips=tips_tmp
	} 
	
	NULL_TIPS=digest(integer(length(tips)))
	if(check.multicore()) f=f=function(X,FUN) mclapply(X,FUN,mc.silent=TRUE) else f=lapply
	
	# store hash keys for trees in order 1:(Node(phy)+Ntip(phy))
	res=lapply(trees, function(phy){
		storage.mode(phy$Nnode) <- "integer"
		# descendants.idx <- .Call("prop_part", list(phy), 1, TRUE, PACKAGE = "ape")
		descendants.idx <- .compile_descendants(phy)
		hashes <- unlist(f(descendants.idx, function(desc) hash.tip(phy$tip.label[desc], tips)))
#		hash_tips <- as.character(unlist(f(phy$tip.label, function(tip) hash.tip(tip, tips))))
#		hashes=c(hash_tips, hash_internal)
		hashes[hashes==NULL_TIPS]=NA
		phy$hash=hashes
		attr(phy,"hashtips")=tips
		return(.uniquify_hashes(phy))
	})
	
	if(is.single) {
		return(res[[1]]) 
	} else {
		class(res)="multiPhylo"
		return(res)
	}
}
	

hash.tip <- function(labels, tips){
	x=as.integer(tips%in%labels)
	return(digest(x))
}

hash.node <- function(node, phy, tips=NULL){
	if(is.null(tips)) tips=phy$tip.label
	if(node<=Ntip(phy)){
		descendants=phy$tip.label[node]
	} else {
		descendants=phy$tip.label[get.descendants.of.node(node, phy, tips=TRUE)]		
	}
	x=as.integer(tips%in%descendants)
	return(digest(x))
}

## SELDOM (NEVER) used? (only if requiring binary matrix of edges)
edges.phylo<-
function(phy, tips=NULL, package=c("phylo", "ape")){
	
	## USE hashes.phylo() if edges matrix is not needed
	## WORKS on single trees 
	
	## REPLACES: 
	#		get.edges()
	#		hash_proppart()
	
	.hash_edges <- 
	function(phy, tips=NULL){
		
	# phy: a phylo object
	# tips: the order of species desired (order is important in hashing)
	# returns an edges matrix with 'hash' attribute
	# same as hash_edges() but slower for trees fewer than 4000 tips
		
		
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
		if(!is.null(tips)){
			mm=match(tips, phy$tip.label)
			out=out[,mm]
		} else {
			tips=phy$tip.label
		}
		out[is.na(out)]=NA
		dimnames(out) = list(NULL, tips)
		hash=apply(out, 1, function(x) digest(as.integer(x)))
		attr(out, "hash")=hash
		return(out)
	} 
	
	.hash_proppart<-
	function(phy, tips=NULL){
		
	# phy: a phylo object
	# tips: the order of species desired (order is important in hashing)
	# returns an edges matrix with 'hash' attribute
	# same as hash_edges() but slower for trees fewer than 4000 tips
		
		
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
		edges[is.na(edges)]=0
		hash=apply(edges, 1, function(x) digest(as.integer(x)))
		attr(edges, "hash")=hash
		return(edges)
	}
	
	
	f=match.arg(package, c("phylo","ape"))
	res=switch(f, "phylo"=.hash_edges(phy, tips=tips), "ape"=.hash_proppart(phy, tips=tips))
	return(res)
}



