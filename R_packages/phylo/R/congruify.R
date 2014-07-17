
## CONGRUIFICATION FUNCTIONS ##
check.multicore=function(){
	tmp=rownames(installed.packages())
	if("multicore"%in%tmp) {
		require(multicore)
		return(TRUE)
	} else {
		return(FALSE)
	}
}


.build_calibrations=function(dat, supertree, tol=0.1){
#	dat:
#	anc des       time                             hash
#	1   0 640 352.234677 3a4adb7cc0d4a51b9012dfb5615b3d71
#	2 640 641 243.269677 33757769ee61bde8dd5574ae35b47053
	
#	supertree: phylo tree with 'hash' object
	
	N=Ntip(supertree)
	backbone_times=dat
	supertree_hash=supertree$hash
	df=data.frame(MRCA=supertree_hash[(N+1):length(supertree_hash)], MaxAge=NA, MinAge=NA, taxonA=NA, taxonB=NA, valid=FALSE, stringsAsFactors=FALSE)
	for(i in 1:nrow(df)){
		if(!is.na(hash.idx<-df$MRCA[i])){
			if(hash.idx%in%backbone_times$hash){
				node.idx=i+N
				df[i,c("MaxAge","MinAge")]<-age.idx<-backbone_times$time[match(hash.idx, backbone_times$hash)]
				df[i,c("taxonA","taxonB")]<-taxa.idx<-spanning_taxa.phylo(supertree, node.idx)
				if(age.idx>tol & all(!is.na(taxa.idx))) df[i,"valid"]=TRUE
			}
		}
	}
	df=df[df$valid,]
	return(df)
}

congruify.phylo=function(backbone, megaphylogeny, taxonomy=NULL, tol=0.1, method="pathd8"){
#	backbone: a time-calibrated phylogeny with tip-labels that can be treated as an exemplar for clades in 'supertree'
#		-- e.g., tip.label in 'backbone' might be "Pinus" while in 'supertree' we might have "Pinus_cembra"
#		-- tips in 'backbone' can be contained in 'supertree' (FIXME: is this true?)
#	megaphylogeny: a rooted phylogeny that is to be time-scaled based on 'backbone'
#	taxonomy: linkage table between tipsets for 'backbone' and 'supertree'; if empty, one is attempted to be built by 'supertree' labels
#		-- if NULL, 'backbone' tips must correspond to tips in 'supertree'... e.g., A, B, C in 'backbone'; A_sp1, B_sp2, C_sp3 in 'supertree'
#		-- rownames of taxonomy must be tips in megaphylogeny
	
	## functions
	method=match.arg(toString(method), c("NA", "pathd8"))
	
	hashes.mapping <- function (phy, taxa, mapping){
	## GENERAL FUNCTION: find hash tag for every edge in tree (using phy$tip.label or predefined set of 'taxa')
	# returns list of hash tags from node 1:(Nnode(phy)+Ntip(phy))
	# taxa: set of species used to create hash keys
	# mapping: named list -- names are tips in 'phy'; elements are tips represented by each tip in 'phy' (and should also be present in 'taxa')
		if(is.null(taxa)) stop("Must supply 'tips'.")
		if(!all(names(mapping)%in%phy$tip.label)) stop("'mapping' must be named list with names in tip labels of 'phy'.")
		
		mapping=mapping[match(names(mapping), phy$tip.label)]
		descendants <- .compile_descendants(phy)
		hashes <- sapply(descendants, function(desc) hash.tip(unlist(mapping[desc]), taxa))
		empty=digest(integer(length(taxa)))
		hashes[hashes==empty]=NA
		phy$hash=hashes
		phy=.uniquify_hashes(phy)

		return(phy)
	}
	
	times.mapping=function(phy, taxa, mapping){
	#	mapping: named list -- names are tips in 'phy'; elements are tips represented by each tip in 'phy' (and should also be present in 'taxa')
	#	taxa: species that are represented by tips in 'backbone' 
		
	# find hash tags for overlapping_edges
		backbone=hashes.mapping(phy, taxa, mapping)
		dat=collect_times(backbone)
		dat$hash=backbone$hash[dat$des] ## key to times and edges
		return(dat)
	}
	
	smooth_supertree=function(backbone, supertree, taxa, spp, tol=0.01, method=c("pathd8", NA)){
		method=match.arg(toString(method), c("NA", "pathd8"))
		if(!is.ultrametric(backbone)) warning("Supplied 'backbone' is non-ultrametric.")
		backbone_dat=times.mapping(backbone, taxa, spp)
		calibration=.build_calibrations(backbone_dat, supertree, tol=tol)
		if(!nrow(calibration)) {
			warning("No concordant branches reconciled between 'backbone' and 'megaphylogeny'")
			return(NA)
		}
		if(method=="pathd8") {
			phy=smooth.phylo(supertree, calibration, base=".tmp_PATHd8", rm=FALSE)
			phy$hash=c(rep("", Ntip(phy)), phy$node.label)
			phy$hash[phy$hash==""]=NA
		} else if(method=="NA"){
			phy=NULL
		}
		backbone$hash=backbone_dat$hash[match(1:(Ntip(backbone)+Nnode(backbone)), backbone_dat$des)]
		backbone$node.label=backbone$hash[(Ntip(backbone)+1):max(backbone$edge)]
		backbone$node.label[is.na(backbone$node.label)]=""
		return(list(phy=phy, cal=calibration, backbone=backbone))
	}
	
	## end functions
	
	## PROCESSING ##
	supertree=megaphylogeny
	classification=taxonomy

	if(class(backbone)=="phylo") backbone=list(backbone)
	if(is.null(classification)) classification=.build_classification(supertree$tip.label)
	
	tips=unique(unlist(lapply(backbone, function(x) x$tip.label)))
	spp=lapply(tips, function(o) {
			   x=rownames(classification)[which(classification==o, arr=TRUE)[,1]]
			   x=x[x%in%supertree$tip.label]
	})
	names(spp)=tips
	taxa=unique(unlist(spp))
	
#	supertree$hash=find_supertree_edges(supertree, taxa)
	supertree=hashes.phylo(supertree, taxa)
	
	f=lapply
	results=f(1:length(backbone), function(i) {phy=backbone[[i]]; smooth_supertree(phy, supertree, taxa, spp, tol=tol, method=method)})
	return(results)
}

## END CONGRUIFICATION FUNCTIONS ##


## CHECKING FUNCTIONS ##

tips.hash=function(phy, hash){
# phy needs 'hash' object 
	idx=which(phy$hash==hash)
	desc=get.descendants.of.node(idx, phy, tips=TRUE)
	if(idx<=Ntip(phy)){
		return(phy$tip.label[idx])	
	} else {
		return(sort(phy$tip.label[desc]))	
	}
}

time.hash=function(phy, hash){
# phy needs 'hash' object
	idx=which(phy$hash==hash)
	desc=get.descendants.of.node(idx, phy, tips=TRUE)
	if(idx<=Ntip(phy)){
		return(phy$tip.label[idx])	
	} else {
		return(sort(phy$tip.label[desc]))	
	}
}

matching.congruified=function(supertree, backbone, hash){
# supertree and backbone need 'hash' function
	s=tips.hash(supertree, hash)
	b=tips.hash(backbone, hash)
	stt=branching.times(supertree)
	st=stt[[which(names(stt)==hash)]]
	btt=branching.times(backbone)
	bt=btt[[which(names(btt)==match(hash, backbone$hash))]]
	return(list(supertree=list(tips=s, time=st), backbone=list(tips=b, time=bt)))	
}

find_splittime=function(tips, phy){
# finds split time of minimal node subtending tips
	anc=getMRCA(phy, sapply(tips, function(x) match(x, phy$tip.label)))
	bt=branching.times(phy)
	bt[[names(bt)==anc]]
}


fetch_hash=function(tips, label, taxonomy, rank="species"){
	tt=taxonomy[which(taxonomy==label, arr=TRUE)[,1],rank]
	hash.tip(labels=tt, tips=tips)
}

## END CHECKING FUNCTIONS -- probably not super useful ##

