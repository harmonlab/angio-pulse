## se_value=0.0345 ## SEE HARMON et al. 2010 Evol
prune_data.polytomy=function(phy, dat, se_value=0.0345){
	tmp=polytomy_sampler(phy)
	ptre=tmp$phy
	mapping=tmp$mapping
	
	composite_data=lapply(names(mapping), function(l) create_composite(l, dat, mapping, se_value))
	names(composite_data)=names(mapping)
	sterrors=sapply(composite_data, "[[", "se")
	means=sapply(composite_data, "[[", "m")
	
	
	res=dropdata(ptre, means, sterrors)
	return(list(phy=res$phy, dat=res$dat, se=res$se, mapping=mapping))
}

create_composite=function(label, dat, mapping, se_value){
	poss=unique(c(label, mapping[[label]]))
	val=dat[poss]
	val=val[!is.na(val)] 
	if(length(val)){
		se=ifelse(length(val)==1, se_value, sem(val))
		return(c(m=mean(val), se=se))
	} else {
		return(c(m=NA, se=NA))
	}
}

dropdata=function(phy, means, sterrors){
	drop.idx=is.na(means)
	m=means[!drop.idx]
	se=sterrors[!drop.idx]
	phy=drop.tip(phy, phy$tip.label[!phy$tip.label%in%names(m)])
	return(list(phy=phy, dat=m, se=se))
}

match_data.polytomy=function(phy, dat, mapping, se_value=0.0345){
## mapping: a named list -- an 'exemplar' mapping between the name (a tip in 'phy') and the tips for which it is representing (due to a polytomy often)
#	finds means for clades represented by an exemplar and prunes tree to match the data
#	treedata() may still be required if there are some names in 'mapping' that cannot be matched to the tree
	
	composite_data=lapply(names(mapping), function(l) create_composite(l, dat, mapping, se_value))
	names(composite_data)=names(mapping)
	sterrors=sapply(composite_data, "[[", "se")
	means=sapply(composite_data, "[[", "m")
	
	
	res=dropdata(phy, means, sterrors)
	return(list(phy=res$phy, dat=res$dat, se=res$se, mapping=mapping))
}
