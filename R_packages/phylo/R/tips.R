tips_plot.phylo=function(phy, tips, col, ...){
	plot(phy, edge.width=0.1, type="f", ...)
	tiplabels(frame="n", pch=19, col=ifelse(phy$tip.label%in%tips, col, NA), cex=0.1)
}

split_taxon=function(label){
	dat=unlist(strsplit(label, " "))
	if(length(dat)>=2){
		if(length(dat)>2){
			s=sapply(list("var.", "subsp.", "f."), function(t) {g=dat==t; if(any(g)) return(max(which(g))) else return(NULL)})
			if(length(unlist(s)->m)){
				m=max(m)
				if(length(dat)>m+1){
					spp=paste(dat[1:(m+1)], collapse="_")
					authority=paste(dat[(m+2):length(dat)], collapse=" ")
				} else if(length(dat)==m+1){
					spp=paste(dat[1:(m+1)], collapse="_")
					authority=NA
				} else {
					stop("Complex species name uninterpretable.")
				}
			} else {
				spp=paste(dat[1:2], collapse="_")
				authority=paste(dat[3:length(dat)], collapse=" ")
			}
		} else {
			spp=paste(dat[1:2], collapse="_")
			authority=NA
		}
		
	} else {
		stop("Uninterpretable species name.")
	}
	return(c(spp, authority))
}

binomial_lab=function(label){
	tmp=unlist(strsplit(label, "_", fixed=TRUE))
	if(length(tmp)>=2) return(paste(tmp[1:2], collapse="_")) else return(NULL)
}


is.ambig=function(label){
	s=split_lab(label)
	if(any(s=="sp" | s=="sp.")) return(TRUE) else return(FALSE)
}

is.aff=function(label){
	s=split_lab(label)
	if(any(s=="cf" | s=="aff")) return(TRUE) else return(FALSE)
}

is.unisexual=function(label){
	s=split_lab(label)
	if(any(s=="unisexual")) return(TRUE) else return(FALSE)
}

split_lab=function(label){
	lab=gsub(".", "_", label, fixed=TRUE)
	lab=gsub(" ", "_", lab, fixed=TRUE)
	lab=unlist(strsplit(lab, "_", fixed=TRUE))
	lab=lab[lab!=""]
	lab
}

# grab two taxa that span the given node
spanning_taxa.phylo=function(phy, node){
	dd=get.desc.of.node(node, phy)
	if(length(dd)) {
		desc=phy$tip.label[sapply(dd[1:2], function(x) {y=get.descendants.of.node(x, phy, tips=TRUE); return(min(x, y[1]))})]
	} else {
		desc=NA
	}
	return(desc)
}	

