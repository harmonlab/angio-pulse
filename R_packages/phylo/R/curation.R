curate.labels=function(labels, taxonomy, novon=NULL, synonyms=NULL){
#labels: a data.frame with label="", match=NA, usable=NA, guess=NA, newlabel=NA, redundant=NA
#taxonomy: a taxonomic database
#novon: a vector of valid taxa not in sys$species 
#synonyms: a named list of synonyms
	
# PRIMARILY used with data from amphibiaweb.org
	sys=taxonomy
	spp=data.frame(label=labels, match=NA, usable=NA, guess=NA, newlabel=NA, redundant=NA) 

	## curation ##
	for(i in 1:nrow(spp)){
		curspp=spp$label[i]
		
	# novon
		x=curspp%in%novon
		if(x) {
			spp$match[i]="novon"
			spp$usable[i]=TRUE
			spp$newlabel[i]=curspp
			next()
		}
		
	# synonomy
		if(!is.null(synonyms)){
			x=sapply(synonyms, function(s) if(curspp%in%s) return(TRUE) else return(FALSE))
			if(any(x)) {
				spp$match[i]="synonym"
				spp$usable[i]=TRUE
				if(all(names(x)[x]%in%c(novon, sys$species))) {
					spp$newlabel[i]=paste(names(x)[x], collapse="_AND_")
				} else {
					print(paste("synonym: ", curspp,  "; ", i, sep=""))
					stop(paste(paste(names(x)[x], collapse="_AND_"), "not valid", sep=" "))
				}
				next()
			}
		}
		
		
	# binomial
		if(curspp%in%c(novon, sys$species)) {
			spp$match[i]="binomial"
			spp$usable[i]=TRUE
			spp$newlabel[i]=curspp
			next()
		}
		
	# unisexual
		x=is.unisexual(curspp)
		if(any(x)) {
			spp$match[i]="unisexual"
			spp$usable[i]=FALSE
			next()
		}
		
	# near match
		x=split_lab(curspp)
		if(length(x>2)){
			gs=paste(x[1:2], collapse="_")
			if(length(unlist(sapply(x[-c(1,2)], grep, sys$species)))==0 & gs%in%c(novon, sys$species)){
				spp$match[i]="nearbinomial"
				spp$usable[i]=TRUE
				spp$newlabel[i]=gs
				print(c(curspp, gs))
				next()
			}
		}
		
	# ambiguous
		x=is.ambig(curspp)
		if(x) {
			spp$match[i]="ambiguous"
			spp$usable[i]=TRUE
			new=paste(split_lab(curspp)[1:2], collapse="_")
			new=gsub(".", "", new, fixed=TRUE)
			spp$newlabel[i]=new
			print(c(curspp, new))
			next()
		}
		
	# affinity
		x=is.aff(curspp)
		if(x) {
			spp$match[i]="affinity"
			spp$usable[i]=TRUE
			new=paste(split_lab(curspp)[c(1,3)], collapse="_")
			if(new%in%c(novon, sys$species)){
				spp$newlabel[i]=new
			} else {
				print(paste("affinity: ", curspp, "; ", i, sep=""))
				stop(paste(new, "not valid", sep=" "))
			}
			next()
		}
				
	# subspecies
		x=paste(split_lab(curspp)[1:2], collapse="_")
		if(x%in%sys$species){
			spp$match[i]="subspecies"
			spp$usable[i]=TRUE
			spp$newlabel[i]=x
			print(c(curspp,x))
			next()
		} 
		
	# fuzzy
		dd=sapply(sys$species, function(x) stringDist(c(curspp, x)))
		ww=which(dd==min(dd))
		if(length(ww)==1) spp$match[i]="fuzzy" else spp$match[i]="multiple"
		spp$guess[i]<-spp$newlabel[i]<-tmp<-paste(names(ww), collapse="_AND_")
		print(c(curspp,tmp))
		next()
		
	}
	
	# adjust systematics database
	add=lapply(novon, retrieve.taxonomy, sys)
	names(add)=names(novon)
	for(i in 1:length(add)){
		sys=rbind(sys, c(add[[i]], novon[i]))
	}
	
	# find missing taxa
	missing=find.missing(spp$newlabel, sys)
	
	return(list(spp=spp, sys=sys, missing=missing))
}

retrieve.taxonomy=function(lab, taxonomy){
	x=split_lab(lab)
	df=taxonomy[taxonomy$genus==x[1],]
	nn=names(taxonomy)
	nn=nn[nn!="species"]
	return(unlist(unique(df[,nn])))
}

acceptabletips.phylo=function(phy, taxonomy, rank="species"){
# converts tip labels and returns tree with only "acceptable" names (those found in taxonomy$rank)
# will further convert species labels: from 'A_b_var._d' to 'A_b'
	checktip=function(tip, taxonomy, rank){
		tip=binomial_lab(tip)
		if(!is.null(tip)){
			if(tip%in%taxonomy[,rank]){
				return(tip)
			} else {
				return("NULL")
			}
		} else {
			return("NULL")
		}
	}
	tips.bv=phy$tip.label%in%taxonomy[,rank]
	phy$tip.label[!tips.bv]=sapply(phy$tip.label[!tips.bv], function(t) checktip(t, taxonomy, rank))
	phy=drop.tip(phy, "NULL")
	return(subset.phylo(phy))
}
