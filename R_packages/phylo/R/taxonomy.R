resolve_taxonomy=function(tips, rank="family"){
	orig_rank=rank
	rank="family"
	
## family-level resolution (ONLY at this point)
	tax=get(data(spermatophyta_APGIII_taxonomy))
	tmp=fetch_gbtaxonomy.above("foo") ## prompt construction of taxdump if necessary (before invoking multicore)
	
	gtips=sapply(tips, function(x) unlist(strsplit(gsub(" ", "_", x), "_"))[1])
	check=sapply(list(tips, gtips), function(x) sum(x%in%tax[,1]))
	if(max(check)==check[2] & length(unique(check))==2) {
		tips=gtips
		warning("Assuming 'tips' supplied as 'Genus_species'")
	} 
	tips=unique(tips)
	
## WIKISPECIES (APGIII) TAXONOMY
	if(check.multicore()) f=mclapply else f=lapply
	if(rank!="family") {
		warning(paste(rank, "search in 'wikispecies' not yet implemented", sep=" "))
		wiki=character(length(tips))
		wiki[]=NA
	} else {
		
		trim_down_wiki=function(result){
			if(length(result)>1) return(NA)
			xx=gsub("(", "SPLIT___SPLIT", gsub(")", "SPLIT___SPLIT", result, fixed=TRUE), fixed=TRUE)
			ss=parse_wiki(xx, "aceae")
			ss
		}
		
		wiki=f(tips, function(x) {
#			   cat(x,"\n")
			   y=try(browse_wikispecies.family(x, open=FALSE), silent=TRUE)
			   if(is.null(y) | inherits(y,"try-error")) return(NA) else return(trim_down_wiki(y))
			   })
	}
	
	
## GENBANK TAXONOMY
	if(check.multicore()) f=mclapply else f=lapply
	gb=f(tips, function(x) {
#		 cat(x,"\n")
		 y=try(fetch_gbtaxonomy.above(x, rank)$all[[rank]],silent=TRUE)
		 if(is.null(y) | inherits(y,"try-error")) {
			return(NA) 
		 } else {
			if(length(grep("idae",y))) return(NA) else return(y)
		 }
	})
	
## PLANTLIST TAXONOMY
	plantlist=get(data(plantlist_classification))
	if(check.multicore()) f=mclapply else f=lapply
	pl=f(1:length(tips), function(i){
		 tip=tips[i]
		 ww=which(tip==plantlist, arr.ind=TRUE)[,"row"]
		 if(length(ww)){
		 ff=unique(plantlist[ww,rank])
		 if(length(ff)==1){
		 return(ff)
		 } else {
		 warning(paste("Problem in 'plantlist' detected: multiple taxa encountered for ", tip, ":\n\t", paste(ff, collapse="\n\t"), sep=""))
		 return(NA)
		 }
		 } else {
		 return(NA)
		 }
		 })
	pl=unlist(pl)
	
## store output
	mm=cbind(tip=tips, NESCENT_wg=NA, plantlist=pl, gb=unlist(gb), wiki=unlist(wiki), notes=NA)
	for(i in 1:nrow(mm)){
		
		tip=tips[i]
		check=which(tax==tip,arr.ind=TRUE)[,"row"]
		
## CHECK NESCENT_wg
		if(length(check)){
			ff=unique(tax[check,rank])
			if(!is.na(ff) & length(ff)==1) {
				mm[i,"NESCENT_wg"]=ff 
			}
		}
	}
	
## CHECK SYNONYMY
	swap_synonym=function(taxon, synonyms){
		if(length(taxon)>1) stop("supply a single 'taxon'")
		tt=which(taxon==synonyms, arr.ind=TRUE)
		if(nrow(tt)==1){
			swap=synonyms[tt[,"row"],2] ## point to 'to'
			return(list(taxon=swap, replaced=!(swap==taxon)))
		} else {
			return(list(taxon=taxon, replaced=FALSE))
		}
	}
	
## SWAP SYNONYMS
	if(rank=="family"){
		family_synonyms=get(data(spermatophyta_APGIII_families))$synonymy
		uu=unique(unlist(family_synonyms))
		for(j in c("NESCENT_wg","gb","plantlist","wiki")){
			cur=mm[,j]
			zz=match(cur,uu)
			
			swaps=lapply(1:length(zz), function(idx) if(!is.na(idx)) return(swap_synonym(cur[idx], family_synonyms)) else return(NA))
			swapped=unname(sapply(swaps, "[[", "replaced"))
			mm[,j]=unname(sapply(swaps, "[[", "taxon"))
			
# document change
			message=ifelse(swapped, paste(mm[,"notes"], paste("synonym",j,sep="_"), sep="_"), mm[,"notes"])
			mm[,"notes"]=message
		}
	}
	
	res=data.frame(mm, stringsAsFactors=FALSE)
	
## check agreement between 'genbank' and 'plantlist'
	res$agreement=0
	res$sum=0
	res$all=0
	for(i in 1:nrow(res)){
		xx=unlist(res[i,c("plantlist","gb","wiki")])
		if(length(unique(xx))==1 & !all(is.na(xx))) {
			if(length(unique(unlist(c(xx, unlist(res[i,"NESCENT_wg"])))))==1) res$all[i]=1
			res$agreement[i]=1
		}
		res$sum[i]=length(xx[!is.na(xx)])
	}
	
## EXPLAIN columns
	cat("'all' is conditional for all labels matching\n'agreement' is conditional for 'gb' and 'plantlist' matching\n'sum' is the amount of data from 'gb', 'plantlist', and 'wiki'\n") 
	
	res=res[order(res$agreement, res$sum),]
	res$mismatch=0
	res$family=NA
	
## DECIDE FAMILY based on corroborating evidence
	for(idx in 1:nrow(res)){
		guesses=unique(unlist(res[idx,c("NESCENT_wg","plantlist","gb","wiki")]))
		guesses_nonNA=guesses[!is.na(guesses)]
		
		
		if(length(guesses_nonNA)==1){ ## use if only one unique family
			res$family[idx]=guesses_nonNA
		} else { ## use if corroboration from external sources
			if(res$sum[idx]>1){
				external=unique(unlist(res[idx,c("plantlist","gb","wiki")]))
				external=external[!is.na(external)]
				if(length(external)) {
					if(length(external)==1){
						res$family[idx]=external
					}
					nescent_fam=res$NESCENT_wg[idx]
					if(!is.na(nescent_fam)){
						if(!nescent_fam %in% external) res$mismatch[idx]=1
					}
				} 
			} 
		}
	}
	
	if(orig_rank=="order"){
		ord=get(data(spermatophyta_APGIII_orders))
		tmp=unlist(sapply(1:length(ord), function(idx) {y=ord[[idx]]; names(y)=rep(names(ord)[idx], length(y)); return(y)}))
		ord=data.frame(family=tmp, order=names(tmp), stringsAsFactors=FALSE)	
		res$order=ord$order[match(res$family, ord$family)]
	}
	
	
	
	res
}


check.taxonomy=function(taxon, position, taxonomy){
# determines whether each nested rank is exclusively nested within a single higher taxon 
	hits=taxonomy[which(taxonomy[,position]==taxon),position-1]
	if(length(unique(hits))==1) return(unique(hits)) else return(NA)	
}

.build_classification=function(species){
	data.frame(genus=sapply(species, function(s) split_lab(s)[1]), species=species, stringsAsFactors=FALSE)
}


fill.taxonomy=function(taxonomy){
	push.taxon=function(taxonomy, indices, column){
		for(n in which(indices)){
			taxonomy[n,column]=paste(taxonomy[n,min((column:ncol(taxonomy))[!is.na(taxonomy[n,column:ncol(taxonomy)])])],paste("R",column,sep=""),sep=".")
		}
		return(taxonomy)
	}
	
	for(c in 1:(ncol(taxonomy)-1)){
		ii=is.na(taxonomy[,c])
		if(any(ii)) taxonomy=push.taxon(taxonomy, ii, c)
	}
	return(taxonomy)	
}


build_phylo.lookup=function(lookup, tips=rownames(lookup), rank=NULL) {
# GENERAL FUNCTION
# lookup is data.frame with ranks as columns
# tips is vector of strings (tip labels), ordered as rows in 'lookup'
# rank corresponds to (numeric) position in rank names (R to L) to use in building taxonomic tree; if null, all ranks are used
# NOTE: taxonomic groups in 'lookup' MUST be ordered by exclusivity from R to L -- e.g., genera must precede families must precede orders
	labels=tips
	taxonomy=lookup
	if(class(taxonomy)=="matrix") taxonomy=data.frame(taxonomy, stringsAsFactors=FALSE)
	taxonomy[taxonomy==""]=NA
	taxonomy=taxonomy[,ncol(taxonomy):1]
	
	# check for 'species' column
	check=apply(taxonomy, 2, function(x) {y=all(x==labels); if(is.na(y)) return(FALSE) else return(y)})
	if(any(check)) taxonomy=as.data.frame(taxonomy[,-which(names(taxonomy)==names(which(check)))])
		
	# check for invariant columns
	check=apply(taxonomy, 2, function(x) length(unique(x)))
	if(any(check==1)) taxonomy=as.data.frame(taxonomy[,-which(names(taxonomy)==names(which(check==1)))])
	
	# initialize taxonomy for missing data
	if(is.null(rank)) rank=1:ncol(taxonomy)

	# find rows lacking information
	t=apply(taxonomy,1,function(t)any(is.na(t)))
	if(any(t)) {
		missing=labels[t]
	} else {
		missing=NULL
	}
	
	# fill in missing information
	cat("Filling taxonomy and checking for internal consistency\n")
	tt=fill.taxonomy(cbind(taxonomy, labels))
	taxonomy=data.frame(tt[,1:(ncol(tt)-1)])
	
	# prepare dataframe for processing
	if(class(taxonomy)=="data.frame" | class(taxonomy)=="matrix"){
		taxonomy=as.matrix(taxonomy)
		if(ncol(taxonomy)>=max(rank)) {
			taxonomy=as.matrix(taxonomy[,rank])
			tt=apply(taxonomy,2,function(x)all(x==labels))
			if(any(tt)) taxonomy=taxonomy[,-which(tt)]
			rank=seq(1,ncol(taxonomy))
			taxonomy=cbind(taxonomy,labels)
			colnames(taxonomy)=c(paste("R",rank,sep=""),"labels")
		} else {
			stop("Supplied taxonomy appears inconsistent with called rank(s).")
		}
	} else {
		stop("Please supply taxonomy as a data.frame.")
	}
	
	### BUILD EDGE MATRIX ###
	cat("Constructing tree representation of taxonomy\n")

	# find unique rank labels
	uranks=lapply(rank, function(x) unique(taxonomy[,x]))
	
	n=nrow(taxonomy)
	N=length(unlist(uranks))+1
	
	# initialize internal and terminal edge matrices
	int.edge=matrix(NA, nrow=N-1, ncol=2)
	tip.edge=matrix(NA, nrow=n, ncol=2)
	
	
	## LABEL INTERNAL NODES ##
	start=n+2
	start.node=0
	for(i in 1:length(rank)) {
		
	# numerically label each unique rank
		curints=length(uranks[[i]])
		end=start+curints-1
		names(uranks[[i]])=start:end
		edge.start=start-n-1
		start=end+1
		
	# indexing for edge matrix
		edge.end=edge.start+curints-1
		
		if(i==1) {
			
	# first rank is bounded by root (n+1)
			int.edge[edge.start:edge.end,]=c(rep(n+1,curints),as.numeric(names(uranks[[i]])))
			
		} else {
			
	# subsequent ranks must be linked to their ancestral node ID
			monophyletic.ranks=sapply(uranks[[i]], function(x) check.taxonomy(x, rank[i], taxonomy))
			
	# record ancestor-descendant relationships for all internal nodes --> int.edge[ancestor,descendant]
			if(!any(is.na(monophyletic.ranks))) {
				start.row<-start.node<-start.node+length(uranks[[i-1]])
				int.nodes=as.numeric(names(uranks[[i-1]][match(monophyletic.ranks, uranks[[i-1]])]))
				int.edge[(start.row+1):(start.row+length(uranks[[i]])),]=cbind(int.nodes, as.numeric(names(uranks[[i]])))
			} else {
				stop(paste("\n\nNOTE: ensure that 'taxonomy' is supplied with most exclusive taxa (e.g., species) in the leftmost columns\n\nThe following taxa appear non-monophyletic: \n\t"),paste(uranks[[i]][which(is.na(monophyletic.ranks))],collapse="\n\t"))
			}
		}
	}
	
	## LABEL TERMINAL NODES ##
	innermost.rank=uranks[[length(rank)]]
	innermost.rank.IDs=as.numeric(names(innermost.rank))
	innernodes=unname(sapply(taxonomy[,max(rank)], function(x) innermost.rank.IDs[which(innermost.rank==x)]))
	tip.edge=cbind(innernodes, 1:n)
	
	# build phylo object
	edge=rbind(int.edge,tip.edge)
	phy=list(edge=edge, 
			 Nnode=N, 
			 tip.label=taxonomy[,ncol(taxonomy)],
			 edge.length=rep(1,nrow(edge)),
			 root.edge=0
			 )
	class(phy)="phylo"
	phy=reorder(collapse.singles(phy))		
	phy$tip.label=unname(phy$tip.label)
	return(list(phy=phy, missing=missing))		
}

### AMPHIBIAWEB ###
fetch.amphibiaweb=function(species=""){
	faw=function(spp){
		label=split_lab(spp)
		if(length(label)>2) warning("Encountered more than the generic and specific epithets.")
		if(length(label)<2) stop("Insufficient information: 'species' should be given as 'Genus_species'.")
		cat(paste("Fetching data for: ", spp, "\n", sep=""))
		url=paste("\"http://amphibiaweb.org/cgi-bin/amphib_query?query_src=aw_search_index&table=amphib&special=one_record&where-genus=", label[1], "&where-species=", label[2], "\"", sep="")
		tmp=".tmpAmphibiaWeb"
		system(paste("curl", url, ">", tmp, sep=" "), wait=TRUE, ignore.stderr=TRUE)
		dat=readLines(tmp)
		if(length(grep("Sorry - no matches.", dat))) dat=c("")
		unlink(tmp)
		dat
	}
	dat=lapply(species, function(x) {t=try(faw(x)); if(inherits(t, "try-error")) t=c(""); t})
	failed=function(html) return(ifelse(length(grep("Sorry, our database can't handle that many queries", html)), TRUE, FALSE))
	
	failures=sapply(dat, failed)
	
	dat=lapply(dat, filter.html)
	names(dat)=species
	return(list(dat=dat, failed=failed))
}

### AMPHIBIAWEB ###
search.amphibiaweb=function(species=""){
	label=split_lab(species)
	if(length(label)>2) warning("Encountered more than the generic and specific epithets.")
	if(length(label)<2) stop("Insufficient information: 'species' should be given as 'Genus_species'.")
	url=paste("\"http://amphibiaweb.org/cgi-bin/amphib_query?query_src=aw_search_index&table=amphib&special=one_record&where-genus=", label[1], "&where-species=", label[2], "\"", sep="")
	system(paste("open", url, sep=" "))
}

### PLANTLIST.ORG ###
split_taxa=function(dat){
	dat$Species=NA
	dat$Authority=NA
	if(all(c("Name", "Status", "Source")%in%names(dat))){
		for(i in 1:nrow(dat)){
			tmp=split_taxon(dat$Name[i])
			dat$Species[i]=tmp[1]
			dat$Authority[i]=tmp[2]
		}
		return(dat[,c("Species", "Authority", "Status", "Source")])
	} else {
		return(NULL)
	}
}

### WIKISPECIES ###
parse_wiki=function(string, type="aceae"){
	yy=unlist(strsplit(string, "SPLIT"))
	if(length(yy)==1 & length(grep(type,yy))) return(yy)
	ww=which(yy=="___")
	ff=grep(type,yy)
	if(length(ff)){
		for(i in 1:length(ff)){
			idx=ff[i]
			if(all(c(idx-1,idx+1)%in%ww)) { ## split on either side of 'aceae' match
				return(yy[idx])
			} 
		}
		return(NA)
	} else {
		return(NA)
	}
}

### WIKISPECIES ###
browse_wikispecies.family=function(taxon, open=FALSE){
	
	## get genus if supplied species
	taxon=split_lab(taxon)[1]
	browse=TRUE
	query_by_browse=browse
	
	pop=function(url) system("open -a Safari",url,intern=FALSE, wait=TRUE, ignore.stderr=TRUE, ignore.stdout=TRUE)
	
	split_wiki=function(string){
		ss=gsub("\"", "SPLIT___SPLIT", gsub("href=\"/wiki/", "SPLIT___SPLIT", string, fixed=TRUE), fixed=TRUE)
		xx=parse_wiki(ss, type="aceae")
		if(!is.na(xx)) {
			yy=gsub("(", "SPLIT___SPLIT", gsub(")", "SPLIT___SPLIT", xx, fixed=TRUE), fixed=TRUE)
			ss=parse_wiki(yy, "aceae")
			return(ss)
		} else {
			return(xx)
		}
	}
	
	## NOT WORKING PROPERLY (restrict browse=TRUE)
	search=function(taxon){
		url=paste("http://species.wikimedia.org/w/index.php?title=Special%3ASearch&search=",taxon,sep="")
		xx=try(system(paste("curl ", url, sep=""), intern=TRUE, wait=TRUE, ignore.stderr=TRUE), silent=TRUE)
		#pop(url)
		return(list(url=url, res=xx))
	}
	
	browse=function(taxon){
		url=paste("http://species.wikimedia.org/wiki",taxon,sep="/")
		xx=try(system(paste("curl ", url, sep=""), intern=TRUE, wait=TRUE, ignore.stderr=TRUE), silent=TRUE)
		#pop(url)
		return(list(url=url, res=xx))
	}
	
	disambig=function(taxon, html){
		contains_link=function(taxon,html){
			link=paste("href=\"/wiki/",taxon,sep="")
			zz=grep(link,html)
			yy=grep("aceae",html)
			return(intersect(zz,yy))
		}
		
		zz=grep(taxon, html)
		yy=grep("may mean:", html)
		ii=intersect(zz,yy)
		if(length(ii)) {
			l=contains_link(taxon,html)
			if(length(l)) return(split_wiki(html[l])) else return(NA)
		} else {
			return(NA)
		}
	}
	
	## SEARCH or BROWSE
	if(query_by_browse) query=browse else query=search
	tmp=suppressWarnings(query(taxon))
	xx=tmp$res
	url=tmp$url
	
	## parse output SPECIFIC to type of QUERY: 'family' in this case
	if(!inherits(xx,"try-error")){
		
		dd=disambig(taxon, xx)
		if(!is.na(dd)) return(dd)
		if(length(f<-grep("Familia",xx,ignore.case=FALSE))){ # BEST
			use=xx[f]
		} else if(length(f<-grep("aceae",xx,ignore.case=TRUE))){ # OKAY
			use=xx[f]
		} else {
			if(open) pop(url)
			return(NA)
		}
	}
	if(open) pop(url)
	return(unname(sapply(use, function(z) suppressWarnings(split_wiki(z)))))
}

### PLANTLIST.ORG ###
fetch_plantlist_taxonomy=function(){
	
	plantlist_families=function(){
		groups=list(A="angiosperms", G="gymnosperms", B="bryophytes", P="pteridophytes")
		
		extractor=function(string){
			gg=gsub("class=\"family\">", "SPLIT", gsub("</i></a>", "SPLIT", string, fixed=TRUE), fixed=TRUE)
			ss=unlist(strsplit(gg, "SPLIT"))
			return(ss[2])
		}
		
		collect_labs=function(html, header){
			a=grep(paste("/browse/",header,sep=""), html)
			b=grep("family", html)
			dat=html[intersect(a,b)]
			labs=unname(sapply(dat, extractor))
			labs
		}
		
		families=lapply(1:length(groups), function(idx){
						header=names(groups)[idx]
						url=paste("http://www.theplantlist.org/browse/",header,"/",sep="")
						html=readLines(url)
						collect_labs(html, header)
						})
		names(families)=unlist(groups)
		families
	} 
	
	plantlist_genera=function(){
		url="http://www.theplantlist.org/browse/-/-/"
		html=readLines(url)
		
		extractor=function(string){
			gg=gsub("class=\"genus\">", "SPLIT", gsub("</i></a> (<i class=\"family\">", "SPLIT", gsub("</i>)</li>", "SPLIT", string, fixed=TRUE), fixed=TRUE), fixed=TRUE)
			ss=unlist(strsplit(gg, "SPLIT"))
			return(ss[2:3])
		}
		
		collect_labs=function(html){
			a=grep("/browse/", html)
			b=grep("family", html)
			dat=html[intersect(a,b)]
			mm=matrix(NA, nrow=length(dat), ncol=2)
			for(i in 1:nrow(mm)){
				mm[i,]=extractor(dat[i])
			}
			colnames(mm)=c("genus", "family")
			return(mm)
		}
		
		genera=collect_labs(html)
		genera
	} 
	
	fam=plantlist_families()
	uu=unlist(sapply(1:length(fam), function(idx) {y=fam[[idx]]; names(y)=rep(names(fam)[idx], length(y)); y}))
	
	gen=plantlist_genera()		
	grp=unname(sapply(gen[,"family"], function(x) names(uu[uu==x])))

	dat=cbind(gen, grp)
	colnames(dat)=c("genus", "family", "group")
	dat
		
}

plantlist_taxonomy=function(group=c("angiosperm", "gymnosperm")){
	
	genus_url=function(genus="", tax, short){
		tmp=match(genus, tax$genus)
		return(paste(short, paste(rev(tax[tmp,]), collapse="/"), sep="/"))
	}
	
	browse_plantlist=function(genus, group, taxlist){
		# get data.frame for genus [using "Browse" in theplantlist.org]
		short=toupper(substring(group, 1, 1))
		base="http://www.theplantlist.org/browse/"
		url=paste(base, genus_url(genus, taxlist, short), sep="/")
		tmp=try(readHTMLTable(url, stringsAsFactors=FALSE), silent=TRUE)
		if(!inherits(tmp, "try-error")){
			if(class(tmp[[1]])=="data.frame"){
				dat=split_taxa(tmp[[1]])
				return(dat)
			} else {
				cat(paste("UNABLE to locate data for: ", genus, "\n", sep=""))
				return(NULL)
			}
		} else {
			return(NULL)
		}
	}
	
## PROCESSING ##
	plural=paste(group, "s", sep="")
	tax=read.table(paste("theplantlist", paste(group, "_genera.txt", sep=""), sep="/"), header=TRUE, stringsAsFactors=FALSE)
	
	families=split(tax, tax$family)
	family_info=lapply(1:length(families), function(fam_idx){
					   cat(paste(toupper(names(families)[fam_idx]), "\n", sep=""))
					   fam_df=families[[fam_idx]]
					   genera=fam_df$genus
					   gg=lapply(genera, function(genus) {
								 cat(paste("fetching data for: ", genus, "\n", sep=""))
								 browse_plantlist(genus, group, fam_df)
								 })
					   names(gg)=genera
					   return(gg)
					   })
	names(family_info)=names(families)
	return(family_info)
}

### PLANTLIST.ORG ###
fetch_plantlist_richnesses=function(genera, richnesses){
	
### FIXME: deal better with monotypic assumption if no match
#	richnesses: "spermatophyta_richnesses.csv" 
#	collects richnesses by genus (as given in richnesses or found online)
	query_plantlist_genus=function(genus){
		url=paste("http://www.theplantlist.org/tpl/search?q=",genus,sep="")
		# only use if genus is not in 'richnesses' (which is complete list of generic richnesses for accepted genera)
		# not an 'Accepted' genus (but possibly in database)
		tt=readHTMLTable(url)
		if(length(tt)){
			dd=tt[[1]]
			return(sum(!dd$Status%in%c("Synonym")))
		} else {
		# assume genus is monotypic
			return(1)
		}
		return(tt)
	}
	rr=sapply(genera, function(g){
			  if(g%in%richnesses$genus){
			  return(richnesses$species[richnesses$genus==g])
			  } else {
			  return(query_plantlist_genus(g))
			  }
			  })
	names(rr)=genera
	rr
}

### PLANTLIST.ORG ###
plantlist_richnesses=function(group=c("angiosperm", "gymnosperm")){
	
	## PROCESSING ##
	rda=paste(group, "_plantlist.rda", sep="")
	dat=get(load(rda))
	tax=read.table(paste("theplantlist", paste(group, "_genera.txt", sep=""), sep="/"), header=TRUE, stringsAsFactors=FALSE)
	tax$species=NA
	poor=data.frame(genus=NA, family=NA)
	for(i in 1:nrow(tax)){
		fam=tax$family[i]
		gen=tax$genus[i]
		dd=dat[[fam]][[gen]]
		tt=nrow(dd)
		if(is.null(tt)) {
			tt=NA
			poor=rbind(poor, c(gen, fam))
			print(c(gen, fam))
		}
		tax$species[i]=tt 
	}
	poor=poor[-1,]
	if(nrow(poor)>0) return(list(tax=tax, poor=poor)) else return(tax)
}

### PLANTLIST.ORG ###
plantlist_compilation=function(group=c("angiosperm", "gymnosperm")){
	
	## PROCESSING ##
	rda=paste(group, "_plantlist.rda", sep="")
	dat=get(load(rda))
	spp=data.frame(Species=NA, Authority=NA, Status=NA, Source=NA, Genus=NA, Family=NA)
	
	for(f in 1:length(dat)){
		if(class(dat[[f]])=="list"){
			for(g in 1:length(dat[[f]])){
				fam=names(dat)[f]
				gen=names(dat[[f]])[g]
				subdat=dat[[f]][[g]]
				subdat$Genus=gen
				subdat$Family=fam
				spp=rbind(spp, subdat)
			}
		} else {
			stop(paste("Unrecognized data format: ", names(dat)[f], sep=""))
		}
	}
	spp=spp[-1,]
	return(spp[,c("Species", "Genus", "Family", "Status", "Authority", "Source")])
}

### PLANTLIST.ORG ###
is.accepted=function(htmlstring){
	string="This name is the <a href=\"/about/#accepted\">accepted</a> name of a species in the genus"
	if(any(grep(string, htmlstring, fixed=TRUE))) return(TRUE) else return(FALSE)
}

### PLANTLIST.ORG ###
fetch_plantlist_synonyms=function(species){
	
	species=gsub("_","+",species)
	url=paste("http://www.theplantlist.org/tpl/search?q=",species,sep="")
	tt=htmlParse(url)
	if(length(tt)){
		syns=readHTMLTable(tt, stringsAsFactors=FALSE)
		if(length(syns)){
			syns=syns[[1]]
			syns=syns[!syns$Status=="Accepted",]
			rownames(syns)=NULL
			if(nrow(syns)>0) return(split_taxa(syns)) else return(NULL)
		} else {
			return(NULL)
		}
	} else {
		stop(paste(species, " is found in theplantlist.org.", sep=""))
	}
}

### PLANTLIST.ORG ###
query_plantlist=function(label){
#	label: should be species or genus
	taxon=gsub("_","+",label)
	url=paste("http://www.theplantlist.org/tpl/search?q=",taxon,sep="")
	htmlParse(url)
}

### PLANTLIST.ORG ###
plantlist.dataframe=function(plantlist){
## plantlist: a list object of records from plantlist
	
#	$Adoxaceae
#	$Adoxaceae$Adoxa
#	Species Authority   Status           Source
#	1 Adoxa_moschatellina        L. Accepted WCSP (in review)
#	2   Adoxa_xizangensis     G.Yao Accepted WCSP (in review)
#	
#	$Adoxaceae$Sambucus
#	Species               Authority   Status           Source
#	1        Sambucus_adnata            Wall. ex DC. Accepted              TRO
#	2  Sambucus_australasica        (Lindl.) Fritsch Accepted              TRO
#	3     Sambucus_australis        Cham. & Schltdl. Accepted              TRO
#	4    Sambucus_canadensis                      L. Accepted WCSP (in review)
#	5     Sambucus_chinensis                  Lindl. Accepted              TRO
#	6         Sambucus_nigra                      L. Accepted              TRO
	
	xx=sum(sapply(plantlist, function(x) sum(sapply(x, function(y) nrow(y)))))
	mm=matrix(NA, nrow=xx, ncol=3)
	ct=0
	for(i in 1:length(plantlist)){
		fam=names(plantlist)[i]
		for(j in 1:length(plantlist[[i]])){
			gen=names(plantlist[[i]])[j]
			spp=plantlist[[i]][[j]][,1]
			mm[(ct+1):(ct+length(spp)),]=cbind(spp, gen, fam)
			ct=ct+length(spp)
		}
	}
	mm=data.frame(mm, stringsAsFactors=FALSE)
	names(mm)=c("species","genus","family")
	mm
}
