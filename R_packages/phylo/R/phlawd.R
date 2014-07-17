


keepfile_phlawd.phylo=function(phy, locus=""){
	
	## GRABS sequencing for species spanning the root of 'phy'
	## proceeds iteratively across tree until data are found
	
	.exemplar_seq=function(seq){
		
		# use width of pairwise hits to determine the best exemplar_seq
		if(length(seq)==1) return(seq)
		dat=blastn(seq,seq)
		dat=dat[dat$query!=dat$subject,]
		if(nrow(dat)){
			widths=c(by(dat$width, dat$query, mean))
			choice=names(widths)[min(which(widths==max(widths)))]
			first=max(dat$qstart[dat$query==choice])
			last=min(dat$qend[dat$query==choice])
			if(first<last) return(subseq(seq[which(names(seq)==choice)], first, last)) else return(NULL)
		} else {
			return(NULL)
		}
	}
	
	.trim_seq=function(seq){
		
		# use width of pairwise hits to determine the best exemplar_seq
		if(length(seq)==1) return(seq)
		dat=blastn(seq,seq)
		dat=dat[dat$query!=dat$subject,]
		if(nrow(dat)){
			seq=seq[names(seq)%in%dat$query]
			nm=names(seq)
			idx=lapply(nm, function(x) {
				   cur=dat[dat$query==x,]
				   first=max(cur$qstart)
				   last=min(cur$qend)
				   return(c(first,last))
				   })
			for(i in 1:length(seq)){
				ss=c(idx[[i]][1], idx[[i]][2])
				if(ss[1]>=ss[2]) return(NULL)
				seq[i]=subseq(seq[i], ss[1], ss[2]) 
			}
			return(seq)
		} else {
			return(NULL)
		}
	}
	
	N=Ntip(phy)+1
	
	if(locus=="") stop("Please supply a 'locus' for which to fetch data")
	
	seqs=sapply(get.desc.of.node(N, phy), function(x) {
		if(x<=Ntip(phy)) y=x else y=get.descendants.of.node(x, phy, TRUE)
		i=1
		while(i <= length(y)){
		   cur=phy$tip.label[y[i]]
		   res=fetch_ncbi(cur, locus=locus, outformat=21345, align=FALSE)
		   if(!is.null(res)){
				res=res[grep(tolower(locus), tolower(names(res)))]
				if(z<-length(res)){
					if(z==1) break()
					contained_taxa=try(fetch_gbtaxonomy.below(cur, "species"),silent=TRUE)
					if(!inherits(contained_taxa, "try-error")) {
						tt=strsplit(names(res), "_")
						check=sapply(tt, function(x) paste(x[3:4], collapse=" ")%in%contained_taxa)
						res=res[check]
						if(length(res)) {
							res=.exemplar_seq(res)
							if(!is.null(res)) break()
						} else {
							res=NULL
						} 
					} else {
						res=NULL
					}
				} else {
					res=NULL
				}
		   }
		   i=i+1
		}
		if(is.null(res)) return(NA) else return(res)	
	})
	
	if(!any(is.na(seqs))){
		ss=sapply(seqs, toString)
		names(ss)=sapply(seqs, names)
		res=DNAStringSet(ss)
		return(.trim_seq(res))
	} else {
		return(NULL)
	}
}


extract.partitions=function(nexusfile.phlawd){
	ss=scan(nexusfile.phlawd, nlines=3, what="char")
	ss=ss[!ss%in%c("#NEXUS","BEGIN","DATA;","]")]
	tt=sapply(ss, function(x) {y=gsub("[", "", x, fixed=TRUE); y=unlist(strsplit(gsub("."," ",y,fixed=TRUE), " ", fixed=TRUE))[1]})
	if(!length(tt)%%2==0) stop("Unknown error.")
	ll=seq(1,length(tt),by=2)
	rr=lapply(ll, function(x) as.numeric(unlist(strsplit(tt[x+1], "-", fixed=TRUE))))
	lims=t(data.frame(rr))
	colnames(lims)=c("L","R")
	lims=cbind(gene=unname(tt[ll]), lims)
	lims=data.frame(lims, stringsAsFactors=FALSE)
	rownames(lims)=NULL
	class(lims[,2])<-class(lims[,3])<-"numeric"
	return(lims)	
}

write.partition.file=function(partitions.df, file=""){
	df=partitions.df
	if(file.exists(file)) unlink(file)
	txt=paste("DNA, ", paste(df[,1], "=", paste(df[,2],df[,3],sep="-"), sep=" "), sep="")
	for(i in 1:length(txt)){
		write(txt[i], file=file, append=TRUE)
	}
}

