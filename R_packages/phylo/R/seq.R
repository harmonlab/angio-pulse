swap.ambiguities=function(dnastring){
	## expecting ambiguities to be enclosed in curly braces ("{","}")
	xx=unlist(strsplit(dnastring,""))
	open<-which(xx=="{")
	close<-which(xx=="}")
	dnastring=gsub("?","N",dnastring,fixed=TRUE)
	if((length(close)==length(open)) && (length(close)>0 & length(open)>0)){
		breaks=rbind(open,close)
		ambiguity=function(string) {
			res=try(mergeIUPACLetters(string),silent=TRUE)
			if(inherits(res,"try-error")){
				warning("some strings in 'DNA' contain non-IUPAC letters; substituting 'N'")
				res="N"
			}
			res
		}
		tot=nchar(dnastring)
		new=sapply(1:ncol(breaks), function(idx) {
			init=breaks[1,idx]
			tail=breaks[2,idx]
			string=subseq(dnastring, init+1, tail-1)  
			amb=ambiguity(string)
			if(idx==1){
				preceding=ifelse(init==1, "", subseq(dnastring, 1,init-1))
			} else {
				preceding=subseq(dnastring, breaks[2,idx-1]+1,init-1)
			}
			res=paste(c(preceding,amb),collapse="",sep="")

			if(idx==ncol(breaks)){
				 ending=ifelse(tail==tot, "", subseq(dnastring, tail+1, tot)) 
				 res=paste(c(res,ending),collapse="",sep="")
			}
			res
		})
		
		dnastring=paste(new, collapse="")
		
	} else if((length(close)>0 | length(open)>0)){
		stop("encountered DNA string with apparent but insoluble ambiguities")
	} 
	
	return(dnastring)
}

read.seq=function(file,format=c("fasta","phylip","nexus")){
#	if(length(format)>1) stop("format of the datafile must be specified")
	datafile=file
	format=match.arg(format, c("fasta","nexus","phylip"))
	if(format=="fasta") return(read.DNAStringSet(datafile))
	dat=readLines(datafile)
	if(format=="phylip") {
		check=as.numeric(unlist(strsplit(dat[1]," ")))
		if(length(check)!=2) stop("Encountered poorly formatted phylip file.")
		taxa=check[1]
		dat=dat[-1]
		blank=which(dat=="")
		if(length(blank)) dat=dat[-blank]
		if(length(dat)%%taxa!=0) {
			stop("Encountered poorly formatted phylip file.")
		} else {
			if(length(dat)>taxa) {
				dat=sapply(1:taxa, function(x) {
						   ii=seq(x,length(dat),by=taxa)
						   paste(dat[ii],collapse="")
						   }
						   )
				
			}
		}
		names=list()
		msa=list()
		for(i in 1:length(dat)) {
			if(dat[i]!=""){
				y=strsplit(dat[i],"\t",fixed=TRUE)[[1]]
				if(length(y)==1) y=strsplit(y," ")[[1]]
				names[[i]]=y[1]
				msa[[i]]=y[length(y)]
			} else {
				names[[i]]=NULL
				msa[[i]]=NULL
			}
		}
	} else if(grep("#NEXUS",dat)){
		flag=FALSE
		if(format!="nexus") warning("Detected 'nexus' formatted data.")
		check=grep("DIMENSIONS", dat, ignore.case=TRUE, value=TRUE)
		if(!length(check)) {
			flag=TRUE
		} else {
			tmp=unlist(strsplit(gsub(";", " ", gsub("=", " ", check)), " ", fixed=TRUE))
			check=c(as.numeric(tmp[grep("NTAX", tmp, ignore.case=TRUE)+1]), as.numeric(tmp[grep("NCHAR", tmp, ignore.case=TRUE)+1]))
			taxa=check[1]
			nchar=check[2]
			blank=dat%in%c("", ";")
			dat=dat[!blank]
			terminus=min(grep("End;", dat, ignore.case=TRUE))
			dat=dat[-c(1:grep("MATRIX", dat, ignore.case=TRUE),terminus:length(dat))]
			if(length(dat)!=taxa) flag=TRUE
		}
		
		if(flag) warning("Encountered poorly formatted nexus file.")
		
		names=list()
		msa=list()
		for(i in 1:length(dat)) {
			if(dat[i]!=""){
				y=strsplit(dat[i],"\t",fixed=TRUE)[[1]]
				y=y[y!=""]
				if(length(y)==1) y=strsplit(y," ",fixed=TRUE)[[1]]
				y=y[y!=""]
				names[[i]]=y[1]
				z=swap.ambiguities(y[length(y)])
				if(nchar(z)==nchar) msa[[i]]=z else stop("Encountered 'nchar' problem in nexus file.")
			} else {
				names[[i]]=NULL
				msa[[i]]=NULL
			}
		}		
	} else {
		stop("format not yet supported")
	}
	dat=NULL
	msa=DNAStringSet(unlist(msa))
	names(msa)=unlist(names)	
	return(msa)
}


write.seq=function(DNAStringSet,format=c("fasta","phylip","nexus"),file){
	format=match.arg(format, c("fasta","nexus","phylip"))
	
	dat=DNAStringSet
	if(format=="phylip") {
		ll=width(dat)
		if(length(unique(ll))==1) {
			ll=unique(ll)
			tt=length(dat)
			write(paste(tt,ll,sep=" "),file)
			for(i in 1:length(dat)) {
				write(paste(names(dat)[i],dat[i],collapse="\t"),file,append=TRUE)
			}
		} else {
			stop("DNAStringSet appears to be a ragged alignment.")
		}	
	} else if(format=="nexus") {
		write.seq.to.nexus(dat, file)
	} else if(format=="fasta"){
		write.XStringSet(dat, filepath=file)
	} else {
		stop("format not yet supported")
	}
}

write.seq.to.nexus=function (DNAStringSet, file) 
{
	x=DNAStringSet
    indent <- "  "
    maxtax <- 5
    defcharsperline <- 80
    defgap <- "-"
    defmissing <- "?"
    ntax <- length(x)
	ll=width(x)
	if(length(unique(ll))==1) {
		nchars=unique(ll)
	} else {
		stop("DNAStringSet appears to be a ragged alignment.")
	}	
    zz <- file(file, "w")
    if (is.null(names(x))) {
        names(x) <- as.character(1:ntax)
    }
    fcat <- function(..., file = zz) {
        cat(..., file = file, sep = "", append = TRUE)
    }
    fcat("#NEXUS\n")
    NCHAR <- paste("NCHAR=", nchars, sep = "")
    NTAX <- paste("NTAX=", ntax, sep = "")
	DATATYPE <- "DATATYPE=DNA"
    GAP <- paste("GAP=", defgap, sep = "")
	MISSING <- paste("MISSING=", defmissing, sep="")
	INTERLEAVE <- "INTERLEAVE=NO"
	fcat("BEGIN DATA;\n")
	fcat(indent, "DIMENSIONS", " ", NTAX, " ", NCHAR, ";\n")
	fcat(indent, "FORMAT", " ", DATATYPE, " ", MISSING, " ", GAP, " ", INTERLEAVE, ";\n")
	fcat(indent, "MATRIX\n")
	for(i in 1:length(x)){
		fcat(paste(names(x)[i], indent, toString(x[[i]]), "\n", sep=""))
	}
	fcat(indent, ";\n")
	fcat("END;\n\n")
   
    close(zz)
}

chimerize.seq=function(partitions, dat, possible.labels){
	# creates chimeric sequence from several 'possible.labels'
	# partitions: must have 'L' and 'R' columns for start and end of a locus
	subdat=dat[names(dat)%in%possible.labels]
	width=unique(width(dat))
	if(length(width)>1) stop("data appear to be of a ragged alignment.")
	subseqdata=lapply(1:nrow(partitions), function(gx){
					  lims=unlist(partitions[gx, c("L","R")])
					  nsubdat=narrow(subdat, lims[1], lims[2])
					  aa=seq.completeness(nsubdat)
					  toString(nsubdat[min(which(aa==max(aa)))])
					  })
	fetched=paste(unlist(subseqdata), collapse="")
	if(nchar(fetched)==width) return(fetched) else stop("Unknown error.")
}

dropdup.seq=function(seqs){
	f=duplicated(seqs,,fromLast=FALSE)
	r=duplicated(seqs,fromLast=TRUE)
	a=apply(alphabetFrequency(seqs), 1, function(x) if(sum(x[1:4])>0) return(TRUE) else return(FALSE))
	return(seqs[!f & !r & a])
}

reduce.seq=function(allgenes, seqs, genes){
	to.narrow=match(genes,allgenes$gene)
	if(all(allgenes$gene%in%genes)) return(seqs)	
	subseqs=sapply(1:length(seqs), function(x) {
				   subseq=seqs[x]
				   ss=c()
				   for(i in 1:length(to.narrow)){
				   ss=c(ss,toString(narrow(subseq, as.numeric(allgenes$L[to.narrow[i]]), as.numeric(allgenes$R[to.narrow[i]]))))
				   }
				   return(paste(ss,collapse=""))
				   })
	newseq=DNAStringSet(subseqs)
	names(newseq)=names(seqs)
	return(newseq)
}

write.partition.file=function(partitions.df, file=""){
	df=partitions.df
	if(file.exists(file)) unlink(file)
	txt=paste("DNA, ", paste(df[,1], "=", paste(df[,2],df[,3],sep="-"), sep=" "), sep="")
	for(i in 1:length(txt)){
		write(txt[i], file=file, append=TRUE)
	}
}


pool.fna=function(dir, coverage=3){
# coverage: trim sites to only those with representatives from at least 'coverage' taxa
	files=dir(dir)
	loci=sapply(strsplit(files, ".", fixed=TRUE), function(s) toupper(s[1]))
	uniq=unique(loci)
	if(coverage>1) fnas="*.trimmed.fna" else fnas="*.pooled.fna"
	
	for(locus in uniq){
		ff=paste(dir, files[which(loci==locus)], sep="/")
		poolfile=paste(dir, paste(locus, "pooled", sep="."), sep="/")
		scat(ff, outfile=poolfile)
		system(paste("mafft --auto", poolfile, " > ", paste(poolfile, "fna", sep="."), sep=" "), ignore.stderr=TRUE, wait=TRUE)
		if(coverage>1){
			tmp=read.DNAStringSet(paste(poolfile, "fna", sep="."))
			dat=narrow.seqs(tmp, coverage=coverage)
			trimfile=paste(poolfile, "trimmed", sep=".")
			write.XStringSet(dat, paste(trimfile, "fna", sep="."))
		}
	}
	system(paste("phyutility -concat -in", paste(dir, fnas, sep="/"), "-out POOLED.nexus", sep=" "))
}

seq.completeness=function(DNAStringSet){
	a=alphabetFrequency(DNAStringSet)
	rowSums(a[,1:4])/width(DNAStringSet)
}

rename.seqs=function(DNAStringSet){
	nn=names(DNAStringSet)
	nn=sapply(strsplit(nn, "_", fixed=TRUE), function(x) paste(x[2:3],collapse="_"))
	nn
}

trim.seqs=function(DNAStringSet, species, rename=TRUE){
# trim DNAStringSet to one exemplar for the given species based on sequence coverage
	dat=DNAStringSet
	if(rename) names(dat)=rename.seqs(dat)
	dat=dat[names(dat)%in%species]
	if(length(dat)){
		uniq=table(names(dat))
		uniq=names(uniq[uniq>1])
		if(length(uniq)){
			for(taxon in uniq){
				qual=seq.completeness(dat[names(dat)==taxon])
				keep=min(which(qual==max(qual)))
				ww=which(names(dat)==taxon)
				qq=(1:length(qual))[-keep]
				drop=ww[qq]
				dat=dat[-drop]
			}
		}
	}
	dat
}

narrow.seqs=function(DNAStringSet, coverage=2){
# remove sites with fewer than 'coverage' individuals with data
	dat=DNAStringSet
	s=width(dat)
	if(length(dat)<coverage) return(DNAStringSet)
	if(length(w<-unique(s))>1){
		if(min(w)==0) stop("Some reads have zero length.")
		warning(paste("Assuming homology for the original", min(w), "bases.", sep=" "))
		dat=narrow(dat, 1, ww)
	}
	ww=min(w)
	mat=as.matrix(dat)
	coverages=apply(mat, 2, function(x) sum(x!="-"))
	mat=mat[,coverages>=coverage]
	strings=apply(mat, 1, paste, collapse="")
	names(strings)=rownames(mat)
	new=DNAStringSet(strings)
	return(new)
}

nexus_to_DNAStringSet=function(nex){
# nex is a list read in by ape:::read.nexus.data()
	x=sapply(nex, paste, collapse="")
	dat=DNAStringSet(x)
	dat
}

nexus_to_fasta=function(nex, file=""){
	dat=nexus_to_DNAStringSet(nex)
	write.XStringSet(dat, file)
}


