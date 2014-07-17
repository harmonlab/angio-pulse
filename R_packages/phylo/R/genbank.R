collect_gbtaxonomy=function(phy, rank=list(from="family", to="phylum")){
	require(phylo)
	
	# phy: tip.labels should be of same rank
	if(check.multicore()) fun=mclapply else fun=lapply
	labels=unique(short<-sapply(orig<-phy$tip.label, function(x) unlist(strsplit(x, "_", fixed=TRUE))[1]))
	tax=fun(phy$tip.label, function(x) {
			   z=try(fetch_gbtaxonomy.above(x, rank=rank$to, update.dmp=FALSE),silent=TRUE)
			   if(!is.null(z) & !inherits(z, "try-error")){
					y=z$short
					ff=which(rownames(y)==rank$from)
					tt=which(rownames(y)==rank$to)
					if(all(sapply(c(ff, tt), length))){
						w=y[ff:tt,]
						names(w)=rownames(y)[ff:tt]
						return(w)
					} else {
						return(NULL)
					}
			   } else {
					return(NULL)
			   }
	})
	xx=sapply(tax, is.null)
	if(any(xx) & rank$from!="species"){
		labels=unique(short<-sapply(orig<-phy$tip.label[xx], function(x) unlist(strsplit(x, "_", fixed=TRUE))[1]))
		for(i in 1:length(labels)){
			z=try(fetch_gbtaxonomy.above(labels[i], rank=rank$to, update.dmp=FALSE),silent=TRUE)
			if(!is.null(z) & !inherits(z, "try-error")){
				y=z$short
				ff=which(rownames(y)==rank$from)
				tt=which(rownames(y)==rank$to)
				if(all(sapply(c(ff, tt), length))){
					w=y[ff:tt,]
					names(w)=rownames(y)[ff:tt]
					tax[[which(xx)[i]]]=w
				} 
			} 
		}
	}	
	
	ll=sapply(tax, length)
	nc=unique(ll[ll>0])
	if(!length(nc)==1) 	stop("gb_taxonomies incongruent")	
	nn=names(tax[[min(which(ll==nc))]])

	res=matrix(NA, nrow=length(tax), ncol=nc)
	for(i in 1:length(tax)) {
		if(!is.null(tax[[i]])){
			res[i,]=tax[[i]]
		}
		
	}
	colnames(res)=nn
	drop=apply(res, 2, function(x) all(is.na(x)))
	rownames(res)=phy$tip.label
	
	return(res)
}

fetch_gbtaxonomy.below=function(taxon, rank="genus", update.dmp=FALSE){
## returns taxonomic information for taxa contained within 'taxon', from the 'rank' given

	gb_path=phylo:::.path_to_gb.backward()
	if(update.dmp) phylo:::.fetch_taxdump.backward()
	cur=system(paste(gb_path[["pl"]], "-idx", gb_path[["idx_dir"]], "-download", gb_path[["download_dir"]], "-taxon", taxon, "-rank", rank, "-outfile .tmp.txt", sep=" "), ignore.stderr=TRUE, intern=TRUE)
	if(inherits(cur, "try-error")) {
		taxa=NA
	} else {
		res=read.delim(".tmp.txt", header=FALSE, stringsAsFactors=FALSE)
		if(nrow(res)){
			taxa=res[,2]
			names(taxa)=res[,1]
		}
		unlink(".tmp.txt")
	}
	taxa
}

fetch_gbtaxonomy.above=function(taxon, rank="genus", update.dmp=FALSE){
	
## returns taxonomic information for a 'taxon' up to the 'rank' given
## requires fetch_genbank.pl and potentially nodes.dmp and names.dmp (if /tmp/idx is not available)
	
	Linnaean=c(
#		       "root",
			   "superkingdom",
			   "kingdom",
			   "superphylum",
			   "phylum",
			   "subphylum",
			   "superclass",
			   "class",
			   "subclass",
			   "infraclass",
			   "superorder",
			   "order",
			   "suborder",
			   "parvorder",
			   "infraorder",
			   "superfamily",
			   "family",
			   "subfamily",
			   "tribe",
			   "subtribe",
			   "genus",
			   "subgenus",
#			   "species group",
#		       "species subgroup",
			   "species",
			   "subspecies",
			   "varietas",
			   "forma"
			   )
	
	if(!is.null(rank)) if(!rank%in%Linnaean) {
		cat(paste("Ensure that supplied 'rank' is in: \n\t", paste(Linnaean, collapse="\n\t"), "\n", sep=""))
		stop("Supplied 'rank' is unrecognized.")
	}
	
	# workhorse function
	grab_ncbi=function(taxon, id=TRUE, gb_path=NULL){
		if(!id) taxon=gsub(" ","_",taxon)
		if(is.null(gb_path)){
			gb_path=phylo:::.path_to_gb.dmp()
			if(!all(sapply(gb_path, file.exists))){
				phylo:::.fetch_taxdump.forward(mv=TRUE)
			}			
		}
		cur=system(paste(gb_path[["pl"]], "--nodes", gb_path[["nodes"]], "--names", gb_path[["names"]], taxon, sep=" "), ignore.stderr=TRUE, intern=TRUE)
		if(length(cur)){
			names(cur)=c("id", "rank", "taxon", "anc")
			return(cur)
		} else {
			return(NULL)
		}
	}
	
	store_ncbi=function(cur_dat){
		res=cur_dat[["taxon"]]
		names(res)=cur_dat[["rank"]]
		res
	}
	
	
	# PROCESSING
	gb_path=phylo:::.path_to_gb.dmp()
	if(!all(sapply(gb_path[c("nodes","names")], file.exists)) | update.dmp==TRUE){
		phylo:::.fetch_taxdump.forward(mv=TRUE)
	}	
	
	cur_dat=grab_ncbi(taxon, id=FALSE)
	if(is.null(cur_dat)) return(cur_dat)
	stored=store_ncbi(cur_dat)
	
	if((init_rank<-match(cur_dat[["rank"]], Linnaean))>(last_rank<-match(rank, Linnaean))){
		while(1){
			new_dat=grab_ncbi(cur_dat[["anc"]], id=TRUE)
			stored=c(stored, store_ncbi(new_dat))
			if(new_dat[["rank"]]==rank) {
				break()
			} else {
				nn=(match(new_dat[["rank"]], Linnaean))<match(rank, Linnaean)
				if(!is.na(nn) & nn==TRUE){
					break()
				}
			}
			cur_dat=new_dat
		}
	} else {
		stop("Supplied 'rank' is inconsistent with supplied 'taxon'.")
	}
	
	tmp=Linnaean[init_rank:last_rank]
	short=rep(NA, length(tmp))
	names(short)=tmp
	m=match(names(stored), names(short))
	short[m[!is.na(m)]]=stored[names(stored)%in%names(short)]
	short=data.frame(short, stringsAsFactors=FALSE)
	names(short)="rank"
	return(list(short=short, all=stored))
}

.path_to_gb.dmp=function(){
	paths=resolve.executable()
	pl=paths[[which(names(paths)=="fetch_genbank")]]
	base_path=gsub("fetch_genbank.pl", "", pl)
	nd=paste(base_path, "nodes.dmp", sep="")
	nm=paste(base_path, "names.dmp", sep="")
	return(c(names=nm, nodes=nd, pl=pl))
}

.path_to_gb.backward=function(){
	paths=resolve.executable()
	pl=paths[[which(names(paths)=="containedtaxa_genbank")]]
	base_path=gsub("containedtaxa_genbank.pl", "", pl)
	return(c(pl=pl, idx_dir=paste(base_path,"taxonomy",sep="/"), download_dir=base_path))
}


fetch_ncbi=function(taxa, locus="", max=c(length=4000, seqs=300), outformat=231, align=FALSE){
	if(!system("which phyutility", ignore.stdout=TRUE)==0) stop("Install 'phyutility' before proceeding.")
	arg=paste(taxa,collapse="+OR+")

	if(locus!="") {
		arg=paste(arg, locus, sep="+AND+")
	} else {
		align=FALSE
	}

	res=".tmpPHYUTILITY"
	if(file.exists(res)) unlink(res)
	unlink(".phyutility.log")
	
	check=as.numeric(unlist(strsplit(system(paste("phyutility -es -db 1 -log .phyutility.log -out .phyutilitySearch -term", arg, sep=" "), intern=TRUE), " = "))[2])
	if(check==0) return(NULL)
	if(check>max[["seqs"]]) stop(paste("Encountered", check, "records; consider increasing 'max' seqs"))
	
	system(paste("phyutility -ef -db 1 -log .phyutility.log -outfor", outformat, "-out", res, "-ll", max[["length"]], "-term", arg, sep=" "))
	unlink(".phyutility.log")

	if(is.empty(res)) {
		unlink(res)
		return(NULL)
	} else {
		if(align){
			mafftfile=".tmpMAFFT"
			system(paste("mafft --auto", res, " > ", mafftfile, sep=" "), ignore.stderr=TRUE, wait=TRUE)
			if(!is.empty(mafftfile)){
				dat=read.DNAStringSet(mafftfile)
				unlink(mafftfile)
			} 
		} else {
			dat=read.seq(res, format="fasta")
			
		}
		unlink(res)
		return(dat)
	}
}

is.DNAStringSet=function(obj){
	any(class(obj)%in%"DNAStringSet")
}

blastn=function(ref.fa="", qry.fa="", outformat=6, ...){
# 'ref.fa' and 'qry.fa': can be file or DNAStringSet
	require(Biostrings)
	
	rm.ref<-rm.qry<-FALSE
	
	if(is.DNAStringSet(ref.fa)) {write.XStringSet(ref.fa, file=".ref.fa"); ref.fa=".ref.fa"; rm.ref=TRUE}
	if(is.DNAStringSet(qry.fa)) {write.XStringSet(qry.fa, file=".qry.fa"); qry.fa=".qry.fa"; rm.qry=TRUE}
	
	if(system("which blastn", ignore.stdout=TRUE)==0){
		tmp=system(paste("blastn",paste("-subject",ref.fa,sep=" "),paste("-query",qry.fa,sep=" "),"-outfmt", outformat," > .tmp",sep=" "))
		if(outformat%in%c(6,7,10)) {
			if(outformat==10) tmp=read.csv(".tmp", as.is=TRUE) else tmp=read.table(".tmp", as.is=TRUE)
			names(tmp)=c("query","subject","pid","width","mismatch","gapopen","qstart","qend","sstart","send","e","bitscore")
		} else {
			tmp=readLines(".tmp")
		}
		unlink(".tmp")
	} else {
		warning("Ensure that ~/.Rprofile and thus 'Sys.getenv()' involves the appropriate PATH.")
		stop("blastn does not appear visible with the system $PATH")
	}
	
	if(rm.ref) unlink(ref.fa)
	if(rm.qry) unlink(qry.fa)
	
	return(tmp)
}


.fetch_taxdump.forward=function(mv=FALSE){
	cat("Please be patient as 'taxdump' is built from NCBI; this should take several minutes...\n")
	if(!system("which wget", ignore.stdout=TRUE)==0) stop("Install 'wget' before proceeding.")
	if(!system("which gunzip", ignore.stdout=TRUE)==0) stop("Install 'gunzip' before proceeding.")
	if(!system("which tar", ignore.stdout=TRUE)==0) stop("Install 'tar' before proceeding.")

	system("wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz", ignore.stderr=TRUE, wait=TRUE)
	system("gunzip taxdump.tar.gz", ignore.stderr=TRUE, wait=TRUE)
	system("tar -xvf taxdump.tar", ignore.stderr=TRUE, wait=TRUE)
	
	# move 'nodes.dmp' and 'names.dmp' to appropriate locations
	if(mv){
		system("rm -rf /tmp/idx")
		gb_path=phylo:::.path_to_gb.dmp()
		if(all(sapply(ff<-c("nodes.dmp", "names.dmp"), file.exists))){
			system(paste("mv nodes.dmp", gb_path[["nodes"]], sep=" "))
			system(paste("mv names.dmp", gb_path[["names"]], sep=" "))
			cleanup=c("taxdump.tar", "readme.txt", "gc.prt", "merged.dmp", "division.dmp", "delnodes.dmp", "citations.dmp", "gencode.dmp")
			cleanup=cleanup[cleanup%in%dir()]
			system(paste("rm -f", paste(cleanup, collapse=" "), sep=" "))
		} else {
			stop("Cannot locate 'nodes.dmp' and (or) 'names.dmp'.")
		}
	}
}

.fetch_taxdump.backward=function(){
	cat("Please be patient as 'taxdump' is built from NCBI; this should take several minutes...\n")
	if(!system("which ftp", ignore.stdout=TRUE)==0) stop("Install 'ftp' before proceeding.")
#	if(!system("which gunzip", ignore.stdout=TRUE)==0) stop("Install 'gunzip' before proceeding.")
#	if(!system("which tar", ignore.stdout=TRUE)==0) stop("Install 'tar' before proceeding.")
	
#	system("wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz", ignore.stderr=TRUE, wait=TRUE)
#	system("gunzip taxdump.tar.gz", ignore.stderr=TRUE, wait=TRUE)
#	system("tar -xvf taxdump.tar", ignore.stderr=TRUE, wait=TRUE)
	
	# move 'nodes.dmp' and 'names.dmp' to appropriate locations
	gb_path=phylo:::.path_to_gb.backward()
	cur=system(paste(gb_path[["pl"]], "-idx", gb_path[["idx_dir"]], "-download", gb_path[["download_dir"]], "-updateonly", sep=" "), ignore.stderr=TRUE, intern=TRUE)

#	if(mv){
#		system("rm -rf /tmp/idx")
#		gb_path=phylo:::.path_to_gb.dmp()
#		if(all(sapply(ff<-c("nodes.dmp", "names.dmp"), file.exists))){
#			system(paste("mv nodes.dmp", gb_path[["nodes"]], sep=" "))
#			system(paste("mv names.dmp", gb_path[["names"]], sep=" "))
#			cleanup=c("taxdump.tar", "readme.txt", "gc.prt", "merged.dmp", "division.dmp", "delnodes.dmp", "citations.dmp", "gencode.dmp")
#			cleanup=cleanup[cleanup%in%dir()]
#			system(paste("rm -f", paste(cleanup, collapse=" "), sep=" "))
#		} else {
#			stop("Cannot locate 'nodes.dmp' and (or) 'names.dmp'.")
#		}
#	}
}

