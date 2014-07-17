

write.r8s=function(phy=NULL, calibrations, base=""){
#	calibrations: dataframe with minimally 'MRCA' 'MaxAge' 'MinAge' 'taxonA' and 'taxonB' from .build_calibrations
#	MRCA							MaxAge     MinAge                                  taxonA                                  taxonB valid
#	c65bacdf65aa29635bec90f3f0447c6e 352.234677 352.234677                          Inga_chartacea             Encephalartos_umbeluziensis  TRUE
#	d4bc6557dccbd4e8e18b867979f34f8e 243.269677 243.269677                          Inga_chartacea                     Nuphar_sagittifolia  TRUE
#	5e30c16531694aff7d94da3b35806677 217.632627 217.632627                          Inga_chartacea                  Schisandra_glaucescens  TRUE
	
	if(file.exists(inp<-paste(base,"infile",sep="."))) unlink(inp)
	
##	check appropriateness of constraints ##
#	check for 'calibrations' and 'phy' congruence
	if(!is.null(phy)){
		check=function(t, phy) all(t%in%phy$tip.label)
		a=check(calibrations$taxonA, phy)
		b=check(calibrations$taxonB, phy)
		
		if(!all(c(a,b))) stop("Some calibrations not encountered in tree.")		
	}
	
##	build PATHd8 file
#	calibrations$fixage=ifelse(calibrations$MinAge==calibrations$MaxAge, TRUE, FALSE)
	constraints<-constraintnames<-character(nrow(calibrations))
	for(i in 1:nrow(calibrations)){
		cal=calibrations[i,]
		taxon=cal$MRCA
		desc=c(cal$taxonA, cal$taxonB)
		
		txt1=paste("min =", taxon, cal$MinAge, sep=" ")
		txt2=paste("max =", taxon, cal$MaxAge, sep=" ")
		txt=paste(txt1,txt2,sep="\n")
		
		constraints[i]=txt
		constraintnames[i]=paste("mrca", taxon, desc[1], desc[2], sep=" ")
	}
#	phy$node.label=NULL
	infile=list(
				# tree=write.tree(phy),
				names=paste(unlist(constraintnames), collapse="\n"),
				mrca=paste(unlist(constraints), collapse="\n")
				)
	
	inp=paste(base,"infile",sep=".")
	writeLines(paste(infile,collapse="\n\n"), con=inp)	
	return(inp)
}


write.pathd8=function(phy, calibrations, base=""){
#	calibrations: dataframe with minimally 'MRCA' 'MaxAge' 'MinAge' 'taxonA' and 'taxonB' 
	
	if(file.exists(inp<-paste(base,"infile",sep="."))) unlink(inp)
	
##	check appropriateness of constraints ##
#	check for 'calibrations' and 'phy' congruence
	check=function(t, phy) all(t%in%phy$tip.label)
	a=check(calibrations$taxonA, phy)
	b=check(calibrations$taxonB, phy)
	
	if(!all(c(a,b))) stop("Some calibrations not encountered in tree.")
	
##	build PATHd8 file
	calibrations$fixage=ifelse(calibrations$MinAge==calibrations$MaxAge, TRUE, FALSE)
	constraints<-constraintnames<-character(nrow(calibrations))
	for(i in 1:nrow(calibrations)){
		cal=calibrations[i,]
		taxon=cal$MRCA
		desc=c(cal$taxonA, cal$taxonB)
		if(cal$fixage) {
			txt=paste("mrca:", paste(desc[1], ", ", desc[2], ", ", sep=""), paste("fixage=", cal$MinAge, ";", sep=""), sep=" ")
		} else {
			txt1=paste("mrca:", paste(desc[1], ", ", desc[2], ", ", sep=""), paste("minage=", cal$MinAge, ";", sep=""), sep=" ")
			txt2=paste("mrca:", paste(desc[1], ", ", desc[2], ", ", sep=""), paste("maxage=", cal$MaxAge, ";", sep=""), sep=" ")
			txt=paste(txt1,txt2,sep="\n")
		}
		constraints[i]=txt
		constraintnames[i]=paste("name of mrca: ", paste(desc[1], ", ", desc[2], ", ", sep=""), paste("name=", cal$MRCA, ";", sep=""), sep=" ")
	}
	phy$node.label=NULL
	infile=list(tree=write.tree(phy),
				mrca=paste(unlist(constraints), collapse="\n"),
				names=paste(unlist(constraintnames), collapse="\n")
				)
	
	inp=paste(base,"infile",sep=".")
	writeLines(paste(infile,collapse="\n\n"), con=inp)	
	return(inp)
}

smooth.phylo=function(phy, calibrations=NULL, base="", rm=TRUE){
#	calibrations: dataframe with minimally 'MRCA' 'MaxAge' 'MinAge' 'taxonA' and 'taxonB' 
#		-- if NULL, simple ultrametricization of 'phy' is performed
	
	phy$node.label=NULL
	if(!is.null(calibrations)){
		infile=write.pathd8(phy, calibrations, base)
	} else {
		infile=paste(base, "infile", sep=".")
		write.tree(phy, infile)
	}
	smooth.file=paste(base, "smoothed.tre", sep=".")
	parsed.outfile=paste(base, "pathd8.out", sep=".")
	outfile=paste(base, "pathd8.orig.out", sep=".")
	if(file.exists(outfile)) unlink(outfile)
	system(paste("pathd8 -n", infile, "-pn >", outfile, sep=" "))
	system(paste("grep \"d8 tree\" ", outfile, ">", parsed.outfile, sep=" "))
	smoothed=read.tree(parsed.outfile)
	if(rm) {
		unlink(parsed.outfile)
		unlink(smooth.file)
		unlink(outfile)
		unlink(infile)
	}
	return(smoothed)
}


