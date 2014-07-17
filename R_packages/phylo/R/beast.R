subset_beast=function(trees="", n=10, type="sequential"){
	x=readLines(trees)
	if(!length(grep("NEXUS", x, ignore.case=TRUE))) stop("'trees' does not appear to be a Nexus-formatted trees file.")
	idx=grep("STATE", x)
	save=switch(type, 
				sequential=idx[round(seq(1, length(idx), length=n))],
				random=idx[sample(1:length(idx), n, replace=FALSE)]
				)
	x=x[c(1:(min(idx)-1), save)]
	outfile=".tmp_TREES"
	writeLines(x, con=outfile)
	write("End;", file=outfile, append=TRUE)
	phy=read.nexus(outfile)
	unlink(outfile)
	return(phy)
}



taxa.beastblock=function(dat){
# rownames are species
# columns are characters
# NA is used for missing
# return string to be used as taxon block in BEAST xml for continuous data
	start="<taxa id=\"taxa\">"
	end="</taxa>"
	nn=paste("\"", names(dat), "\"", sep="")
	tt=lapply(1:nrow(dat), function(r) {
			  s=sapply(1:ncol(dat), function(c){
					   paste("\t\t<attr name=", nn[c], ">", dat[r,c], "</attr>", sep="")
					   })
			  taxstart=paste("\t", "<taxon id=", "\"", rownames(dat)[r], "\"", ">", sep="")
			  taxend="\t</taxon>"
			  paste(taxstart, paste(s, collapse="\n"), taxend, sep="\n")
			  })
	return(paste(start, paste(tt,collapse="\n"), end, collapse="\n", sep="\n"))
}

treemodel.beastblockQ=function(character=""){
# for quantitative data
	nn=paste(c("all", "root", "internal", "tip"), character, sep=".")
	types=list(c(1,1,1),c(1,0,0),c(0,1,0),c(0,0,1))
	names(types)=nn
	start="<treeModel id=\"treeModel\">"
	end="</treeModel>"
	tf=c("false", "true")
	char=quotef(character)
	tt=sapply(1:length(types), function(x){
			  typ=types[[x]]
			  typ=paste("\"", tf[match(typ, c(0,1))], "\"", sep="")
			  nodestart=paste(paste("\t", "<nodeTraits name=", char, sep=""), " rootNode=", typ[1], " internalNodes=", typ[2], " leafNodes=", typ[3], " traitDimension=1>", sep="")
			  nodeend="\t</nodeTraits>"
			  nt=paste("\"",names(types)[x],"\"",sep="")
			  paste(nodestart, paste("\t\t<parameter id=", nt, "/>", sep=""), nodeend, sep="\n")
			  })
	return(paste(start, paste(tt, collapse="\n"), end, sep="\n"))
}

quotef=function(char){paste("\"", char, "\"", sep="")}

startend.xml=function(parm, tabs=0){
	s=paste(ifelse(tabs>0, paste(rep("\t",tabs),collapse=""),""), "<",parm,">", sep="")
	e=paste(ifelse(tabs>0, paste(rep("\t",tabs),collapse=""),""), "</",parm,">", sep="")
	return(c(s,e))
}

parameter_value.xml=function(val, tabs=0){
	paste(ifelse(tabs>0, paste(rep("\t",tabs),collapse=""),""), "<parameter value=",quotef(val),"/>", sep="")
}

id.xml=function(parm, id, tabs=0){
	paste(ifelse(tabs>0, paste(rep("\t",tabs),collapse=""),""), "<", parm, " id=", id,">", sep="")
}

end.xml=function(parm, tabs=0){
	paste(ifelse(tabs>0, paste(rep("\t",tabs),collapse=""),""), "</", parm, ">", sep="")
}

start.xml=function(parm, tabs=0, eol=">", ...){
	paste(ifelse(tabs>0, paste(rep("\t",tabs),collapse=""),""), "<", parm, paste(..., sep=""), eol, sep="")
}


diffusion.beastblock=function(character="", root=0, prior=list(root=0.001, rate=2, wishart=100)){
	se.precisionP=startend.xml("precisionParameter",2)
	se.meanP=startend.xml("meanParameter",2)
	se.data=startend.xml("data",2)
	se.precisionM=startend.xml("precisionMatrix",2)
	se.matrixP=startend.xml("matrixParameter",3)
	se.scaleM=startend.xml("scaleMatrix",2)	
	
	## ROOT STATE 
	mvn_root=c(
			   start.xml("multivariateNormalPrior", 1, ">", paste(" id=", quotef(paste("rootprior",character,sep=".")), sep="")),
			   se.meanP[1],
			   parameter_value.xml(root, 3),
			   se.meanP[2],
			   se.precisionP[1],
			   se.matrixP[1],
			   parameter_value.xml(prior$root, 4),
			   se.matrixP[2],
			   se.precisionP[2],
			   se.data[1],
			   start.xml("parameter", 3, eol="/>", paste(" idref=", paste("root", character, sep="."), sep="")),
			   se.data[2],
			   end.xml("multivariateNormalPrior", 1)
			   )
	
	## DIFFUSION MODEL
	mvn_diffusion=c(
					start.xml("multivariateDiffusionModel", 1, ">", paste(" id=", quotef(paste("diffusionModel",character,sep=".")), sep="")),
					se.precisionM[1],
					se.matrixP[1],
					start.xml("parameter",4,eol="/>", paste(" id=", paste("prec", character, sep="."), " value=", quotef(prior$rate), sep="")),
					se.matrixP[2],
					se.precisionM[2],
					end.xml("multivariateDiffusionModel", 1)
					)
## WISHART PRIOR (Cauchy with 1 df)
## TRAIT LIKELIHOOD under DIFFUSION
}