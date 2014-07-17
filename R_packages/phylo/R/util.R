fix.nonASCII=function(dat){
	sub.nonASCII=function(string){
		require(tools)
		ss=unlist(strsplit(string,""))
		sink(file=".tmp") 
		test=showNonASCII(ss)
		sink()
		unlink(".tmp")
		ss=paste(ss[!ss%in%test],collapse="")
		ss
		
	}
	
	dat=as.matrix(dat)
	sink(file=".tmp") 
	nn=lapply(1:ncol(dat), function(x) {v=showNonASCII(dat[,x]); if(!length(v)) return(NA) else return(v)})
	sink()
	unlink(".tmp")
	for(i in 1:length(nn)){
		sub=unname(nn[[i]])
		if(!all(is.na(sub))){
			dat[,i]<-tmp<-unname(sapply(dat[,i], function(cur) if(cur%in%sub) return(sub.nonASCII(cur)) else return(cur)))
		}
	}
	dat
}

add.transparency <- function (col, alpha) 
{
    tmp <- col2rgb(col)/255
    rgb(tmp[1, ], tmp[2, ], tmp[3, ], alpha = alpha)
}

withinrange <- function (x, min, max) 
{
    a = sign(x - min)
    b = sign(x - max)
    if (abs(a + b) == 2) 
	return(FALSE)
    else return(TRUE)
}

basename.noext=function(path=""){
	return(sub("[.][^.]*$", "", basename(path), perl=TRUE))
	
}


grep.list=function(dat, term="", strict=FALSE, ...){
# dat: list of character vectors
	res=sapply(dat, function(d) {
			   g=sapply(term, function(t) {
						if(strict) {
						grep(t, d, ignore.case=FALSE, ...)
						} else {
						grep(paste(" ", t, " ", sep=""), d, ignore.case=TRUE, ...)
						}
						})
			   d[unique(unlist(g))]
			   })
	names(res)=names(dat)
	res
}	

scat <- function(files, outfile){
	system(paste("cat", paste(files, collapse=" "), ">>", outfile, sep=" "))
}

is.empty=function(file=""){
	length(scan(file, quiet=TRUE, what="char", nlines=1))==0
}

sem=function(val){
	sd(val, na.rm=TRUE)/sqrt(length(val))
}

resolve.executable=function(){
	packagedir=system.file(package="phylo")
	execs=lapply(d<-dir(paste(packagedir,"exec",sep="/")), function(x) {paste(packagedir, "exec", x, sep="/")})
	names(execs)=basename.noext(d)
	return(execs)
}

