\name{add_fossils.phylo}
\alias{add_fossils.phylo}
\alias{add.fossiltip}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Grafting (fossil) taxa to phylogenetic trees}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
add_fossils.phylo(phy, fossils.df)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
%%     ~~Describe \code{phy} here~~
}
  \item{fossils.df}{
%%     ~~Describe \code{fossils.df} here~~
}
}
\details{
	
# fossils.df: data.frame of defining taxa, fossil_age, and name...
#	   t1 t2  time      name
#	1 t42 t53   10 fossil_a
#	2 t33 t34   10 fossil_b
#	3 t31 t55   10 fossil_c
#	4 t48 t38   10 fossil_d
#	5 t54 t55   10 fossil_e
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
scaletree=function(phy, time){
	depth = max(diag(vcv.phylo(phy)))
    phy$edge.length = (phy$edge.length/depth) * time
    phy
}

# generates an example data.frame to be used with add_fossils.phylo()
example.df=function(nodes, phy, tol=1){
	
# nodes: number of fossils to add
# tol: defines threshold internal edge to which to bind fossil
	
	times=get.btimes(phy)
	tips=phy$tip.label
	df=matrix(NA, nrow=nodes, ncol=2)
	split=numeric(nodes)
	anc=numeric(nodes)
	root=Ntip(phy)+1
	i=0
	while(i<nodes){
		tt=sample(tips, 2)
		aa=getMRCA(phy, match(tt,tips))
		a=times$anc[which(times$des==aa)]
		start=times$time[which(times$des==a)]
		end=times$time[which(times$des==aa)]
		if(aa>root) {
			if((start-end)>tol){
				i=i+1
				df[i,1:2]=tt
				split[i]=runif(1, 0, start)
				anc[i]=aa
			}
		}
	}
	df=data.frame(df, stringsAsFactors=FALSE)
	names(df)=c("t1", "t2")
	df$time=split
	df$ancestor=anc
	df$name=paste("fossil", letters[1:nrow(df)], sep="_")
	return(df)
}


target=20
fossils=5
while(1){
	phy=try(scaletree(rcoal(target),100),silent=TRUE)
	if(!inherits(phy, "try-error")) if(Ntip(phy)==target) break()
}
orig=phy
phy$tip.label=paste("t",1:Ntip(phy),sep="")
fossils.df=example.df(fossils,phy)
system.time(phy<-add_fossils.phylo(phy, fossils.df))
plot(ladderize(phy),cex=0.5,lwd=0.1,edge.color="gray")
mtext("assembled fossils")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
