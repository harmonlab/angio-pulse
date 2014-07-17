\name{glomogram.phylo}
\alias{glomogram.phylo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Grafting subtrees to a backbone}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
glomogram.phylo(phy, subtrees)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
%%     ~~Describe \code{phy} here~~
}
  \item{subtrees}{
%%     ~~Describe \code{subtrees} here~~
}
}
\details{
# subtrees: a named list of subtrees (each can be a 'multiPhylo' object)
#		-- minimally has 'subtree' and 'age' elements
#		-- if 'tip' is missing, the name of the 'subtrees' object is taken as 'tip'
# phy: the skeleton tree with tip labels corresponding to names in 'subtrees'
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
	
choose.tip=function(phy, subtree){
	while(1){
		split_age=runif(1, 0, max(phy$edge.length[phy$edge[,2]<=Ntip(phy)]))
		subtree=scaletree(subtree, split_age)
		tipedges=phy$edge[,2]<=Ntip(phy)
		tt=phy$edge.length[tipedges]
		names(tt)=phy$tip.label[phy$edge[tipedges,2]]
		tt=tt[tt>split_age]
		if(length(tt)){
			tip=names(tt[sample(1:length(tt),1)])
			break()
		}
	}
	return(list(age=split_age, tip=tip, subtree=subtree))	
}

phy<-orig<-scaletree(rcoal(80),100)
st=3
subtrees=lapply(1:st, function(x) {y=rcoal(sample(3:100, 1)); y$tip.label=paste(letters[x], 1:Ntip(y), sep="_"); choose.tip(phy, y)})
names(subtrees)=letters[1:st]
phy=glomogram.phylo(phy, subtrees)
plot.phylo(phy, cex=0.5, type="f", tip.color=ifelse(is.element(phy$tip.label,orig$tip.label),"black","red"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
