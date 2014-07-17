\name{build_phylo.lookup}
\alias{build_phylo.lookup}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Constructing a phylogenetic tree from taxonomic classification}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
build_phylo.lookup(lookup, tips, rank = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{taxonomy}{
%%     ~~Describe \code{taxonomy} here~~
}
  \item{labels}{
%%     ~~Describe \code{labels} here~~
}
  \item{rank}{
%%     ~~Describe \code{rank} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
tmp=get(data(spermatophyta_APGIII_orders))
dat=attributes(tmp)$clades
phy=build_phylo.lookup(dat[,-ncol(dat)], rownames(dat))$phy
plot(ladderize(phy, right=FALSE), type="fan", cex=0.5, edge.width=0.1)

## AN EXAMPLE of building a full family-level topology for seed plants
orders=get(data(spermatophyta_APGIII_orders))
orders.tax=attributes(orders)$clades
tmp=build_phylo.lookup(orders.tax, rownames(orders.tax))
orders.phy=nodelabel.phylo(tmp$phy, orders.tax)

brlen=min(orders.phy$edge.length[which(orders.phy$edge[,2]<=Ntip(orders.phy))])
family.phy=lapply(1:length(orders), function(idx){
	oo=orders[[idx]]
	return(list(age=brlen/2, tip=names(orders)[idx], subtree=polytomy.phylo(unlist(oo), brlen/2)))
})
names(family.phy)=names(orders)
if(length(ww<-which(names(family.phy)=="UNPLACED"))){
	unplaced=family.phy[[ww]]
	family.phy=family.phy[-ww]
} else {
	unplaced=NULL
}

phy=orders.phy
for(f in 1:length(family.phy)){
	cur_family=family.phy[[f]]
	if(Ntip(cur_family$subtree)==1){
		phy$tip.label[which(phy$tip.label==cur_family$tip)]=unname(cur_family$subtree$tip.label)
	} else {
		phy=glomogram.phylo(phy, family.phy[f])
	}
}

tmp=unlist(unname(sapply(family.phy, function(x) {y=x$subtree$tip.label; names(y)=rep(x$tip, length(y)); y})))
family.tax=data.frame(family=tmp, order=names(tmp), stringsAsFactors=FALSE)
phy=compute.brlen(phy, method="Grafen")
cld=get(data(spermatophyta_APGIII_clades))
lkp_tmp=build_lookup.phylo(phy, family.tax, cld)
lkp=data.frame(as.matrix(lkp_tmp[,-1]), stringsAsFactors=FALSE)
phy=nodelabel.phylo(phy, lkp)

# place unplaced taxa
for(i in 1:length(cld)){
	cur_cld=cld[[i]]
	if(!is.null(cur_tips<-attributes(cur_cld)$unplaced)){
		label=names(cld)[i]
		nd=nodelookup.phylo(phy, label)
		pl=polytomy.phylo(cur_tips, brlen/2)
		phy=bind.tree(phy, pl, where=nd)
	}
}

phy=compute.brlen(ladderize(phy, right=FALSE), method="Grafen")
tax=matrix("", Ntip(phy), Nnode(phy))
rownames(tax)<-tips<-phy$tip.label
nodes=(Ntip(phy)+1):(Ntip(phy)+Nnode(phy))
labels=phy$node.label
for(idx in 1:Nnode(phy)){
	nd=nodes[idx]
	dd=tips[get.descendants.of.node(nd, phy, tips=TRUE)]
	tax[is.element(rownames(tax),dd),idx]=labels[idx]
}
colnames(tax)=labels
tax=tax[,ncol(tax):1]
spermatophyta_APGIII_families=list(phy=phy, clades=tax)
plot(ladderize(phy, right=FALSE), type="f", cex=0.25, show.node=TRUE, edge.width=0.1, tip.col="gray")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
