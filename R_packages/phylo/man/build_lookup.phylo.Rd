\name{build_lookup.phylo}
\alias{build_lookup.phylo}
\alias{nodelabel.phylo}
\alias{rank_prune.phylo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Constructing a 'lookup' table for tips in a phylogeny based on clade definitions and a nomenclatural table}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
build_lookup.phylo(phy, taxonomy, clades = NULL)
nodelabel.phylo(phy, taxonomy, rank = NULL, strict = FALSE)
rank_prune.phylo(phy, taxonomy, rank = "family")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
a phylogenetic tree in 'phylo' format
}
  \item{taxonomy}{
a linkage table between tips of the phylogeny and clades represented in the tree; rownames of 'taxonomy' should be tips found in the phylogeny
}
  \item{clades}{
a named list of clade definitions (i.e., spanning taxa; see Examples); spanning taxa may invoke other definitions found within the 'clades' list
}
  \item{rank}{
a clade level of interest ('rank' must appear in the column names of 'taxonomy')
}
}
\details{
## taxonomy expected to have first column at same level as tip labels in phy
## first row in taxonomy is most exclusive
## clade_defs are phylogenetic trees of a clade representation

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
	orders=get(data(spermatophyta_APGIII_orders))
	families=get(data(spermatophyta_APGIII_families))
	cld=get(data(spermatophyta_APGIII_clades))
	print(cld)
	tmp=unlist(unname(sapply(1:length(orders), function(idx) {cur=orders[[idx]]; names(cur)=rep(names(orders)[idx], length(cur)); cur})))
	tax=data.frame(family=tmp, order=names(tmp), stringsAsFactors=FALSE)
	tax[tax=="UNPLACED"]=""
	
	phy=families$phy
	lookup=build_lookup.phylo(phy, tax, cld)
	label_phy=nodelabel.phylo(phy, lookup)
	plot.phylo(label_phy, show.tip=TRUE, show.node.label=FALSE, type="f", cex=0.15, edge.width=0.15, tip.color="gray", label.offset=0.025)
	nodelabels(label_phy$node.label, frame="n", cex=0.5, col="red")
	
	rownames(tax)=tax[,"family"]
	pruned=rank_prune.phylo(phy, tax, rank="order")
	otax=attributes(orders)$clades
	olabel_phy=nodelabel.phylo(pruned, otax)
	plot.phylo(olabel_phy, cex=0.5, type="f", show.node=FALSE, edge.width=0.1)
	nodelabels(olabel_phy$node.label, col="red", cex=0.5, frame="n")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
