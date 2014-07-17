\name{congruify.phylo}
\alias{congruify.phylo}
\alias{nodelookup.phylo}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ultrametricization of megaphylogenies from a backbone timetree}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
congruify.phylo(backbone, megaphylogeny, taxonomy = NULL, tol = 0.1)
nodelookup.phylo(phy, label)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{backbone}{
a 'donor' tree of dates; used to time-scale the 'megaphylgeny'; backbone should be ultrametric
}
  \item{megaphylogeny}{
a phylogram that is sought to be ultrametricized based on the backbone
}
  \item{phy}{
a phylogenetic tree in 'phylo' format
}

  \item{taxonomy}{
a linkage table between tips of the phylogeny and clades represented in the tree; rownames of 'taxonomy' should be tips found in the phylogeny
}
  \item{tol}{
branching time in 'backbone' above which constraints will be applied for ultrametricization of 'megaphylogeny'
}
  \item{label}{
a particular node label for which the numeric identifier is sought ('label' must appear in the 'node.label' object of the phylogeny)
}

}
\details{
#	backbone: a time-calibrated phylogeny with tip-labels that can be treated as an exemplar for clades in 'supertree'
#		-- e.g., tip.label in 'backbone' might be "Pinus" while in 'supertree' we might have "Pinus_cembra"
#		-- tips in 'backbone' can be contained in 'supertree' (FIXME: is this true?)
#	megaphylogeny: a rooted phylogeny that is to be time-scaled based on 'backbone'
#	taxonomy: linkage table between tipsets for 'backbone' and 'supertree'; if empty, one is attempted to be build by 'supertree' labels
#		-- if NULL, 'backbone' tips must correspond to tips in 'supertree'... e.g., A, B, C in 'backbone'; A_sp1, B_sp2, C_sp3 in 'supertree'
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
	data(AToL)
	backbone=AToL$backbone
	megatree=AToL$megaphylogeny

	# construct a lookup table for the megaphylogeny
	tax=get(data(spermatophyta_APGIII_taxonomy))
	megatree_tax=matrix(NA, nrow=Ntip(megatree), ncol=3)
	gen=sapply(megatree$tip.label, function(x) split_lab(x)[1])
	megatree_tax=tax[match(gen, tax$genus),]
	rownames(megatree_tax)=megatree$tip.label
	cld=get(data(spermatophyta_APGIII_clades))
	lkp=build_lookup.phylo(megatree, megatree_tax, cld)

	# label nodes in megaphylogeny based on its lookup table
	reduced_lkp=lkp[-which(is.element(names(lkp),c("species","genus","family")))]
	megatree_labeled=nodelabel.phylo(megatree, reduced_lkp)

	# extract gunnerid subtree
	gunnerid_nd=nodelookup.phylo(megatree_labeled, "Gunneridae")
	gunnerid_subtree=ladderize(extract.clade(megatree, gunnerid_nd), right=FALSE)
	plot.phylo(gunnerid_subtree, show.tip=FALSE, type="f", edge.width=0.1)

	# congruify gunnerid subtree with pathd8
	tmp=congruify.phylo(backbone, gunnerid_subtree, taxonomy=lkp)
	gunnerid_timetree=tmp[[1]]$phy
	dev.new()
	plot.phylo(gunnerid_timetree, show.tip=FALSE, type="f", edge.width=0.1)
	
	# create node-labeled family-level phylogeny from the backbone tree
	backbone_family=rank_prune.phylo(backbone, tax, rank="family")
	tax_fam=unique(tax[,-1])
	backbone_tax=matrix(NA, nrow=Ntip(backbone_family), ncol=2)
	backbone_tax=tax_fam[match(backbone_family$tip.label, tax_fam$family),]
	rownames(backbone_tax)=backbone_family$tip.label	
	lkp_backbone=build_lookup.phylo(backbone_family, backbone_tax, cld)
	backbone_labeled=nodelabel.phylo(backbone_family, lkp_backbone)
	dev.new()
	plot.phylo(backbone_labeled, edge.width=0.1, type="f", cex=0.1)
	nodelabels(backbone_labeled$node.label, col="red", cex=0.5, frame="n")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
