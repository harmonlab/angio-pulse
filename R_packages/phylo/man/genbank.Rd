\name{fetch_gbtaxonomy.below}
\alias{fetch_gbtaxonomy.below}
\alias{fetch_gbtaxonomy.above}
\alias{collect_gbtaxonomy}
\alias{genbank}
\alias{ncbi}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fetching NCBI genbank taxonomy}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
	fetch_gbtaxonomy.below(taxon, rank = "genus", update.dmp = FALSE)

	fetch_gbtaxonomy.above(taxon, rank = "genus", update.dmp = FALSE)

	collect_gbtaxonomy(phy, rank = list(from = "family", to = "phylum"))

}
\details{
%- maybe also 'usage' for other objects documented here.
\itemize{
  \item{phy}{
%%     ~~Describe \code{phy} here~~
}
  \item{rank}{
%%     ~~Describe \code{rank} here~~
}
%- maybe also 'usage' for other objects documented here.
  \item{taxon}{
%%     ~~Describe \code{taxon} here~~
}
  \item{rank}{
%%     ~~Describe \code{rank} here~~
}
  \item{update.dmp}{
%%     ~~Describe \code{update.dmp} here~~
}
%- maybe also 'usage' for other objects documented here.
  \item{taxon}{
%%     ~~Describe \code{taxon} here~~
}
  \item{rank}{
%%     ~~Describe \code{rank} here~~
}
  \item{update.dmp}{
%%     ~~Describe \code{update.dmp} here~~
}
}
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(taxon, rank="genus", update.dmp=FALSE){
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
