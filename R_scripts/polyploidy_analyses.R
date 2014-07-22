require(geiger)

assess_nearness_in_rate_shift=function(phy, nodes, tips=TRUE){
	cache=geiger:::.cache_tree(phy)
	N=Ntip(phy)
	n=Nnode(phy)
	nn=N+n
	others=(1:nn)[-(N+1)]
    if(!tips) others=others[-c(1:N)]
	dD=function(node){
		if(node<=N) return(NA)
		count=0
		
		zfun=function(n){
			any(n%in%nodes)
		}
		
		dd=node
		while(1){
			if(zfun(dd)) break()
			count=count+1
			dd=dd[dd>N]
			if(!length(dd)){
				count=NA
				break()
			}
			dd=unlist(cache$fdesc[dd])
		}
		return(count)
	}
	
	aA=function(node){
		anc=cache$anc[[node]]
		count=0
		if(node%in%nodes) return(count)
		count=NA
		if(any(anc%in%nodes)){
			count=-min(match(nodes, anc), na.rm=TRUE)
		}
		return(count)
	}

	tmp=lapply(others, function(x){
        if(x%in%nodes) return(0)
		aa=aA(x)
		dd=dD(x)
		f=c(aa, dd)
		f=f[!is.na(f)]
		if(!length(f)) return(NA)
		if(length(unique(f))==1) return(f)
		z=abs(c(dd,aa))
		return(f[which(z==min(z))])
	})
	
    ndists=unlist(tmp)
	ndists=ndists[!is.na(ndists)]
	return(ndists)
}

## PLOT node distance from random nodes to nearest shift against polyploidization nodes from nearest shift
## USING MLE tree
dat=get(load("~/Dropbox/angiospermsnescent/analyses/medusa/macgyver_turbo_0.92-3/results/spermatophyta_AToL_639_PL_MEDUSA_BATCH.effectsize.rda"))$summary
dat=dat[which(dat$r>0),]
phy=get(load("~/Dropbox/angiospermsnescent/analyses/medusa/macgyver_turbo_0.92-3/results/spermatophyta_AToL_639_PL_MEDUSA_BATCH.familial_MLE.rda"))$phy

nodes=match(rownames(dat), c(phy$tip.label, phy$node.label))

rand_ndist=assess_nearness_in_rate_shift(phy, nodes)

polyploidy=read.csv("polyploidization_events.csv", header=TRUE)
p_ndist=polyploidy$ndist
p_ndist=p_ndist[!is.na(p_ndist)]

## INCLUDING TIPS
cat("INCLUDING TIP LINEAGES\n\n")
## BINOMIAL TESTS
# perfect correspondence
x=sum(p_ndist==0)
n=length(p_ndist)
binom.test(x, n, p=sum(rand_ndist==0)/length(rand_ndist))

# delayed correspondence
x=sum(p_ndist>0)
n=length(p_ndist)
binom.test(x, n, p=sum(rand_ndist>0)/length(rand_ndist))

## EXCLUDING TIPS
cat("EXCLUDING TIP LINEAGES\n\n")
rand_tdist=assess_nearness_in_rate_shift(phy, nodes, tips=FALSE)

t_ndist=polyploidy$ndist
t_ndist=t_ndist[!is.na(t_ndist) & !polyploidy$tip]
x=sum(t_ndist>0)
n=length(t_ndist)
binom.test(x, n, p=sum(rand_tdist>0)/length(rand_tdist))


pdf("internode_distance_to_nearest_shift.pdf")
hist(p_ndist, breaks=c(-4,-2,0,2,4), include.lowest=TRUE, col=geiger:::.transparency("black", 0.85), freq=FALSE, xlim=c(-20,20), xlab="internode distance", ylab="density", main="")
hist(rand_ndist, breaks=10,  col=geiger:::.transparency("lightgray", 1), freq=FALSE, xlim=c(-20,20), xlab="internode distance", ylab="density", main="", add=TRUE)
legend("topright", legend=c("polyploidization", "random"), pt.bg=c(geiger:::.transparency("black", 0.85), geiger:::.transparency("lightgray", 1)), pch=22, bty="n")
dev.off()


pdf("internode_distance_to_nearest_shift -- excluding tips.pdf")
hist(t_ndist, breaks=c(-4,-2,0,2,4), include.lowest=TRUE, col=geiger:::.transparency("black", 0.85), freq=FALSE, xlim=c(-20,20), xlab="internode distance", ylab="density", main="")
hist(rand_tdist, breaks=10,  col=geiger:::.transparency("lightgray", 1), freq=FALSE, xlim=c(-20,20), xlab="internode distance", ylab="density", main="", add=TRUE)
legend("topright", legend=c("polyploidization", "random"), pt.bg=c(geiger:::.transparency("black", 0.85), geiger:::.transparency("lightgray", 1)), pch=22, bty="n")
dev.off()
