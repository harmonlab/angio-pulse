require(phylo)

syn=get(data(spermatophyta_APGIII_families))$synonymy
ord=get(data(spermatophyta_APGIII_orders))
tmp=unlist(lapply(1:length(ord), function(x){y=ord[[x]]; names(y)=rep(names(ord)[x], length(y)); y}))
oo=data.frame(family=tmp, order=names(tmp), stringsAsFactors=FALSE)
oo=oo[oo$order!="UNPLACED",]
rownames(oo)=oo$family
pl=get(data(plantlist_classification))
famr=split(pl, pl$family)
frichnesses=data.frame(n.taxa=sapply(famr, nrow))
#frichnesses[frichness
swap=rownames(frichnesses)%in%syn$from
if(any(swap)){
	rownames(frichnesses)[swap]=syn$to[match(rownames(frichnesses)[swap], syn$from)]
}
tmp=cbind(oo, n.taxa=frichnesses[match(oo$family, rownames(frichnesses)),"n.taxa"])
ordr=split(tmp, tmp$order)
orichnesses=data.frame(n.taxa=sapply(ordr, function(x) sum(x$n.taxa, na.rm=TRUE)))

phy=read.tree("out_dates.tre")
drop_genera=list(Escalloniaceae="Platyspermation", Olacaceae="Ximenia", Phytolaccaceae="Rivina", Gesneriaceae="Saintpaulia", Aristolochiaceae="Saruma")
phy=drop.tip(phy, unlist(drop_genera))

tax=read.csv("spermatophyta_AToL_639_PL_MEDUSA_BATCH.MLE_0.92.phy-taxonomy.csv",as.is=TRUE,row=1)
ordphy=rank_prune.phylo(phy, tax, "order")

prn=function(phy, dat){
	common=intersect(rownames(dat),phy$tip.label)
	prnphy=drop.tip(phy, phy$tip.label[!phy$tip.label%in%common])
	prndat=dat[which(rownames(dat)%in%common), ]
	names(prndat)=rownames(dat)[rownames(dat)%in%common]
	cld=lapply(1:length(prndat), function(x) paste(names(prndat)[x],1:prndat[x], sep="_") )
	names(cld)=names(prndat)
	prnphy$clades=cld
	class(prnphy)=c("clade.tree", class(prnphy))
	prnphy
}

richnesses_and_ages.phylo=function(phy){
	if(!any("clade.tree"%in%class(phy))) stop("supply 'phy' as a cladetree")
	if(!is.ultrametric(phy)) {
		phy=ultrametricize.phylo(phy)
	}
	rr=sapply(phy$clades, length)
	ww=phy$edge[,2]<=Ntip(phy)
	aa=phy$edge.length[ww]
	names(aa)=phy$tip.label[phy$edge[ww,2]]
	dat=merge(data.frame(rr), data.frame(aa), by=0)
	rownames(dat)=dat$Row.names
	dat=dat[,-1]
	names(dat)=c("richness","age")
	dat
}


#prnord=prn(ordphy, orichnesses)
prnfam=get(load("spermatophyta_AToL_639_PL_MEDUSA_BATCH.familial_MLE.rda"))$phy


#orddat=richnesses_and_ages.phylo(prnord)
#olm=lm(log(orddat$richness)~orddat$age)

famdat=richnesses_and_ages.phylo(prnfam)
flm=lm(log(famdat$richness)~famdat$age)

gety=function(m,b,x){
	m*x+b
}

pdf("richness_age.spermatophyta.pdf", width=6, height=12)
#layout(matrix(1:2, nrow=2))
xx=range(famdat$age)
yy=gety(coef(flm)[2], coef(flm)[1], xx)
plot(famdat$richness~famdat$age,  log="y", xlab="stem age",  ylab="species richness", bty="n", pch=19, xlim=c(0, max(famdat$age)))
mtext("families")
lines(xx,exp(yy),lty=2)

#xx=range(orddat$age)
#yy=gety(coef(olm)[2], coef(olm)[1], xx)
#plot(orddat$richness~orddat$age,  log="y", xlab="stem age",  ylab="species richness", bty="n", pch=19, xlim=c(0, max(orddat$age)))
#mtext("orders")
#abline(a=coef(olm)[1], b=coef(olm)[2], lty=2)

#lines(xx,exp(yy),lty=2)

dev.off()
