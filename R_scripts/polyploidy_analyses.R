require(geiger)
require(phylo)


## USING MLE tree
setwd("~/Documents/tank/angiosperms/angio-pulse") # change to match local directory
dat=get(load("spermatophyta_AToL_639_PL_MEDUSA_BATCH.effectsize.rda"))$summary

# only consider upshifts
dat=dat[which(dat$r>0),]

# get the tree
phy=get(load("spermatophyta_AToL_639_PL_MEDUSA_BATCH.familial_MLE.rda"))$phy

shiftNode=match(rownames(dat), c(phy$tip.label, phy$node.label)) # nodes where shifts occur

# polyploidy events from tank
polyploidy=read.csv("data_files/Polyploidization.Dec2014.LJH.csv", header=TRUE)


# this returns all combinations of alternate placements for 
# the polyploidy events - there are 6 x 2 x 3 = 36 combinations
return_all_combinations<-function(polyploidy) {
	res<-list()
	counter<-1
	for(i in 3:8) # loop through WGD 6
		for(j in 11:12) # loop through WGD 6
			for(k in 15:17) { # loop through WGD 9 
				rowSelections<-c(1, 2, i, 9, 10, j, 13, 14, k)
				res[[counter]]<-polyploidy[rowSelections,]
				counter <- counter + 1
			}
	
	res			
}

# get them
wgdComb<-return_all_combinations(polyploidy)

# randomly assign n polyploidy events to branches in the tree -
# this is the null model
randomPolyploidyNodes<-function(phy, n) {
	sample(phy$edge[,2], size=n, replace=FALSE)
}

# test for exact matches
exactMatches<-numeric(length(wgdComb))
randomExactMatches<-matrix(nrow=1000, ncol=length(wgdComb))

# loop through all combinations
for(i in 1:length(wgdComb)) {
	
	# get the ith one
	pp<-wgdComb[[i]]
	
	# left and right tips used to get the nodes
	pl<-as.character(pp$left_tip)
	pr<-as.character(pp$right_tip)

	polyNode<-numeric(length=nrow(pp))
	for(j in 1:nrow(pp)) {
		if(pl[j]==pr[j]) { # this is a tip
			ww<-which(phy$tip.label==pl[j])
			xx<-which(phy$edge[,2]==ww)
			polyNode[j]<-phy$edge[xx,1]
		} else { # this is an internal node
			pair<-c(pl[j], pr[j])
			polyNode[j]<-getMRCA(phy, pair)
		}
	}
	
	exactMatches[i]<-sum(!is.na(match(polyNode, shiftNode))) # count exact matches
	
	for(j in 1:1000) { # repeat for 1000 random assignments of polyploidy
		randomShifts<-randomPolyploidyNodes(phy, length(shiftNode))
		randomExactMatches[j,i] <- sum(!is.na(match(polyNode, randomShifts)))
	}
}

# calculate p-values
pExact<-numeric(length(exactMatches)) 
for(i in 1:length(exactMatches)) {
	pExact[i]<-(sum(randomExactMatches[,i] >= exactMatches[i])+1)/1001
}

# six are significant
pExact < 0.05
	
# now repeat using downstream distances

# test for exact matches
tipDistances<-numeric(length(wgdComb))
randomTipDistances <-matrix(nrow=1000, ncol=length(wgdComb))

#tcache<-geiger:::.cache.descendants(phy)

library(phytools)
p2<-phy
p2$edge.length[]<-1
lookupDistance<-dist.nodes(p2)

closeCount<-function(pn, sn, cutoff=3) {
	
	mdists<-numeric(length(pn))
	# gather descendants of poly nodes pn
	for(i in 1:length(pn)) {
		pd<-getDescendants(phy, pn[i])
		targets<-pd[which(pd %in% sn)]
		targetDistance<-lookupDistance[pn[i],targets]
		if(length(targetDistance)!=0) {
			mdists[i]<-min(targetDistance)
		} else mdists[i]<-Inf
	}
	testStat<-sum(mdists <= cutoff)
	testStat
}

# test for close matches
testStat<-numeric(length(wgdComb))
nullDist<-matrix(nrow=1000, ncol=length(wgdComb))

# loop through all combinations
for(i in 17:length(wgdComb)) {
	
	# get the ith one
	pp<-wgdComb[[i]]
	
	# left and right tips used to get the nodes
	pl<-as.character(pp$left_tip)
	pr<-as.character(pp$right_tip)

	polyNode<-numeric(length=nrow(pp))
	for(j in 1:nrow(pp)) {
		if(pl[j]==pr[j]) { # this is a tip
			ww<-which(phy$tip.label==pl[j])
			xx<-which(phy$edge[,2]==ww)
			polyNode[j]<-phy$edge[xx,1]
		} else { # this is an internal node
			pair<-c(pl[j], pr[j])
			polyNode[j]<-getMRCA(phy, pair)
		}
	}
	
	testStat[i]<-closeCount(polyNode, shiftNode)
		
	for(j in 1:1000) { # repeat for 1000 random assignments of polyploidy
		randomShifts<-randomPolyploidyNodes(phy, length(shiftNode))
		nullDist[j,i] <- closeCount(polyNode, randomShifts)
	}
}

# calculate p-values
pClose<-numeric(length(testStat)) 
for(i in 1:length(testStat)) {
	pClose[i]<-(sum(nullDist[,i] >= testStat[i])+1)/1001
}

