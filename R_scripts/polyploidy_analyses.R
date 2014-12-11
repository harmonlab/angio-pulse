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

shiftNode<-numeric(dim(dat)[1])

for(i in 1:dim(dat)[1]) {
	nameToMatch<-rownames(dat)[i]
	
	if(nameToMatch %in% phy$node.label) {
		theMatch<-which(phy$node.label==nameToMatch)
		shiftNode[i]<-theMatch + length(phy$tip.label) # aligns phy$node.label with numbers in edge matrix
	} else {
		theMatch<-which(phy$tip.label==nameToMatch)
		shiftNode[i]<-theMatch
	}
	
}



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
			polyNode[j]<-ww
		} else { # this is an internal node
			pair<-c(pl[j], pr[j])
			polyNode[j]<-getMRCA(phy, pair)
		}
	}
	
	exactMatches[i]<-sum(!is.na(match(polyNode, shiftNode))) # count exact matches
	
	for(j in 1:1000) { # repeat for 1000 random assignments of polyploidy
		randomPolys<-randomPolyploidyNodes(phy, length(polyNode))
		randomExactMatches[j,i] <- sum(!is.na(match(shiftNode, randomPolys)))
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

tipDistances<-numeric(length(wgdComb))
randomTipDistances <-matrix(nrow=1000, ncol=length(wgdComb))


library(phytools)
p2<-phy
p2$edge.length[]<-1
lookupDistance<-dist.nodes(p2)

# eliminate tip node WDGs
polyploidy2<-polyploidy[c(-9, -13, -14),]

# this returns all combinations of alternate placements for 
# the polyploidy events - there are 6 x 2 x 3 = 36 combinations
# this one deals with polyploidy2, which has tip WGDs left out
return_all_combinations2<-function(polyploidy2) {
	res<-list()
	counter<-1
	for(i in 3:8) # loop through WGD 6
		for(j in 10:11) # loop through WGD 6
			for(k in 12:14) { # loop through WGD 9 
				rowSelections<-c(1, 2, i, 9, j, k)
				res[[counter]]<-polyploidy2[rowSelections,]
				counter <- counter + 1
			}
	
	res			
}

# get them
wgdComb2<-return_all_combinations2(polyploidy2)

closeCount<-function(pn, sn, cutoff=3, print=F) {
	
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
	if(print) cat(pn[which(mdists <= cutoff)], "\n")
	testStat
}

# test for close matches
testStat<-numeric(length(wgdComb2))
nullDist<-matrix(nrow=1000, ncol=length(wgdComb2))

# loop through all combinations
for(i in 1:length(wgdComb2)) {
	
	# get the ith one
	pp<-wgdComb2[[i]]
	
	# left and right tips used to get the nodes
	pl<-as.character(pp$left_tip)
	pr<-as.character(pp$right_tip)

	polyNode<-numeric(length=nrow(pp))
	for(j in 1:nrow(pp)) {
		if(pl[j]==pr[j]) { # this is a tip
			ww<-which(phy$tip.label==pl[j])
			polyNode[j]<-ww
		} else { # this is an internal node
			pair<-c(pl[j], pr[j])
			polyNode[j]<-getMRCA(phy, pair)
		}
	}
	
	testStat[i]<-closeCount(polyNode, shiftNode, print=T)
		
	for(j in 1:1000) { # repeat for 1000 random assignments of polyploidy
		randomPolys<-randomPolyploidyNodes(phy, length(polyNode))
		nullDist[j,i] <- closeCount(randomPolys, shiftNode)
	}
	cat(i, "\n")
}

# calculate p-values
pClose<-numeric(length(testStat)) 
for(i in 1:length(testStat)) {
	pClose[i]<-(sum(nullDist[,i] >= testStat[i])+1)/1001
}

pClose


# retest with different cutoff values

allPValues<-matrix(nrow=5, ncol=length(wgdComb2))
for(ct in 1:5) {
	testStat<-numeric(length(wgdComb2))
	nullDist<-matrix(nrow=1000, ncol=length(wgdComb2))

	# loop through all combinations
	for(i in 1:length(wgdComb2)) {
	
	# get the ith one
	pp<-wgdComb2[[i]]
	
	# left and right tips used to get the nodes
	pl<-as.character(pp$left_tip)
	pr<-as.character(pp$right_tip)

	polyNode<-numeric(length=nrow(pp))
	for(j in 1:nrow(pp)) {
		if(pl[j]==pr[j]) { # this is a tip
			ww<-which(phy$tip.label==pl[j])
			polyNode[j]<-ww
		} else { # this is an internal node
			pair<-c(pl[j], pr[j])
			polyNode[j]<-getMRCA(phy, pair)
		}
	}
	
	testStat[i]<-closeCount(polyNode, shiftNode, cutoff= ct)
		
	for(j in 1:1000) { # repeat for 1000 random assignments of polyploidy
		randomPolys<-randomPolyploidyNodes(phy, length(polyNode))
		nullDist[j,i] <- closeCount(randomPolys, shiftNode, cutoff= ct)
	}
	cat(i, " ")
	}

	# calculate p-values
	for(i in 1:length(testStat)) {
		allPValues[ct,i]<-(sum(nullDist[,i] >= testStat[i])+1)/1001
	}
	cat("\n", ct, " is done\n")
}

rownames(allPValues)<-1:5
colnames(allPValues)<-paste("set", 1:36)
write.csv(allPValues, file="allThePValues.csv")

# time distance calculations for table

# these are the crown distances
lookupTimeDistance<-dist.nodes(phy)
p3<-phy
p3$node.label<-NULL
lookupBt<-branching.times(p3)


getBranchPairDistance<-function(phy, node1, node2) {
	# crown-to-crown: that's easy
	midDist<-lookupTimeDistance[node1, node2]
	
	age1<-lookupBt[which(names(lookupBt)==node1)]
	age2<-lookupBt[which(names(lookupBt)==node2)]
	
	#fixing problem when these are tips
	if(length(age1)==0) age1<-0
	if(length(age2)==0) age2<-0
	
	ww<-which(phy$edge[,2]==node1)
	n1anc<-phy$edge[ww,1]
	ww<-which(phy$edge[,2]==node2)
	n2anc<-phy$edge[ww,1]

	ancAge1<-lookupBt[which(names(lookupBt)==n1anc)]
	ancAge2<-lookupBt[which(names(lookupBt)==n2anc)]

	if(age1 > age2) {
		oldCrown<-age1
		youngCrown<-age2
		oldStem<-ancAge1
		youngStem<-ancAge2
	} else {
		oldCrown<-age2
		youngCrown<-age1
		oldStem<-ancAge2
		youngStem<-ancAge1
	}
	
	
	# minumum possible: old crown to young stem
	minDist<-oldCrown - youngStem
	
	# maximum possible: old stem to young crown
	maxDist<-oldStem - youngCrown
	
	res<-c(midDist, minDist, maxDist)
	names(res)<- c("midDist", "minDist", "maxDist")	
	res
}


getNodePairDistance<-function(phy, node1, node2) {
	
	nodeDist<-lookupDistance[node1, node2]
	
	age1<-lookupBt[which(names(lookupBt)==node1)]
	age2<-lookupBt[which(names(lookupBt)==node2)]
	
	#fixing problem when these are tips
	if(length(age1)==0) age1<-0
	if(length(age2)==0) age2<-0
	
	ww<-which(phy$edge[,2]==node1)
	n1anc<-phy$edge[ww,1]
	ww<-which(phy$edge[,2]==node2)
	n2anc<-phy$edge[ww,1]

	ancAge1<-lookupBt[which(names(lookupBt)==n1anc)]
	ancAge2<-lookupBt[which(names(lookupBt)==n2anc)]

	if(age1 > age2) {
		oldCrown<-age1
		youngCrown<-age2
		oldStem<-ancAge1
		youngStem<-ancAge2
	} else {
		oldCrown<-age2
		youngCrown<-age1
		oldStem<-ancAge2
		youngStem<-ancAge1
	}
	
	
	# minumum possible: old crown to young stem
	minDist<-oldCrown - youngStem
	
	# maximum possible: old stem to young crown
	maxDist<-oldStem - youngCrown
	
	res<-c(midDist, minDist, maxDist)
	names(res)<- c("midDist", "minDist", "maxDist")	
	res
}

# test
gg<-getBranchPairDistance(phy, 379, 381)

allDistances<-matrix(nrow=dim(polyploidy)[1], ncol=5)
colnames(allDistances)<- c("match","midDist", "minDist", "maxDist", "nodeDist")	

findClosestShift<-function(pn, snodes) {
	
	pd<-getDescendants(phy, pn)
	targets<-pd[which(pd %in% sn)]
	targetDistance<-lookupTimeDistance[pn,targets]
	if(length(targetDistance)!=0) {
			closeDist<-min(targetDistance)
			names(closeDist)<-targets[which(targetDistance==min(targetDistance))]
	} else closeDist<-Inf
	
	closeDist
}


for(i in 1:dim(polyploidy)[1]) {
	pl<-as.character(polyploidy[i,"left_tip"])
	pr<-as.character(polyploidy[i,"right_tip"])

	if(pl==pr) { # this is a tip
		ww<-which(phy$tip.label==pl)
		polyNode<-ww
	} else { # this is an internal node
		pair<-c(pl, pr)
		polyNode<-getMRCA(phy, pair)
	}
	
	cs<-findClosestShift(polyNode, shiftNode)
	if(cs < Inf) {
		targetNode<-as.numeric(names(cs))
		thisResult<-getBranchPairDistance(phy, polyNode, targetNode)
		thisNodeDis<-lookupDistance[polyNode, targetNode]
		
		allDistances[i,]<-c(targetNode, thisResult,thisNodeDis)
	}
	else {
		allDistances[i,]<-c(NA, Inf, Inf, Inf, NA)
	}
}

timeResult<-cbind(polyploidy, allDistances)
write.csv(timeResult, file="timeresult.csv")