##NG Swenson (swenson@umd.edu)
##22 August 2018

setwd("/Users/nateswenson/data/ctfs/jenny.overlap/output")
library(picante)

allom = read.csv("serc_allometry.csv")

sec = read.csv("serc.second.mat.csv",row.names=1)
sec[is.na(sec)] = 0


tmp = as.data.frame(allom[,"Mnemonic"])
tmp2 = cbind(rownames(allom),tmp)
colnames(tmp2) = c("num","sp")

# a function that calculates jaccard distance between heterospecific neighborhoods
# lapply'd to list
# performs randomization of species names (i.e. completely randomized forest)
rand.jac = function(xx){

	xx$num = sample(xx$num)
	spp = split(xx,xx$sp)

	#generates species level community data matrix (CDMs)
	com.list = list()
		for(i in 1:length(spp)){
		nam = names(spp)[i]
		ind.mat = sec[spp[[i]][,1]  , -i]
		com.list[[nam]] = ind.mat
	}
	
	jac.fun = function(x){
		jacs = vegdist(x,method="jaccard")
		mean.jacs = mean(jacs,na.rm=T)
		median.jacs = median(jacs,na.rm=T)
		sd.jacs = sd(jacs,na.rm=T)
		output = cbind(mean.jacs,median.jacs,sd.jacs)
		return(output)
	}
	return(lapply(com.list, jac.fun)) #applying across list of species neighborhood CDMs
	
	}

#replicating the rand.jac function 999 times each time with a random forest
rand.output = replicate(999,rand.jac(tmp2))

means = matrix(NA, nrow=nrow(rand.output),ncol=999)
medians = matrix(NA, nrow=nrow(rand.output),ncol=999)
sds = matrix(NA, nrow=nrow(rand.output),ncol=999)

#summarizing output
for(i in 1:nrow(rand.output)){
	for(j in 1:999){
	means[i,j] =	rand.output[i,][[j]][1,1]
	medians[i,j] =	rand.output[i,][[j]][1,2]
	sds[i,j] =	rand.output[i,][[j]][1,3]	
	}
}



write.table(means,"serc.jac.rand.means.csv",sep=",",row.names=F,col.names=T,quote=F)
write.table(medians,"serc.jac.rand.medians.csv",sep=",",row.names=F,col.names=T,quote=F)
write.table(sds,"serc.jac.rand.sds.csv",sep=",",row.names=F,col.names=T,quote=F)


