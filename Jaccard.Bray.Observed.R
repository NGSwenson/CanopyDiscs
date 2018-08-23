##NG Swenson (swenson@umd.edu)
##22 August 2018

setwd("/Users/nateswenson/data/ctfs/jenny.overlap/output")
library(picante)

allom = read.csv("serc_allometry.csv")

sec = read.csv("serc.second.mat.csv",row.names=1) #rows individuals cols species
sec[is.na(sec)] = 0


tmp = as.data.frame(allom[,"Mnemonic"])
tmp2 = cbind(rownames(allom),tmp)
colnames(tmp2) = c("num","sp") #making a dataframe with first col is rowname indicator and second is species name
spp = split(tmp2,tmp2$sp) #splitting by species

beta1.out = matrix(NA, ncol = 4, nrow=length(spp))
beta2.out = matrix(NA, ncol = 4, nrow=length(spp))


for(i in 1:length(spp)){
	print(i)
	ind.mat = sec[spp[[i]][,1]  , -i] # get data matrix for all individuals of species i excluding own species
	jac.out = vegdist(ind.mat,method="jaccard") # calculate jaccard 
	bray.out = vegdist(ind.mat,method="bray") # calculate bray

	beta1.out[i,1] = nrow(ind.mat) #number of individuals of species i
	beta1.out[i,2] = mean(jac.out,na.rm=T) # mean jaccard between heterospecific neighborhoods
	beta1.out[i,3] = median(jac.out,na.rm=T) # median jaccard between heterospecific neighborhoods
	beta1.out[i,4] = sd(jac.out,na.rm=T) # st dev jaccard between heterospecific neighborhoods
	
	beta2.out[i,1] = nrow(ind.mat) #number of individuals of species i
	beta2.out[i,2] = mean(bray.out,na.rm=T) # mean bray between heterospecific neighborhoods
	beta2.out[i,3] = median(bray.out,na.rm=T) # median bray between heterospecific neighborhoods
	beta2.out[i,4] = sd(bray.out,na.rm=T) # st dev bray between heterospecific neighborhoods
}

rownames(beta1.out) = names(spp)
colnames(beta1.out) = c("num.indiv","mean.Jac","median.Jac","sd.Jac")

rownames(beta2.out) = names(spp)
colnames(beta2.out) = c("num.indiv","mean.Bray","median.Bray","sd.Bray")

write.table(beta1.out,"serc.jac.hetero.csv",sep=",",row.names=T,col.names=T,quote=F)
write.table(beta2.out,"serc.bray.hetero.csv",sep=",",row.names=T,col.names=T,quote=F)

