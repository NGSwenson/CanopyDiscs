##NG Swenson (swenson@umd.edu)
##21 August 2018

#setwd("/Users/ngswenson/tmp")
myData = read.csv("amacayacu_allometry.csv")
#setwd("/Users/ngswenson/tmp/output_amacayacu")

listy = list.files()

first.mat = matrix(NA, nrow=length(listy),ncol=4)

for(i in 1:length(listy)){
  print(i)
  tmp = read.csv(paste(i,".tmp.csv",sep=""),header=T,row.names=1)
  colnames(tmp) = gsub("X","",colnames(tmp))
  first.mat[i,1] = as.matrix(myData[colnames(tmp),"dbh"]) #dbh
  first.mat[i,2] = nrow(tmp)-1 #number of overlaps
  first.mat[i,3] = sum(myData[rownames(tmp),"Species"]==myData[colnames(tmp),"Species"])-1 # num conspecifics
  first.mat[i,4] = first.mat[i,2] - (first.mat[i,3]+1) # num heterospecifics
  }

first.mat = cbind(as.data.frame(myData$Species), first.mat)


second.mat = matrix(NA, nrow=length(listy),ncol=length(unique(myData$Species)))
colnames(second.mat) = unique(myData$Species)
rownames(second.mat)= gsub(".tmp.csv","",listy)

for(i in 1:length(listy)){
  print(i)
  tmp = read.csv(paste(i,".tmp.csv",sep=""),header=T,row.names=1)
  colnames(tmp) = gsub("X","",colnames(tmp))
  if(nrow(tmp)==1){
    second.mat[i,] = NA
    }
  else{
    tmp2 = subset(tmp,rownames(tmp)!=colnames(tmp))
    tmp3= cbind(tmp2, myData[rownames(tmp2),"Species"])
    second.mat[i,] = tapply(tmp3[,1],tmp3[,2],sum)/myData[colnames(tmp),"canopy_radius"]
    }
}

