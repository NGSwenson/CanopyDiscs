#setwd("/home/nswenson/jenny.overlap")
#install.packages(c("spatstat","PBSmapping","sp"))
require(spatstat)
require(PBSmapping)
require(sp)
require(dplyr)
require(parallel)
library(foreach)
library(doParallel)

cores=detectCores()
cl <- makeCluster(cores[1]-1)
registerDoParallel(cl)

euc.dist <- function(data, point) {
  apply(data, 1, function (row) sqrt(sum((point - row) ^ 2)))
  }

#Function to calculate ellipse
ellipse <- function(a, b=a, xc=0, yc=0, theta=0, n=20){
  # a - length of long axis
  # b - length of short axis
  # xc,yc - center point of ellipse
  # theta - angle of long axis (in radians)
  # number of points on "polygon" of ellipse
  
  t <- seq(-pi,pi,length=n)
  x <- xc + a*cos(t)*cos(theta)/2 - b*sin(t)*sin(theta)/2
  y <- yc + a*cos(t)*sin(theta)/2 + b*sin(t)*cos(theta)/2
  return(data.frame(x = x, y = y))
  }

getOverlap <- function(e1, e2){
  n1 <- nrow(e1)
  n2 <- nrow(e2)
  names(e1) <- c("X","Y")
  names(e2) <- c("X","Y")
  E1 <- data.frame(PID = 1, POS = 1:n1, e1)
  E2 <- data.frame(PID = 2, POS = 1:n2, e2)
  E.both <- joinPolys(E1, E2)
  
  if(length(E.both) < 1){  
    return(0)
    }
  else{
    with(E.both, Polygon(cbind(X,Y)))@area
    }
}

myData = read.csv("bci_allometry.csv")
max.rad = max(myData$canopy_radius)

noms = rownames(myData)

foreach(i=1:nrow(myData),.packages=c("spatstat","PBSmapping","sp","dplyr")) %dopar% {

  print(date())
  tree.X = myData[i, "x"]
  tree.Y = myData[i, "y"]
  first.sub = subset(myData,  myData$x>=(tree.X-(max.rad+myData[i,"canopy_radius"])) & myData$x <=(tree.X+max.rad+myData[i,"canopy_radius"]) )
  second.sub = subset(first.sub,  first.sub$y>=(tree.Y-(max.rad+myData[i,"canopy_radius"])) & first.sub$y <=(tree.Y+max.rad+myData[i,"canopy_radius"]) )
  diss = euc.dist(second.sub[,6:7], second.sub[as.character(noms[i]),6:7])
  max.rad2 = max(second.sub$canopy_radius)
  neighs = names(diss[diss<= (max.rad2+myData[i,"canopy_radius"])])
  smData = second.sub[neighs,]
  #setwd("/home/nswenson/jenny.overlap/output_bci")
  cc= apply(smData, 1, function(tt){tt = as.numeric(tt);ellipse(tt[11],xc=tt[6],yc=tt[7])})
  myOutput = lapply(cc , function(x) lapply(cc, function(y) getOverlap(x,y)))
  myOutput2 = do.call(rbind,lapply(myOutput,unlist))
  myOutput3 = as.matrix(myOutput2[as.character(noms[i]),myOutput2[as.character(noms[i]),]>0])
  colnames(myOutput3) = i
  write.table(myOutput3,paste(as.character(noms[i]),".tmp.csv",sep=""),sep=",",row.names=T,col.names=T,quote=F)
  print(i)
  print(date())
  }
stopCluster(cl)
