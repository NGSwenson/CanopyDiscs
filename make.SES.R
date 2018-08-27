##NG Swenson (swenson@umd.edu)
##26 August 2018

#setwd("/Users/ngswenson/data/jenny.overlap/output")

obs = read.csv("palanan.jac.hetero.csv",header=T,row.names=1)

jac.means = read.csv("palanan.jac.mean.nulls.csv",header=T,row.names=1)

	null.mean = apply(jac.means,1,mean,na.rm=T)
	null.sd = apply(jac.means,1,sd,na.rm=T)
	ses = as.matrix((obs[,2]-null.mean)/null.sd)
	
	rownames(ses) = rownames(obs)
	colnames(ses)= "SES"
	p.value = t(apply(cbind(obs[,2],jac.means),1,rank))[,1]
	output = cbind(ses,p.value)


write.table(output,"palanan.jac.means.ses.csv",sep=",",col.names=T,row.names=T,quote=F)



jac.medians = read.csv("palanan.jac.median.nulls.csv",header=T,row.names=1)

	null.mean = apply(jac.medians,1,mean,na.rm=T)
	null.sd = apply(jac.medians,1,sd,na.rm=T)
	ses = as.matrix((obs[,3]-null.mean)/null.sd)
	
	rownames(ses) = rownames(obs)
	colnames(ses)= "SES"
	p.value = t(apply(cbind(obs[,3],jac.means),1,rank))[,1]
	output = cbind(ses,p.value)


write.table(output,"palanan.jac.medians.ses.csv",sep=",",col.names=T,row.names=T,quote=F)

