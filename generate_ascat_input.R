generate_ascat_input<-function(tumor_count,normal_count,outputfile,cov_cut){
	tn=read.table(tumor_count,as.is=T)
	nn=read.table(normal_count,as.is=T)
	
	tn=tn[which(tn[,7]>=cov_cut),]
	nn=nn[which(nn[,7]>=cov_cut),]
	
	tmp1=paste(tn[[1]],tn[[2]],sep="_")
	tmp2=paste(nn[[1]],nn[[2]],sep="_")
	
	indt=match(intersect(tmp1,tmp2),tmp1)
	indn=match(intersect(tmp1,tmp2),tmp2)
	
	tn=tn[indt,]
	nn=nn[indn,]
	
	tlogr=tn[,7]/nn[,7]
	tlogr=log2(tlogr/mean(tlogr))

	tbaf=apply(tn[,3:6],1,max)/tn[,7]
	nbaf=apply(nn[,3:6],1,max)/nn[,7]
	tbaf[which(tbaf==1)]=0
	nbaf[which(nbaf==1)]=0
	
	res=cbind(tn[,1:2],tlogr,matrix(0,nrow=nrow(tn),ncol=1),tbaf,nbaf)
	colnames(res)=c("chrs","pos","tlogr","nlogr","tbaf","nbaf")
	res[,1]=sapply(res[,1],function(z) strsplit(z,"chr")[[1]][2])
    	write.table(res,outputfile,sep="\t",row.names=F,col.names=T,quote=F)	
}
