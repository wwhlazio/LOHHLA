run_ascat<-function(a,dir,out,sample){

setwd(dir)
tlogr=a[,c(1:3)]
row.names(tlogr)=paste("snp",c(1:nrow(tlogr)),sep='')
colnames(tlogr)[3]="s1"
nlogr=a[,c(1:2,4)]
row.names(nlogr)=paste("snp",c(1:nrow(nlogr)),sep='')
colnames(nlogr)[3]="s1"
tbaf=a[,c(1:2,5)]
row.names(tbaf)=paste("snp",c(1:nrow(tbaf)),sep='')
colnames(tbaf)[3]="s1"
nbaf=a[,c(1:2,6)]
row.names(nbaf)=paste("snp",c(1:nrow(nbaf)),sep='')
colnames(nbaf)[3]="s1"
write.table(tlogr,'Tumor_logR.txt',row.names=T,col.names=T,quote=F,sep="\t")
write.table(nlogr,'Normal_logR.txt',row.names=T,col.names=T,quote=F,sep="\t")
write.table(tbaf,'Tumor_baf.txt',row.names=T,col.names=T,quote=F,sep="\t")
write.table(nbaf,'Normal_baf.txt',row.names=T,col.names=T,quote=F,sep="\t")
library(ASCAT)
ascat.bc = ascat.loadData("Tumor_logR.txt","Tumor_baf.txt","Normal_logR.txt","Normal_baf.txt")
ascat.bc = ascat.aspcf(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc,gamma=1)
res=matrix(c(ascat.output$ploidy,ascat.output$aberrantcellfraction),nrow=1)
colnames(res)=c("tumorPloidy","tumorPurity")
row.names(res)=paste(sample,"_tumor",sep="")
write.table(res,out,col.names=T,row.names=T,sep="\t",quote=F)
}
