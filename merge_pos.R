merge_pos<-function(file_dir,outputfile){

 files=system(paste(paste("ls ",file_dir,sep=""),"*",sep=""),intern=T)		
 nf=length(files)
 res=NULL
 for (i in 1:nf){
       a=read.table(files[i],as.is=T)			
       a=a[which(rowSums(a[,3:7])>0),]
       if(nrow(a)>0){
	       res=rbind(res,a)
	}
 }
 write.table(res,outputfile,quote=F,row.name=F,col.name=F)
}
