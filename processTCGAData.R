names=c("BRCA","COAD","ESCA","LIHC","LUAD","LUSC","READ","STAD")
for (seq in names[1:length(names)]){
count=data.frame()
path=paste("/home/user_25/cancer/rawData/TCGA/",seq,"/htseq_counts",sep="")
files=list.files(path)
count=read.table(gzfile(paste(path,'/',files[1],sep=""))) # read in first colume as dataframe
colnames(count)=c("gene",files[1])
for (i in files[2:length(files)]){
    temp=read.table(gzfile(paste(path,'/',i,sep=""))) # read in .gz data downloaded from TCGA
    colnames(temp)=c("gene",i)
    count=merge(count,temp,by='gene')
    print(ncol(count)/length(files))
}
write.table(count,paste("/home/user_25/cancer/rawData/TCGA/",seq,"/counts.txt",sep=""),sep='\t',quote=F,row.names=F)
}
