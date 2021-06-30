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

counts=read.table("/home/user_25/cancer/rawData/TCGA/ESCA/counts.txt",header = T,row.names = 1)
labels=read.table("/home/user_25/cancer/rawData/TCGA/ESCA/ESCA_sample_sheet.2019-06-09.tsv",header = T,sep="\t")
counts=counts[grep("ENSG",rownames(counts)),]
name=colnames(counts)
labels$File.Name=gsub("-",".",labels$File.Name)
name=gsub("X","",name,ignore.case = FALSE)
labels=labels[labels$File.Name %in% name,]
labels$Sample.Type=gsub("Solid Tissue Normal","NC",labels$Sample.Type)
labels$Sample.Type=gsub("Primary Tumor","ESCA",labels$Sample.Type)
labels=cbind(labels$File.Name,labels$Sample.Type)
colnames(labels)=c("sample","label")
write.table(counts,"/home/user_25/cancer/rawData/TCGA/ESCA/count.txt",quote=F)
write.table(labels,"/home/user_25/cancer/rawData/TCGA/ESCA/label.txt",quote=F,row.names=F)

seq=c("ESCA","STAD","LUAD","COAD","READ","LIHC")
all_count=read.table("/home/user_25/cancer/rawData/TCGA/ESCA/counts.txt",header = T,row.names = 1)
all_label=read.table("/home/user_25/cancer/rawData/TCGA/ESCA/label.txt",header = T,sep="\t")
for (i in 2:length(seq)){
    c_path=paste("/home/user_25/cancer/rawData/TCGA/",seq[i],"/count.txt",sep="")
    l_path=paste("/home/user_25/cancer/rawData/TCGA/",seq[i],"/label.txt",sep="")
    c=read.table(c_path,header = T,row.names = 1)
    l=read.table(l_path,header = T,sep="\t")
    all_count=cbind(all_count,c)
    all_label=rbind(all_label,l)
}
