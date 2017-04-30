library("DESeq2")
args <- commandArgs(trailingOnly = TRUE)
count_table		<- args[1]
conditions		<- args[2]
out			<- args[3]

countData		<- as.matrix(read.table(count_table,header=TRUE, sep = "\t",row.names=1,as.is=TRUE))
colData 		<- read.csv(conditions,sep = "\t",row.names=1)

dds 			<- DESeqDataSetFromMatrix(countData=countData, colData= colData,design=~condition)
dds 			<- DESeq(dds)
res 			<- results(dds)

write.table(res,file=out,sep=",")
