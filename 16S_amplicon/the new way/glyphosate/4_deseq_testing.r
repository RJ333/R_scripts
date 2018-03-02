## deseq new approach
##tcast otu tables selected for glyph data times 4 to 10 (day 44 to 83), 3 before and 4 after addition
##reinstall DESeq2 (htmltools reinstalled)
library(DESeq2)


##################################################################################dna 5 bis 11

coldata_dna5_11<-read.csv(file.choose(),row.names=1,sep=";") #meta_data, samples in rows, meta data in columns, same order as count data
countdata_dna5_11<-read.csv(file.choose(),row.names=1,sep=";")
names(countdata_dna5_11)
countdata_dna5_11<-countdata_dna5_11[,(-22)]
names(countdata_dna5_11) #removing sums column
rownames(coldata_dna5_11)
all(rownames(coldata_dna5_11) %in% colnames(countdata_dna5_11)) #check übereinstimmung
all(rownames(coldata_dna5_11) == colnames(countdata_dna5_11))

dds_dna5_11 <- DESeqDataSetFromMatrix(countData = countdata_dna5_11,
                              colData = coldata_dna5_11,
                              design = ~ condition)
dds_dna5_11
dds_dna5_11$condition <- relevel(dds_dna5_11$condition, ref = "untreated")

dds_dna5_11<-DESeq(dds_dna5_11)
res_dna5_11<-results(dds_dna5_11,alpha=0.05)
write.csv(res_dna5_11,file="res_dna5_11.csv")
summary(res_dna5_11)
sum(res_dna5_11$padj < 0.1, na.rm=TRUE) #26
sum(res_dna5_11$padj < 0.05, na.rm=TRUE) #24
sum(res_dna5_11$padj < 0.01, na.rm=TRUE) #20
sum(res_dna5_11$padj < 0.001, na.rm=TRUE) #15


################################################################cdna 5 bis 5_11

coldata_cdna5_11<-read.csv(file.choose(),row.names=1,sep=";") #meta_data, samples in rows, meta data in columns, same order as count data
countdata_cdna5_11<-read.csv(file.choose(),row.names=1,sep=";")
names(countdata_cdna5_11)
countdata_cdna5_11<-countdata_cdna5_11[,(-22)]
names(countdata_cdna5_11) #removing sums column
rownames(coldata_cdna5_11)
all(rownames(coldata_cdna5_11) %in% colnames(countdata_cdna5_11)) #check übereinstimmung
all(rownames(coldata_cdna5_11) == colnames(countdata_cdna5_11))

dds_cdna5_11 <- DESeqDataSetFromMatrix(countData = countdata_cdna5_11,
                              colData = coldata_cdna5_11,
                              design = ~ condition)
dds_cdna5_11
dds_cdna5_11$condition <- relevel(dds_cdna5_11$condition, ref = "untreated")

dds_cdna5_11<-DESeq(dds_cdna5_11)
res_cdna5_11<-results(dds_cdna5_11,alpha=0.05)
write.csv(res_cdna5_11,file="res_cdna5_11.csv")
summary(res_cdna5_11)
sum(res_cdna5_11$padj < 0.1, na.rm=TRUE) #41
sum(res_cdna5_11$padj < 0.05, na.rm=TRUE) #38
sum(res_cdna5_11$padj < 0.01, na.rm=TRUE) #26
sum(res_cdna5_11$padj < 0.001, na.rm=TRUE) #19