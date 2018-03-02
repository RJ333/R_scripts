#Deseq nutzt die count table und eine metadaten table als input. 
#wichtig: 
#metatable: samples in rows, meta data in colums; 
#count table: samples in columns, otus als rows



####################################################################
#wenn mich nur die transformierten counts interessieren, also von DESeq aufgrund meiner Proben schätzt, wie die tatsächliche Verteilung von OTUS/Genen aussieht:

#create meta_data file in excel
#read in meta_data
#integer to factor
meta_data<-read.csv(file.choose(),header=T,row.names=1,sep=";") 										#samples in rows, meta data in colums
str(meta_data) 
meta_data$time<-as.factor(meta_data$time)
meta_data$days<-as.factor(meta_data$days)
meta_data$parallel<-as.factor(meta_data$parallel)
str(meta_data)

library(DESeq2)
library(ggplot2)

#giving data to deseq: wichtige variable am ende bei "~x + y" design
dds_start <- DESeqDataSetFromMatrix(countData = tcast_start_otu,colData = meta_data,design = ~time)
#running a LRT likelihood ratio test for testing several factors at once, reduced: reducing the degrees of freedom by 1, local dispersion fit
dds_start<-DESeq(dds_start,test="LRT",reduced= ~ 1, fitType='local') 

sizeFactors(dds_start)
write.csv(sizeFactors(dds_start),file="factor_start.csv")
counts_start<-counts(dds_start, normalized=T) #das sind die transformierten counts
counts_start
write.csv(counts_start,file="norm_start.csv")

#########################################################################
#wenn ich richtige analysen durchführen will (hier sind die metadaten wichtig)

#transponieren (proben als columns, otus als rows)
tcast_raw_control_dna_biofilm<-t(cast_raw_control_dna_biofilm)
head(tcast_raw_control_dna_biofilm)

#metadateien einlesen und zahlen zu faktoren:
meta_control_dna_biofilm<-read.csv(file.choose(),header=T,row.names=1,sep=";") #proben als rows, infos in colums
str(meta_control_dna_biofilm)
meta_control_dna_biofilm$time<-as.factor(meta_control_dna_biofilm$time)
meta_control_dna_biofilm$days<-as.factor(meta_control_dna_biofilm$days)
meta_control_dna_biofilm$parallel<-as.factor(meta_control_dna_biofilm$parallel)
str(meta_control_dna_biofilm)

#relevel für time 6 als reference (vorsicht beim plotten, dafür kann "days" in originalreihenfolge genutzt werden): 
meta_control_dna_biofilm$time <- relevel(meta_control_dna_biofilm$time, "6")
str(meta_control_dna_biofilm$time) 

library(DESeq2)
library(pheatmap)
library(ggplot2)

#giving data to deseq: wichtige variable am ende bei "~x + y" design, hier zum Beispiel "time", der erste Eintrag wird als Referenz genutzt, daher oben der "relevel"-Schritt
dds_control_dna_biofilm <- DESeqDataSetFromMatrix(countData = tcast_raw_control_dna_biofilm,colData = meta_control_dna_biofilm,design = ~time)
#running a LRT likelihood ratio test for testing several factors at once, reduced: reducing the degrees of freedom by 1, local dispersion fit
dds_control_dna_biofilm<-DESeq(dds_control_dna_biofilm,test="LRT",reduced= ~ 1, fitType='local') 

#checking and saving the results from the overall test and from the several 1:1-comparisons
#mit "results" wird bei mehr als zwei levels für den vergleichenden faktor (hier time mit 16 levels) nur die letzte analyse angezeigt -> t19 gegen t6
resultsNames(dds_control_dna_biofilm)
res_condb<-results(dds_control_dna_biofilm) 
res_condb5<-results(dds_control_dna_biofilm,name="time_5_vs_6")
res_condb8<-results(dds_control_dna_biofilm,name="time_8_vs_6")
res_condb11<-results(dds_control_dna_biofilm,name="time_11_vs_6")
res_condb13<-results(dds_control_dna_biofilm,name="time_13_vs_6")
res_condb15<-results(dds_control_dna_biofilm,name="time_15_vs_6")
res_condb17<-results(dds_control_dna_biofilm,name="time_17_vs_6")
res_condb19<-results(dds_control_dna_biofilm,name="time_19_vs_6")

#sum kann immer den gleichen wert an signifikanten geben, aber evtl unterschiedlich anzahl hoch/runter --> summary
sum(res_condb19$padj < 0.1, na.rm=TRUE)
#55
sum(res_condb5$padj < 0.1, na.rm=TRUE)
#55

summary(res_condb5)
summary(res_condb8)
summary(res_condb11)
summary(res_condb13)
summary(res_condb15)
summary(res_condb17)
summary(res_condb19)


#We subset the results table and sort it by the log2 fold change estimate, then we save the ordered significants table and a complete table:

otuSig <-subset(res_condb, padj < 0.1)
otuSig <-otuSig[order(otuSig$log2FoldChange),]
otuAll <-res_condb
write.csv(as.data.frame(otuSig),file="con_db.csv")
write.csv(as.data.frame(otuAll),file="con_db_all.csv")

otuSig5 <-subset(res_condb5, padj < 0.1)
otuSig5 <-otuSig5[order(otuSig5$log2FoldChange),]
otuAll5 <-res_condb5
write.csv(as.data.frame(otuSig5),file="con_db5.csv")
write.csv(as.data.frame(otuAll5),file="con_db5_all.csv")

otuSig8 <-subset(res_condb8, padj < 0.1)
otuSig8 <-otuSig8[order(otuSig8$log2FoldChange),]
otuAll8 <-res_condb8
write.csv(as.data.frame(otuSig8),file="con_db8.csv")
write.csv(as.data.frame(otuAll8),file="con_db8_all.csv")

otuSig11 <-subset(res_condb11, padj < 0.1)
otuSig11 <-otuSig11[order(otuSig11$log2FoldChange),]
otuAll11 <-res_condb11
write.csv(as.data.frame(otuSig11),file="con_db11.csv")
write.csv(as.data.frame(otuAll11),file="con_db11_all.csv")

otuSig13 <-subset(res_condb13, padj < 0.1)
otuSig13 <-otuSig13[order(otuSig13$log2FoldChange),]
otuAll13 <-res_condb13
write.csv(as.data.frame(otuSig13),file="con_db13.csv")
write.csv(as.data.frame(otuAll13),file="con_db13_all.csv")

otuSig15 <-subset(res_condb15, padj < 0.1)
otuSig15 <-otuSig15[order(otuSig15$log2FoldChange),]
otuAll15 <-res_condb15
write.csv(as.data.frame(otuSig15),file="con_db15.csv")
write.csv(as.data.frame(otuAll15),file="con_db15_all.csv")

otuSig17 <-subset(res_condb17, padj < 0.1)
otuSig17 <-otuSig17[order(otuSig17$log2FoldChange),]
otuAll17 <-res_condb17
write.csv(as.data.frame(otuSig17),file="con_db17.csv")
write.csv(as.data.frame(otuAll17),file="con_db17_all.csv")

otuSig19 <-subset(res_condb19, padj < 0.1)
otuSig19 <-otuSig19[order(otuSig19$log2FoldChange),]
otuAll19 <-res_condb19
write.csv(as.data.frame(otuSig19),file="con_db19.csv")
write.csv(as.data.frame(otuAll19),file="con_db19_all.csv")