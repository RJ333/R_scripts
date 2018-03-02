#tnt-daten-schlosser
#read relative abundances into R
tnt_cast_rel_tax<-read.csv(file.choose(),row.names=1,sep=";")
tnt_cast_rel_tax[is.na(tnt_cast_rel_tax)] <- 0
ttnt<-t(tnt_cast_rel_tax)
meta_tnt<-read.csv(file.choose(),header=T,row.names=1,sep=";") 

#create meta_data file in excelfor merging
#read in meta_data
meta_tnt<-read.csv(file.choose(),header=T,row.names=1,sep=";") 										#samples in rows, meta data in colums
str(meta_tnt)

#merge:zwei tabellen, die eine gemeinsame spalte aufweisen (hier: row.names), werden zusammengesetzt. soviele zeilen, wie in x (der erstgenannten tabelle enthalten sind, kommen in die finale datei)
tnt_rel_abu_with_meta<-merge(ttnt,meta_tnt,by="row.names",all.x=TRUE)			
#after merging the row.names need to be adjusted
head(tnt_rel_abu_with_meta)
row.names(tnt_rel_abu_with_meta)<-tnt_rel_abu_with_meta$Row.names
head(tnt_rel_abu_with_meta)
tnt_rel_abu_with_meta<-tnt_rel_abu_with_meta[,(-1)] #getting rid of the first column
names(tnt_rel_abu_with_meta)
head(tnt_rel_abu_with_meta)
write.csv(tnt_rel_abu_with_meta,file="tnt_rel_abu_with_meta.csv")

#melt into final table with all samples on your chosen taxonomic level and all metadata combinations
library(reshape2)
tnt_rel_abu_molten<-melt(tnt_rel_abu_with_meta, id=c("habitat","treatment"))
head(tnt_rel_abu_molten,50)
tail(tnt_rel_abu_molten,50) #check if the end of list is correct, otherwise you probably missed a column name in the melting command
str(tnt_rel_abu_molten)			#new melt function from reshape2 turns value into character format? not always
#tnt_rel_abu_molten$value<-as.numeric(tnt_rel_abu_molten$value)
#str(tnt_rel_abu_molten)

#add taxonomic levels
tnt_tax<-read.csv(file.choose(),row.names=1,sep=";")
tnt_rel_abu_molten_tax<-merge(tnt_rel_abu_molten,tnt_tax,by.x="variable", by.y="row.names",all.x=TRUE)
str(tnt_rel_abu_molten_tax)
write.csv(tnt_rel_abu_molten_tax,file="tnt_rel_abu_molten_tax.csv")

#gemeinschaftsplot  
test<-ggplot(tnt_rel_abu_molten_tax[which(tnt_rel_abu_molten_tax$value>0.005),],aes(x ="", y = value, fill=genus))+
	facet_wrap( ~habitat*treatment,nrow=1)+
	geom_bar(width = 1, stat = "identity")+
	geom_text(data=subset(tnt_rel_abu_molten_tax,value>0.005),aes(label=genus, vjust=1), position="stack",size=2.5)+
	theme(legend.position='none')
test

###diversity with vegan

##otus in spalten/proben als rownames 
h_ttnt<-diversity(ttnt)
s_ttnt<-specnumber(ttnt)
#j_ttnt<-h_ttnt/ln(ttnt) that in excel
write.csv(h_ttnt,file="h_ttnt.csv")
write.csv(s_ttnt,file="s_ttnt.csv")

#nmds
tnt_nmds_orig<-metaMDS(ttnt,try=100,autotransform=FALSE)
tnt_nmds_transformed<-metaMDS(ttnt,try=100,autotransform=TRUE)
plot(tnt_nmds_transformed,type="t",display="sites")
text(tnt_nmds_transformed,labels=row.names(tnt_nmds_transformed))
plot(tnt_nmds_orig,type="t",display="sites")
text(tnt_nmds_orig,labels=row.names(tnt_nmds_orig))
##########################################################################################################################plot ideas





