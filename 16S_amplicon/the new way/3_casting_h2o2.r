library(reshape2) 																						#load package reshape for cast (not reshape2!)

#calculate sum of reads per sample in excel, then copy new table the relative abundances
#copy wholetax into second column, split one wholetax into taxonomic levels again (in excel)
#read into R
tcast_h2o2_rel<-read.csv(file.choose(),row.names=1,sep=",")

#wholetax also row names
tcast_h2o2_whole<-tcast_h2o2_rel[,c(1:(ncol(tcast_h2o2_rel)-6))] #takes all rows and the columns 1 to (max number of columns - 6) into new object
head(tcast_h2o2_whole)
#transpose back for merging
cast_h2o2_whole<-t(tcast_h2o2_whole)
#adjust names in excel for merging
write.csv(cast_h2o2_whole,file="cast_h2o2_whole.csv")
#create meta_data file in excelfor merging
#read in meta_data, sort stations in meta data
meta_h2o2<-read.csv(file.choose(),header=T,row.names=1,sep=";") 										#samples in rows, meta data in colums
str(meta_h2o2) 
#meta_h2o2$time<-as.factor(meta_h2o2$time)																#depending on the aims it might be necessary to change the classes of the "numeric" columns to "factor" columns
#meta_h2o2$days<-as.factor(meta_h2o2$days)
#meta_h2o2$parallel<-as.factor(meta_h2o2$parallel)
#str(meta_h2o2)

#merge:zwei tabellen, die eine gemeinsame spalte aufweisen (hier: row.names), werden zusammengesetzt. soviele zeilen, wie in x (der erstgenannten tabelle enthalten sind, kommen in die finale datei)
h2o2_wholetax_meta<-merge(cast_h2o2_whole,meta_h2o2,by="row.names",all.x=TRUE)			
#after merging the row.names need to be adjusted
head(h2o2_wholetax_meta)
row.names(h2o2_wholetax_meta)<-h2o2_wholetax_meta$Row.names
head(h2o2_wholetax_meta)
h2o2_wholetax_meta<-h2o2_wholetax_meta[,(-1)] #getting rid of the first column
head(h2o2_wholetax_meta)
#melt into final table with all samples on your chosen taxonomic level and all metadata combinations
final_h2o2_whole<-melt(h2o2_wholetax_meta, id=c("station","time","treatment","parallel","nucleic_acid","DV"))
head(final_h2o2_whole,50)
str(final_h2o2_whole)			#new melt function from reshape2 turns value into character format? no
#final_h2o2_whole$value<-as.numeric(final_h2o2_whole$value)
str(final_h2o2_whole)

#to aggregate means out of paralleles you can use:
final_h2o2_whole_mean<-aggregate(value~variable+time+treatment+station+DV+nucleic_acid, data = final_h2o2_whole, mean)					#possible step to combine parallels as mean or sum, all columns except the to-be-combined must be named

h2o2_otu_tax<-read.csv(file.choose(),row.names=1,sep=";")
final_h2o2_tax<-merge(final_h2o2_whole,h2o2_otu_tax,by.x="variable", by.y="row.names",all.x=TRUE)
final_h2o2_tax_mean<-merge(final_h2o2_whole_mean,h2o2_otu_tax,by.x="variable", by.y="row.names",all.x=TRUE)

#don't calculate from rel. abundance (counts smaller than 1)
head(cast_h2o2)

## for Shannon Index (H) you need Species richness (S) and Pielou's evenness (J):
S <- specnumber(table) 					## rowSums(table > 0) does the same...
J <- H/log(S)

##otus in spalten/proben als rownames 
h_cast_h2o2<-diversity(cast_h2o2)
s_cast_h2o2<-specnumber(cast_h2o2)
#falsch: j_cast_h2o2<-h_cast_h2o2/log(cast_h2o2)
write.csv(h_cast_h2o2,file="h_cast_h2o2.csv")
write.csv(s_cast_h2o2,file="s_cast_h2o2.csv")
write.csv(j_cast_h2o2,file="j_cast_h2o2.csv")
#j korrekt?
#über sample name mit meta daten kombinieren, in plots einbauen, z.B. mit gallaeci
#16s plot über z.b. pseudomonas bin plot?

##########################################################################################################################plot ideas
#you can further subset table or use everything for plotting. the plotting itself also provides subsetting options
sample_subset<-subset(final_h2o2_tax_mean,treatment=="h2o2")
more_cell_counts<-read.csv(file.choose(),sep=";")
#specific otu
#gallaeci
gtest<-subset(final_table_tax,habitat == "water" & grepl("Gallaeci",variable))
meta_h2o2_test<-subset(meta_h2o2,habitat =="water")
gg<-ggplot(gtest,aes(x=days))+
	geom_point(aes(y=value,colour=nucleic_acid))+
	geom_point(data=meta_h2o2_test,aes(y=shannon*5,colour=nucleic_acid,shape=nucleic_acid))+
	facet_wrap(~treatment,nrow=2)+
	ggtitle("gtest")+
	stat_summary(data=gtest,fun.y="mean",geom="line",aes(y=value,colour=nucleic_acid))+
	

##only plot diversity
divtest<-ggplot(meta_h2o2,aes(x=days))+
	geom_point(data=meta_h2o2,aes(y=meta_h2o2$shannon*100,colour=nucleic_acid))+
	geom_point(data=meta_h2o2,aes(y=meta_h2o2$richness,shape=nucleic_acid))+
	facet_wrap(~treatment*habitat,nrow=2,ncol=2)
divtest	
ggplot(df2, aes(x=x), y=y)) + stat_summary(fun.y="mean", geom="line", aes(group=factor(grouping)))	
#phnM_iso<-ggplot(molten_phnM_gene_isoforms, aes(x=day))+
	#geom_line(data=more_cell_counts, aes(y=glyph_mg_L*(cellfactor/15),colour="glyphosate concentration"),alpha=0.8,linetype="solid", size=1)+
	#geom_line(data=more_cell_counts, aes(y=glyph_theor*(cellfactor/15),colour="glyphosate dilution"),alpha=0.5,linetype="F1", size=1)+
	
##hier etwas ausgefeilter: ich will nur einen organismus darstellen, in einem punkt diagramm (geom_point). auf der x achse könnte ich zB die parallelen darstellen, auf der y-achse natürlich die read counts
#mit ggtitle gebe ich dem plot einen namen
#facet_wrap ermöglicht mir viele verschiedene plots gleichzeitig, die sich in den faktoren "day" und "treatment" unterschieden. diese vielen plots möchte ich in zwei reihen und 4 spalten dargestellt haben (wenn das ausreicht)
#die "aes"-eingabe ist ziemlich kompliziert, siehe ggplot2-manual, ich glaube, hier geht es um die legende

#gemeinschaftsplot

test<-ggplot(sample_subset, aes(x ="", y = value, fill=family))+		
	facet_wrap( ~station*time,nrow=2)+				
	geom_bar(width = 2, stat = "identity")+
	theme(legend.position='none')



#gemeinschaftsplot mit mindestreadsanzahl 
test_groesser_0.2<-ggplot(final_h2o2_tax_mean[which(final_h2o2_tax_mean$value>0.2),], aes(x = "", y = value, fill=variable))+
	facet_wrap( ~days*treatment,nrow=2,ncol=7)+				
	geom_bar(width = 1, stat = "identity")+
	theme(legend.position='none')



#ohne legende: +theme(legend.position='none')	

#um aus säulendiagramm ein tortendiagramm zu machen:
torten_test<-test + coord_polar("y", start=0)
