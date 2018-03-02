otu_table$wholetax<-do.call(paste, c(otu_table[c("genus","family","order","class","phylum")],sep="_")) 	#merge taxonomic classes to one name in a column called "wholetax"


library(reshape2) 																						#load package reshape for cast (not reshape2!)
#renaming
otu_table<-read.csv(file.choose(),row.names=1,sep=",")
otu_table<-otu_table[,c(2,8,1)]																			#select columns and reorder
head(otu_table)
cast_otu_table<-dcast(otu_table,sample_name~wholetax)
cast_otu_table[is.na(cast_otu_table)]<-0 																#replace NA with 0
names(cast_otu_table)
row.names(cast_otu_table)<-cast_otu_table$sample_name													#turns "sample_name" into row names
cast_otu_table<-cast_otu_table[,-1]																		#remove "samples_names"
nrow(cast_otu_table)																					#number of rows
ncol(cast_otu_table)																					#number of columns
otu_library_sizes<-rowSums(cast_start_otu) 																				#shows reads per sample
write.csv(otu_library_sizes,file="otu_library_sizes.csv")

names(cast_otu_table)
cast_otu_table<-as.data.frame(cast_otu_table)
write.csv(cast_otu_table,file="cast_otu_table.csv")
tcast_otu_table<-t(cast_otu_table)
write.csv(tcast_otu_table,file="tcast_otu_table.csv")													#ones to zeros in excel, then deleting otus that are not appearing any more (create (otu sum column/row and remove all zero lines/rows)
#calculate sum of reads per sample in excel, then copy new table the relative abundances
#copy wholetax into second column, split one wholetax into taxonomic levels again (in excel)
#read into R
tcast_otu_table2<-read.csv(file.choose(),row.names=1,sep=",")

#choose phylogenetic level to work upon, here wholetax
tcast_otu_table2_whole<-tcast_otu_table2[,c(6:ncol(tcast_otu_table2))] #takes all rows and the columns 6 to (max number of columns) into new object
head(tcast_otu_table2_whole)
#transpose back for merging
cast_otu_table3_whole<-t(tcast_otu_table2_whole)
#adjust names in excel for merging
write.csv(cast_otu_table3_whole,file="cast_otu_table3_whole.csv")
#create meta_data file in excelfor merging
#read in meta_data
meta_data<-read.csv(file.choose(),header=T,row.names=1,sep=";") 										#samples in rows, meta data in colums
str(meta_data) 
#meta_data$time<-as.factor(meta_data$time)																#depending on the aims it might be necessary to change the classes of the "numeric" columns to "factor" columns
#meta_data$days<-as.factor(meta_data$days)
#meta_data$parallel<-as.factor(meta_data$parallel)
#str(meta_data)

#merge:zwei tabellen, die eine gemeinsame spalte aufweisen (hier: row.names), werden zusammengesetzt. soviele zeilen, wie in x (der erstgenannten tabelle enthalten sind, kommen in die finale datei)
meta_counts_whole<-merge(cast_cut_all,meta_data,by="row.names",all.x=TRUE)			
#after merging the row.names need to be adjusted
head(meta_counts_whole)
row.names(meta_counts_whole)<-meta_counts_whole$Row.names
head(meta_counts_whole)
meta_counts_whole<-meta_counts_whole[,(-1)] #getting rid of the first column
head(meta_counts_whole)


#melt into final table with all samples on your chosen taxonomic level and all metadata combinations
final_table_whole<-melt(meta_counts_whole, id=c("time","days","treatment","parallel","nucleic_acid","habitat"))
head(final_table_whole,50)
str(final_table_whole)			#new melt function from reshape2 turns value into character format?
final_table_whole$value<-as.numeric(final_table_whole$value)
str(final_table_whole)

#to aggregate means out of paralleles you can use:
final_table_whole_mean<-aggregate(value~variable+time+days+treatment+habitat+nucleic_acid, data = final_table_whole, mean)					#possible step to combine parallels as mean or sum, all columns except the to-be-combined must be named

otu_tax<-read.csv(file.choose(),row.names=1,sep=";")
final_table_tax_mean<-merge(final_table_whole_mean,otu_tax,by.x="variable", by.y="row.names",all.x=TRUE)

#don't calculate from rel. abundance (counts smaller than 1)
head(cast_cut_all)

## for Shannon Index (H) you need Species richness (S) and Pielou's evenness (J):
S <- specnumber(table) 					## rowSums(table > 0) does the same...
J <- H/log(S)

##otus in spalten/proben als rownames 
h_cast_cut_all<-diversity(cast_cut_all)
s_cast_cut_all<-specnumber(cast_cut_all)
#falsch: j_cast_cut_all<-h_cast_cut_all/log(cast_cut_all)
write.csv(h_cast_cut_all,file="h_cast_cut_all.csv")
write.csv(s_cast_cut_all,file="s_cast_cut_all.csv")
write.csv(j_cast_cut_all,file="j_cast_cut_all.csv")
#j korrekt?
#über sample name mit meta daten kombinieren, in plots einbauen, z.B. mit gallaeci
#16s plot über z.b. pseudomonas bin plot?

##########################################################################################################################plot ideas
#you can further subset table or use everything for plotting. the plotting itself also provides subsetting options
sample_subset<-subset(final_table_tax_mean,habitat=="water"&nucleic_acid=="dna")
more_cell_counts<-read.csv(file.choose(),sep=";")
#specific otu
#gallaeci
gtest<-subset(final_table_tax,habitat == "water" & grepl("Gallaeci",variable))
meta_data_test<-subset(meta_data,habitat =="water")
gg<-ggplot(gtest,aes(x=days))+
	geom_point(aes(y=value,colour=nucleic_acid))+
	geom_point(data=meta_data_test,aes(y=shannon*5,colour=nucleic_acid,shape=nucleic_acid))+
	geom_point(data=more_cell_counts,aes(x=day, y=cells_ml/5000000,colour="cell counts per mL"))+
	facet_wrap(~treatment,nrow=2)+
	ggtitle("gtest")+
	stat_summary(data=gtest,fun.y="mean",geom="line",aes(y=value,colour=nucleic_acid))+
	stat_summary(data=meta_data_test, fun.y="mean",geom="line",aes(y=shannon*5,colour=nucleic_acid))

##only plot diversity
divtest<-ggplot(meta_data,aes(x=days))+
	geom_point(data=meta_data,aes(y=meta_data$shannon*100,colour=nucleic_acid))+
	geom_point(data=meta_data,aes(y=meta_data$richness,shape=nucleic_acid))+
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
	facet_wrap( ~treatment*days,nrow=2,ncol=18)+				
	geom_bar(width = 1, stat = "identity")+
	theme(legend.position='none')



#gemeinschaftsplot mit mindestreadsanzahl 
test_groesser_0.2<-ggplot(final_table_tax_mean[which(final_table_tax_mean$value>0.2),], aes(x = "", y = value, fill=variable))+
	facet_wrap( ~days*treatment,nrow=2,ncol=2)+				
	geom_bar(width = 1, stat = "identity")+
	theme(legend.position='none')



#ohne legende: +theme(legend.position='none')	

#um aus säulendiagramm ein tortendiagramm zu machen:
torten_test<-test + coord_polar("y", start=0)
