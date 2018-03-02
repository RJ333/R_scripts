otu_table$wholetax<-do.call(paste, c(otu_table[c("genus","family","order","class","phylum")],sep="_")) 	#merge taxonomic classes to one name in a column called "wholetax"


library(reshape2) 																						#load package reshape for cast (not reshape2!)
#renaming
otu_table<-read.csv(file.choose(),row.names=1,sep=",")
otu_table<-otu_table[,c(2,8,1)]																			#select columns and reorder
head(otu_table)
cast_otu_table<-dcast(otu_table,sample_name~wholetax)
cast_otu_table[is.na(cast_otu_table)]<-0 																#replace NA with 0
names(cast_otu_table)
head(cast_otu_table)
row.names(cast_otu_table)<-cast_otu_table$sample_name													#turns "sample_name" into row names
cast_otu_table<-cast_otu_table[,-1]																		#remove "samples_names"
names(cast_otu_table)
head(cast_otu_table)
nrow(cast_otu_table)																					#number of rows
ncol(cast_otu_table)																					#number of columns
otu_library_sizes<-rowSums(cast_otu_table) 																#shows reads per sample, but missing the removed singletons!
write.csv(otu_library_sizes,file="otu_library_sizes.csv")

names(cast_otu_table)
#cast_otu_table<-as.data.frame(cast_otu_table)
write.csv(cast_otu_table,file="cast_otu_table.csv")
tcast_otu_table<-t(cast_otu_table)
write.csv(tcast_otu_table,file="tcast_otu_table.csv")
#calculate sum of reads per sample in excel, then copy new table the relative abundances
#read relative abundances into R
tcast_otu_table2<-read.csv(file.choose(),row.names=1,sep=";")

#transpose back for merging
cast_otu_table3<-t(tcast_otu_table2)
head(cast_otu_table3)
#adjust names in excel for merging
write.csv(cast_otu_table3,file="cast_otu_table_rel.csv")
#create meta_data file in excelfor merging
#read in meta_data
meta_data<-read.csv(file.choose(),header=T,row.names=1,sep=";") 										#samples in rows, meta data in colums
str(meta_data)
levels(meta_data$cellcounts) <- c(levels(meta_data$cellcounts),"NA")
meta_data$cellcounts[is.na(meta_data$cellcounts)] <- "NA"
str(meta_data) #"NA" in cell counts are now strings, which is necessary for aggregate function later
#meta_data sollte keine richness beinhalten: sonst können die Parallelen nicht gemittelt werden


#meta_data$time<-as.factor(meta_data$time)																#depending on the aims it might be necessary to change the classes of the "numeric" columns to "factor" columns
#meta_data$days<-as.factor(meta_data$days)
#meta_data$parallel<-as.factor(meta_data$parallel)
#str(meta_data)

#merge:zwei tabellen, die eine gemeinsame spalte aufweisen (hier: row.names), werden zusammengesetzt. soviele zeilen, wie in x (der erstgenannten tabelle enthalten sind, kommen in die finale datei)
rel_abu_with_meta<-merge(cast_otu_table3,meta_data,by="row.names",all.x=TRUE)			
#after merging the row.names need to be adjusted
head(rel_abu_with_meta)
row.names(rel_abu_with_meta)<-rel_abu_with_meta$Row.names
head(rel_abu_with_meta)
rel_abu_with_meta<-rel_abu_with_meta[,(-1)] #getting rid of the first column
names(rel_abu_with_meta)
head(rel_abu_with_meta)
write.csv(rel_abu_with_meta,file="rel_abu_with_meta.csv")
rel_abu_with_meta_mean<-aggregate(c(1~2:659,661:665),data = rel_abu_molten, FUN=mean)

#melt into final table with all samples on your chosen taxonomic level and all metadata combinations
rel_abu_molten<-melt(rel_abu_with_meta, id=c("time","days","treatment","parallel","nucleic_acid","habitat","richness","cellcounts","disturbance"))
head(rel_abu_molten,50)
tail(rel_abu_molten,50) #check if the end of list is correct, otherwise you probably missed a column name in the melting command
str(rel_abu_molten)			#new melt function from reshape2 turns value into character format? not always
rel_abu_molten$value<-as.numeric(rel_abu_molten$value)
str(rel_abu_molten)

#to aggregate means out of paralleles you can use (possible step to combine parallels as mean or sum, all columns except the to-be-combined must be named):
rel_abu_molten_mean<-aggregate(value~time+days+treatment+nucleic_acid+habitat+variable+cellcounts+disturbance, data = rel_abu_molten, mean) 
write.csv(rel_abu_molten_mean,file="rel_abu_molten_mean.csv")
#add taxonomic levels
otu_tax<-read.csv(file.choose(),row.names=1,sep=";")
rel_abu_molten_tax<-merge(rel_abu_molten,otu_tax,by.x="variable", by.y="row.names",all.x=TRUE)
str(rel_abu_molten_tax)
rel_abu_molten_mean_tax<-merge(rel_abu_molten_mean,otu_tax,by.x="variable", by.y="row.names",all.x=TRUE)
write.csv(rel_abu_molten_tax,file="rel_abu_molten_tax.csv")
write.csv(rel_abu_molten_mean_tax,file="rel_abu_molten_mean_tax.csv")

#-------------------------------------------------------------------------------
#absolute values by using cell counts


testcell<-subset(rel_abu_with_meta,habitat=="water")
tmp <- merge(testcell,more_cell_counts, by=c("time","treatment"), all.x=TRUE, all.y=FALSE)
#row.names(tmp)<-row.names(testcell)#not in right order
ttestcell<-t(tmp)
ttestcell<-as.data.frame(ttestcell)
write.csv(ttestcell,file="ttestcell.csv")
#in excel relabu durch cell counts teilen, steigende otus suchen!extra durch 100 teilen (prozent)

t_abs_abu_water_otus<-read.csv(file.choose(),sep=";",row.names=1) #ttestcell_ready
abs_abu_water_otus<-t(t_abs_abu_water_otus)
abs_abu_water_otus<-as.data.frame(abs_abu_water_otus)
abs_abu_water_otus$time<-as.numeric(as.character(abs_abu_water_otus$time))
abs_abu_water_otus2<-subset(abs_abu_water_otus,time > 2)
abs_abu_molten<-melt(abs_abu_water_otus2, id=c("time","days","treatment","parallel","nucleic_acid","habitat","richness","cells_ml","cells_se","glyph_mg_L","glyph_se","glyph_theor"))
#warnmeldung because values are seen as factors, coerced to characters
head(abs_abu_molten,100)
str(abs_abu_molten)
abs_abu_molten$value<-as.numeric(abs_abu_molten$value)
abs_abu_molten$parallel<-as.integer(as.character(abs_abu_molten$parallel))
str(abs_abu_molten)
write.csv(abs_abu_molten,file="abs_abu_molten.csv")
#abs_abu_molten wholetax column split in excel to taxonomic levels, saved as abs_abu_molten_tax
abs_abu_molten_tax<-read.csv(file.choose(),row.names=1, sep=";")
#trying different things because mean (below) is not calculated, but after removal of columns it works!
abs_abu_molten2<-abs_abu_molten[,c(1:5,13,14)] 
abs_abu_molten_mean<-aggregate(value~time+days+treatment+nucleic_acid+variable, data = abs_abu_molten2, mean) 
write.csv(abs_abu_molten_mean,file="abs_abu_molten_mean.csv")
#abs_abu_molten_mean wholetax column split in excel to taxonomic levels, saved as abs_abu_molten_mean_tax
abs_abu_molten_mean_tax<-read.csv(file.choose(),row.names=1, sep=";")
#in excel concatenate to unique names except for treatment, sort by new name, then substract every second line
abs_abu_substract<-read.csv(file.choose(),sep=";")
abs_abu_substract_tax<-merge(abs_abu_substract,otu_tax,by.x="variable", by.y="row.names",all.x=TRUE)




#not sure if below is correct (mean calculated?)
#abs_abu_molten_tax<-merge(abs_abu_molten,otu_tax,by.x="variable", by.y="row.names",all.x=TRUE)
#abs_abu_molten_tax$days<-as.numeric(as.character(abs_abu_molten_tax$days))
#str(abs_abu_molten)
#abs_abu_molten_mean_tax<-merge(abs_abu_molten_mean,otu_tax,by.x="variable", by.y="row.names",all.x=TRUE)

write.csv(abs_abu_molten_mean_tax,file="abs_abu_molten_mean_tax.csv")
#-----------------------------------------------------------------------------------------
#substract control values from treatment: test:
df <- data.frame(id1=gl(2, 3, labels=c("a", "b")),
                 id2=rep(gl(3, 1, labels=c("live1", "live2", "killed")), 2), 
                 id3=gl(3, 2, labels=c("hello", "howdy","goodbye")),
				 y=c(10, 10, 1, 12, 12, 2),
                 otherFactor = gl(3, 2))
library(plyr)
df2 <- ddply(df, .(id1), transform, y = y-y[id2=="killed"&id3=="howdy"])
df2[-which(df2$id2=="killed"),]
----------------------------------------------------------------------------------------------------

#don't calculate from rel. abundance (counts smaller than 1)
head(cast_otu_table3)

## for Shannon Index (H) you need Species richness (S) and Pielou's evenness (J):
S <- specnumber(table) 					## rowSums(table > 0) does the same...
J <- H/log(S)

##otus in spalten/proben als rownames 
h_cast_otu_table3<-diversity(cast_otu_table3)
s_cast_otu_table3<-specnumber(cast_otu_table3)
#falsch: j_cast_otu_table3<-h_cast_otu_table3/log(cast_otu_table3)
write.csv(h_cast_otu_table3,file="h_cast_otu_table3.csv")
write.csv(s_cast_otu_table3,file="s_cast_otu_table3.csv")
write.csv(j_cast_otu_table3,file="j_cast_otu_table3.csv")
#j korrekt?
#über sample name mit meta daten kombinieren, in plots einbauen, z.B. mit gallaeci
#16s plot über z.b. pseudomonas bin plot?

##########################################################################################################################plot ideas
#you can further subset table or use everything for plotting. the plotting itself also provides subsetting options
sample_subset<-subset(abs_abu_molten_mean_tax,days>40&value>100000)
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

test<-ggplot(testdata2, aes(x ="", y = value, fill=variable))+		
	facet_wrap( ~treatment*days,nrow=2,ncol=16)+				
	geom_bar(width = 1, stat = "identity")+
	geom_text(aes(label=class,vjust=0.5),position=c("stack","jitter"),size=2)+
	theme(legend.position='none')
test


#gemeinschaftsplot mit mindestreadsanzahl
sample_subset<-subset(abs_abu_molten_mean_tax,days>40&value>80000)
test_groesser_0.5<-ggplot(sample_subset, aes(x = new_days, y = value, fill=order))+
	scale_fill_manual(values = cols)+
	facet_wrap( ~treatment,nrow=2,ncol=1)+
	geom_bar(data=subset(sample_subset,nucleic_acid=="dna"),aes(x=new_days-0.5), width = 1, stat = "identity")+
	geom_bar(data=subset(sample_subset,nucleic_acid=="cdna"),aes(x=new_days+0.5),width = 1, stat = "identity")
	#geom_text(data=subset(sample_subset,nucleic_acid=="dna"), aes(x=new_days-0.5,label=order),vjust=0.5,position="stack",size=1.8)+
	#geom_text(data=subset(sample_subset,nucleic_acid=="cdna"), aes(x=new_days+0.5,label=order),vjust=0.5,position="stack",size=1.8)
	#theme(legend.position='none')
test_groesser_0.5
ggsave(file="absolute_community_order_level.png", width=14, height=8)


#shuffled colours
#To shuffle the order of the colors, you could use
sample_subset<-subset(abs_abu_molten_mean_tax,days>40&value>50000)
sample_subset<-droplevels(sample_subset)
require(scales)
n <- length(levels(sample_subset$order)) # number of colors
cols <- hue_pal(h = c(0, 360) + 15, 
                c = 100, l = 65, 
                h.start = 0, direction = 1)(n)[order(sample(1:n, n))] # color palette in random order

test_groesser_0.5<-ggplot(sample_subset, aes(x = new_days))+
	scale_fill_manual(values = cols)+
	facet_wrap( ~treatment,nrow=2,ncol=1)+
	geom_bar(data=subset(sample_subset,nucleic_acid=="dna"),aes(x=new_days-0.5, y = value, fill=order), width = 1, stat = "identity")+
	geom_bar(data=subset(sample_subset,nucleic_acid=="cdna"),aes(x=new_days+0.5, y = value, fill=order), width = 1, stat = "identity")+
	guides(fill=guide_legend(ncol=1))+
	geom_point(data=more_cell_counts_44,aes(x=new_day,y=glyph_mg_L*5000000))+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 20))
	#theme(legend.position='none')
test_groesser_0.5
ggsave(file="absolute_community_order_level.png", width=14, height=8)				


#geom_text(data=subset(sample_subset,nucleic_acid=="dna"), aes(x=days-0.5,label=phylum),vjust=0.5,position="stack",size=1.8)
#geom_text(data=subset(sample_subset,nucleic_acid=="cdna"), aes(x=days+0.5,label=phylum),vjust=0.5,position="stack",size=1.8)

#unique labels
geom_text(aes(x = longitude, y = latitude, label = as.character(state)), 
agg.data <- aggregate(cbind(days,treatment,value)~order, data = sample_subset, FUN=mean)
 geom_text(data = agg.data, 
            aes(x = longitude, y = latitude, label = as.character(state)))
#ohne legende: +theme(legend.position='none')	

#um aus säulendiagramm ein tortendiagramm zu machen:
torten_test<-test + coord_polar("y", start=0)

#sample for stack overflow
testdata<-subset(abs_abu_molten_mean_tax,nucleic_acid=="cdna"&days>40&value>50000)
testdata2<-read.csv(file.choose(),sep=";")
test<-ggplot(testdata2, aes(x =days, y = value, fill=order))+		
	geom_bar(width = 1, stat = "identity")+
	geom_text(aes(label=order,vjust=0.5),position="stack",size=2)+
	theme(legend.position='none')
test