#otu table R:
#open file with excel
#make columns wide and create 6 empty columns after sample name 
#data -> text in columns -> separated via tabstop
#name path columns

#clean up table: 
#remove non-bacteria, sort for "kingdom", remove "no relative"
#replace " " and "." with "-" in OTU names
#remove singletons
#second table, remove all <5 reads
#adjust sample name: extra column, =Teil(cell;startcharacter;numberofcharacters) --> =Teil(A3;6;7)
#save as csv

#open R
#set working folder
#load workspace
#download packages: vegan, ggplot2, reshape2


#first prepare data with only singletons removed
tnt_otu1_table<-read.csv(file.choose(),sep=";",header=TRUE)													#read in csv-table to R
tnt_otu1_table																								#display the object "tnt_otu1_table"
head(tnt_otu1_table,6)																						#show the top 6 lines (compare "tail")
names(tnt_otu1_table)																						#show the different columns ("headers")
str(tnt_otu1_table)																							#show the structure of "tnt_otu1_table"
levels(tnt_otu1_table$sample_name)																			#show the different values of the column "sample_name" in the object "tnt_otu1_table"
													
tnt_otu1_table$wholetax<-do.call(paste, c(tnt_otu1_table[c("genus","family","order","class","phylum","kingdom")],sep="_")) 	#merge taxonomic classes to one name in a column called "wholetax"
#agg_tnt_otu1_table<-aggregate(total_seq~sample_name+wholetax, data = tnt_otu1_table, sum) 					#possible step to combine parallels as mean or sum, all columns except the to-be-combined must be named

library(reshape2) 																						

tnt_otu1_table<-tnt_otu1_table[,c(2,9,1)]	 																#select columns and reorder
#forgot to shorten sample_names, samplenames should not start with numbers, also removed controls as reads were minimal (Delftia 43)																
write.csv(tnt_otu1_table,file="tnt_otu1_table.csv")
tnt_otu1_table<-read.csv(file.choose(),sep=";",header=TRUE)

cast_tnt_otu1_table<-dcast(tnt_otu1_table,sample_name~wholetax)
cast_tnt_otu1_table[is.na(cast_tnt_otu1_table)]<-0 																#replace NA with 0
names(cast_tnt_otu1_table)
nrow(cast_tnt_otu1_table)																					#number of rows
ncol(cast_tnt_otu1_table)																					#number of columns
row.names(cast_tnt_otu1_table)<-cast_tnt_otu1_table$sample_name													#turns "sample_name" into row names
cast_tnt_otu1_table<-cast_tnt_otu1_table[,-1]																		#remove "samples_names"
names(cast_tnt_otu1_table)
rowSums(cast_tnt_otu1_table) 																				#shows reads per sample
otu_library_sizes<-rowSums(cast_tnt_otu1_table) 																#shows reads per sample, but missing the removed singletons!
write.csv(otu_library_sizes,file="tnt_library_sizes.csv")
names(cast_tnt_otu1_table)
write.csv(cast_tnt_otu1_table,file="cast_tnt_otu1_table.csv")
tcast_tnt_otu1_table<-t(cast_tnt_otu1_table)
head(tcast_tnt_otu1_table)
write.csv(tcast_tnt_otu1_table,file="tcast_tnt_otu1_table.csv")
#calculate sum of reads per sample in excel, then copy new table the relative abundances
#read relative abundances into R
tcast_tnt_otu1_table2<-read.csv(file.choose(),row.names=1,sep=";")

#transpose back for merging
cast_tnt_otu1_table3<-t(tcast_tnt_otu1_table2)
head(cast_tnt_otu1_table3)
#create meta_data file in excel for merging
#read in meta_data
meta_data_tnt<-read.csv(file.choose(),header=T,row.names=1,sep=";") 										#samples in rows, meta data in colums
str(meta_data_tnt)
levels(meta_data_tnt$cellcounts) <- c(levels(meta_data_tnt$cellcounts),"NA")
meta_data_tnt$cellcounts[is.na(meta_data_tnt$cellcounts)] <- "NA"
str(meta_data_tnt) #"NA" in cell counts are now strings, which is necessary for aggregate function later

#meta_data_tnt$time<-as.factor(meta_data_tnt$time)																#depending on the aims it might be necessary to change the classes of the "numeric" columns to "factor" columns
#meta_data_tnt$days<-as.factor(meta_data_tnt$days)
#meta_data_tnt$parallel<-as.factor(meta_data_tnt$parallel)
#str(meta_data_tnt)

#merge:zwei tabellen, die eine gemeinsame spalte aufweisen (hier: row.names), werden zusammengesetzt. soviele zeilen, wie in x (der erstgenannten tabelle enthalten sind, kommen in die finale datei)
tnt_rel_abu_with_meta<-merge(cast_tnt_otu1_table3,meta_data_tnt,by="row.names",all.x=TRUE)			
#after merging the row.names need to be adjusted
head(tnt_rel_abu_with_meta)
row.names(tnt_rel_abu_with_meta)<-tnt_rel_abu_with_meta$sample_name
head(tnt_rel_abu_with_meta)
tnt_rel_abu_with_meta<-tnt_rel_abu_with_meta[,(-1)] #getting rid of the first column "Row.names"
tnt_rel_abu_with_meta<-tnt_rel_abu_with_meta[,(-1712)] #getting rid of "sample_name"
names(tnt_rel_abu_with_meta)
head(tnt_rel_abu_with_meta)
write.csv(tnt_rel_abu_with_meta,file="tnt_rel_abu_with_meta.csv")
#tnt_rel_abu_with_meta_mean<-aggregate(c(1~2:659,661:665),data = tnt_rel_abu_molten, FUN=mean)

#melt into final table with all samples and all metadata combinations
tnt_rel_abu_molten<-melt(tnt_rel_abu_with_meta, id=c("udemm_nr","station","station_ord","mooring","area","lat","long","cruise","cruise_date","parallel","sampling","sample","extraction_yield","nucleic_acid","library_size","shannon","max_shannon","species_richness","pielou_evenness"))
head(tnt_rel_abu_molten,50)
tail(tnt_rel_abu_molten,50) #check if the end of list is correct, otherwise you probably missed a column name in the melting command
str(tnt_rel_abu_molten)			#new melt function from reshape2 turns value into character format? not always

#add taxonomic levels
tnt_tax<-read.csv(file.choose(),row.names=1,sep=";")
tnt_rel_abu_molten_tax<-merge(tnt_rel_abu_molten,tnt_tax,by.x="variable", by.y="row.names",all.x=TRUE)
str(tnt_rel_abu_molten_tax)
#remove exraction controls
tnt_rel_abu_molten_tax<-subset(tnt_rel_abu_molten_tax,area != "lab_control")
write.csv(tnt_rel_abu_molten_tax,file="tnt_rel_abu_molten_tax.csv")




##########################################################################################################################plot ideas
#you can further subset table or use everything for plotting. the plotting itself also provides subsetting options
sample_subset<-subset(abs_abu_molten_mean_tax,days>40&value>100000)
more_cell_counts<-read.csv(file.choose(),sep=";")
#specific otu
#gallaeci
gtest<-subset(tnt_rel_abu_molten_tax, grepl("Photobacterium",genus))
gg<-ggplot(gtest,aes(x=parallel,y=value,colour=nucleic_acid))+
	geom_bar(width = 1, stat = "identity")+
	facet_wrap(~nucleic_acid*station_ord,nrow=2)
gg	

##only plot diversity
divtest<-ggplot(meta_data_tnt,aes(x=days))+
	geom_point(data=meta_data_tnt,aes(y=meta_data_tnt$shannon*100,colour=nucleic_acid))+
	geom_point(data=meta_data_tnt,aes(y=meta_data_tnt$richness,shape=nucleic_acid))+
	facet_wrap(~treatment*habitat,nrow=2,ncol=2)
divtest	


#gemeinschaftsplot
test<-ggplot(tnt_rel_abu_molten_tax, aes(x =parallel, y = value, fill=family))+		
	geom_bar(width = 1, stat = "identity")+
	geom_point(aes(y=library_size/2500),colour="red",size=3)+
	geom_point(aes(y=shannon*15),colour="black")+
	geom_point(aes(y=species_richness/12),colour="yellow")+
	facet_wrap(~nucleic_acid*station_ord*cruise_date,nrow=2)+
	theme(legend.position='none')
test
#gemeinschaftsplot mit mindestreadsanzahl
sample_subset<-subset(tnt_rel_abu_molten_tax,value>0.02)
test<-ggplot(sample_subset, aes(x =parallel, y = value, fill=genus))+		
	geom_bar(width = 1, stat = "identity")+
	geom_point(aes(y=library_size/2500),colour="red",size=3)+
	geom_point(aes(y=shannon*15),colour="black")+
	geom_point(aes(y=species_richness/12),colour="yellow")+
	facet_wrap(~nucleic_acid*station_ord*cruise_date,nrow=2)+
	theme(legend.position='none')
test

#gemeinschaftsplot mit mindestreadsanzahl
sample_subset<-subset(tnt_rel_abu_molten_tax,value>80000)
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