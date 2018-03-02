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
tnt_otu10_table<-read.csv(file.choose(),sep=";",header=TRUE)													#read in csv-table to R
tnt_otu10_table																								#display the object "tnt_otu10_table"
head(tnt_otu10_table,6)																						#show the top 6 lines (compare "tail")
names(tnt_otu10_table)																						#show the different columns ("headers")
str(tnt_otu10_table)																							#show the structure of "tnt_otu10_table"
levels(tnt_otu10_table$sample_name)																			#show the different values of the column "sample_name" in the object "tnt_otu10_table"
													
tnt_otu10_table$wholetax<-do.call(paste, c(tnt_otu10_table[c("genus","family","order","class","phylum","kingdom")],sep="_")) 	#merge taxonomic classes to one name in a column called "wholetax"
#agg_tnt_otu10_table<-aggregate(total_seq~sample_name+wholetax, data = tnt_otu10_table, sum) 					#possible step to combine parallels as mean or sum, all columns except the to-be-combined must be named

library(reshape2) 																						

tnt_otu10_table<-tnt_otu10_table[,c(2,9,1)]	 																#select columns and reorder
#forgot to shorten sample_names, samplenames should not start with numbers, also removed controls as reads were minimal (Delftia 43)																
write.csv(tnt_otu10_table,file="tnt_otu10_table.csv")
tnt_otu10_table<-read.csv(file.choose(),sep=";",header=TRUE,row.names=1)

cast_tnt_otu10_table<-dcast(tnt_otu10_table,sample_name~wholetax)
cast_tnt_otu10_table[is.na(cast_tnt_otu10_table)]<-0 																#replace NA with 0
names(cast_tnt_otu10_table)
nrow(cast_tnt_otu10_table)																					#number of rows
ncol(cast_tnt_otu10_table)																					#number of columns
row.names(cast_tnt_otu10_table)<-cast_tnt_otu10_table$sample_name													#turns "sample_name" into row names
cast_tnt_otu10_table<-cast_tnt_otu10_table[,-1]																		#remove "samples_names"
names(cast_tnt_otu10_table)
rowSums(cast_tnt_otu10_table) 																				#shows reads per sample
otu_library_sizes<-rowSums(cast_tnt_otu10_table) 																#shows reads per sample, but missing the removed singletons!
write.csv(otu_library_sizes,file="tnt_library_sizes.csv")
names(cast_tnt_otu10_table)
write.csv(cast_tnt_otu10_table,file="cast_tnt_otu10_table.csv")
tcast_tnt_otu10_table<-t(cast_tnt_otu10_table)
head(tcast_tnt_otu10_table)
write.csv(tcast_tnt_otu10_table,file="tcast_tnt_otu10_table.csv")
#calculate sum of reads per sample in excel, then copy new table the relative abundances
#read relative abundances into R
tcast_tnt_otu10_table2<-read.csv(file.choose(),row.names=1,sep=";")

#transpose back for merging
cast_tnt_otu10_table3<-t(tcast_tnt_otu10_table2)

#hier zuerst vegan einbauen, um diversity indices in meta datei einzubauen?
head(cast_tnt_otu10_table3)
#diversity analysis (pielou in excel)
h_cast_tnt_otu10_table3<-diversity(cast_tnt_otu10_table3)
s_cast_tnt_otu10_table3<-specnumber(cast_tnt_otu10_table3)
write.csv(h_cast_tnt_otu10_table3,file="h_cast_tnt_otu10_table3.csv")
write.csv(s_cast_tnt_otu10_table3,file="s_cast_tnt_otu10_table3.csv")
#create meta_data file in excel for merging
#read in meta_data
meta_data_tnt10<-read.csv(file.choose(),header=T,row.names=1,sep=";") 										#samples in rows, meta data in colums


#merge:zwei tabellen, die eine gemeinsame spalte aufweisen (hier: row.names), werden zusammengesetzt. soviele zeilen, wie in x (der erstgenannten tabelle enthalten sind, kommen in die finale datei)
tnt10_rel_abu_with_meta<-merge(cast_tnt_otu10_table3,meta_data_tnt10,by="row.names",all.x=TRUE)			
#after merging the row.names need to be adjusted
head(tnt10_rel_abu_with_meta)
row.names(tnt10_rel_abu_with_meta)<-tnt10_rel_abu_with_meta$sample_name
head(tnt10_rel_abu_with_meta)
tnt10_rel_abu_with_meta<-tnt10_rel_abu_with_meta[,(-1)] #getting rid of the first column "Row.names"
tnt10_rel_abu_with_meta<-tnt10_rel_abu_with_meta[,(-922)] #getting rid of "sample_name"
#tnt10_rel_abu_with_meta<-tnt10_rel_abu_with_meta[,c(-937:-940)] #getting rid of old diversity values from otu1 table
names(tnt10_rel_abu_with_meta)
head(tnt10_rel_abu_with_meta)
write.csv(tnt10_rel_abu_with_meta,file="tnt10_rel_abu_with_meta.csv")
#tnt10_rel_abu_with_meta_mean<-aggregate(c(1~2:659,661:665),data = tnt10_rel_abu_molten, FUN=mean)

#melt into final table with all samples and all metadata combinations
tnt10_rel_abu_molten<-melt(tnt10_rel_abu_with_meta, id=c("udemm_nr","station","station_ord","mooring","area","lat","long","cruise","cruise_date","parallel","sampling","sample","extraction_yield","nucleic_acid","library_size","shannon","max_shannon","species_richness","pielou_evenness"))
head(tnt10_rel_abu_molten,50)
tail(tnt10_rel_abu_molten,50) #check if the end of list is correct, otherwise you probably missed a column name in the melting command
str(tnt10_rel_abu_molten)			#new melt function from reshape2 turns value into character format? not always

#add taxonomic levels
tnt_tax<-read.csv(file.choose(),row.names=1,sep=";")
tnt10_rel_abu_molten_tax<-merge(tnt10_rel_abu_molten,tnt_tax,by.x="variable", by.y="row.names",all.x=TRUE)
str(tnt10_rel_abu_molten_tax)
write.csv(tnt10_rel_abu_molten_tax,file="tnt10_rel_abu_molten_tax.csv")


##########################################################################################################################plot ideas
#you can further subset table or use everything for plotting. the plotting itself also provides subsetting options
#specific otu
#vibrio
gtest<-subset(tnt10_rel_abu_molten_tax, grepl("Vibrio",family))
gg<-ggplot(gtest,aes(x=parallel,y=value,colour=nucleic_acid,fill=genus))+
	geom_bar(width = 1, stat = "identity",fill="black")+
	labs(x="parallels", y="relative abundance [%]")+
	scale_colour_manual(values=c("dna"="blue2","cdna"="brown"),
						name="Vibrio amplicons ",
						breaks=c("cdna","dna"),
						labels=c("16S rRNA","16S rRNA gene"))+
	theme(legend.position="bottom")+
	facet_wrap(~nucleic_acid*station_ord,nrow=2)
	#theme(legend.position='none')
gg	
ggsave(file="Clostridia_STV.png", width=20, height=12)

##only plot diversity
jo<-subset(meta_data_tnt10,area!="lab_control")
divtest<-ggplot(jo,aes(x=as.factor(parallel)))+
	geom_point(aes(y=shannon,colour="Shannon Diversity"),colour="black",size=3)+
	#geom_point(aes(y=species_richness,colour="Species Richness"))+
	theme(legend.text=element_text(size=11))+
	theme(panel.grid.major=element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor=element_line(colour = NA, size = 0.5))+
	theme(legend.position="bottom")+
	theme(legend.background = element_rect(fill="grey90",linetype="solid"))+
	labs(x="Parallels", y="Shannon Diversity")+	
	facet_wrap(~nucleic_acid*station_ord,nrow=2)+
	theme(legend.position='none')
divtest	
ggsave(file="shannon_diversity_otu10.png", width=14, height=8)

#gemeinschaftsplot
test_min<-ggplot(tnt10_rel_abu_molten_tax, aes(x =parallel, y = value, fill=phylum))+		
	geom_bar(width = 1, stat = "identity")+
	#geom_text(aes(x = parallel, y =value, label = class),position="stack",vjust=0.5,size=1.2)+
	#geom_point(aes(y=library_size/2500),colour="red",size=3)+
	#geom_point(aes(y=shannon*15),colour="black")+
	#geom_point(aes(y=species_richness/12),colour="yellow")+
	facet_wrap(~nucleic_acid*station_ord*cruise_date,nrow=2)+
	theme(legend.position='none')
test_min
ggsave(file="rel_community_phylum_dna_min10_2.png", width=20, height=12)	

sample_subset<-subset(tnt_rel_abu_molten_tax,value>3)
test_min<-ggplot(sample_subset, aes(x =parallel, y = value, fill=genus))+		
	geom_bar(width = 1, stat = "identity")+
	#geom_text(aes(x = parallel, y =value, label = class),position="stack",vjust=0.5,size=1.2)+
	#geom_point(aes(y=library_size/2500),colour="red",size=3)+
	#geom_point(aes(y=shannon*15),colour="black")+
	#geom_point(aes(y=species_richness/12),colour="yellow")+
	facet_wrap(~nucleic_acid*station_ord*cruise_date,nrow=2)
	#theme(legend.position='none')
test_min
ggsave(file="rel_community_genus_groesser_3_dna_min10.png", width=20, height=12)

#gemeinschaftsplot stationauswahl
station_subset<-subset(tnt10_rel_abu_molten_tax,grepl("Mo",station))
test<-ggplot(station_subset, aes(x =parallel, y = value, fill=variable))+		
	geom_bar(width = 1, stat = "identity")+
	#geom_point(aes(y=library_size/2500),colour="red",size=3)+
	#geom_point(aes(y=shannon*15),colour="black")+
	#geom_point(aes(y=species_richness/12),colour="yellow")+
	facet_wrap(~nucleic_acid*station_ord*cruise_date,nrow=1)+
	theme(legend.position='none')
test
ggsave(file="rel_community_otu_Mo_transect_min10.png", width=14, height=8)	

#geom_point(data=subset(meta_nmds_tnt10_cdna,cruise=="L17_07"),aes(x=MDS1,y=MDS2,shape=as.factor(parallel),colour=cruise),size=2.5,colour="blue") +
 # geom_point(data=subset(meta_nmds_tnt10_cdna,cruise=="L17_08"),aes(x=MDS1,y=MDS2,shape=as.factor(parallel),colour=cruise),size=4,colour="red")+
  #geom_point(data=subset(meta_nmds_tnt10_cdna,cruise=="L16_14"),aes(x=MDS1,y=MDS2,shape=as.factor(parallel),colour=cruise),size=4,colour="green")+
  #geom_point(data=subset(meta_nmds_tnt10_cdna,station=="Mo7")

#gemeinschaftsplot mit mindestreadsanzahl
sample_subset<-subset(tnt_rel_abu_molten_tax,value>2)
test<-ggplot(sample_subset, aes(x =parallel, y = value, fill=))+		
	geom_bar(width = 1, stat = "identity")+
	geom_point(aes(y=library_size/2500),colour="red",size=3)+
	geom_point(aes(y=shannon*15),colour="black")+
	geom_point(aes(y=species_richness/12),colour="yellow")+
	facet_wrap(~nucleic_acid*station_ord*cruise_date,nrow=2)+
	theme(legend.position='none')
test


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