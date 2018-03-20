# open otu tax breakdown file with excel
# make columns wide and create 6 empty columns after sample name 
# data -> text in columns -> separated via tabstop
# name tax path columns with genus, family and so on

# clean up table:
# from the count columns, only the summed reads are important
# remove non-bacteria, sort for "kingdom", remove "no relative"
# replace " " and "." with "-" in OTU names
# remove singletons (e.g. sort by reads)
# save as csv

# open R
# set working folder
# load workspace
# download packages: vegan, ggplot2, reshape2


# read in csv-table to R	
otu_table<-read.csv(file.choose(),sep=";",header=TRUE)	# or sep=","  depending on excel												
																		
													
otu_table$wholetax<-do.call(paste, c(otu_table[c("genus","family","order","class","phylum","kingdom")],sep="_")) 	#merge taxonomic classes to one name in a column called "wholetax"
#agg_otu_table<-aggregate(total_seq~sample_name+wholetax, data = otu_table, sum) 					#possible step to combine parallels as mean or sum, all columns except the to-be-combined must be named

library(reshape2) 																						


#select columns with counts, name, and sample name and reorder
otu_table<-otu_table[,c(2,9,1)]	 	

#turn table into wide format															
cast_otu_table<-dcast(otu_table,sample_name~wholetax)

#replace NA with 0
cast_otu_table[is.na(cast_otu_table)]<-0 																

# check column names
names(cast_otu_table)

# turns "sample_name" into row names
row.names(cast_otu_table)<-cast_otu_table$sample_name													

# remove column "sample name"
cast_otu_table<-cast_otu_table[,-1] 