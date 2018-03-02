#plot single contigs
####################compare functional diversity
#START HERE (if clean workspace)

##read required files and libraries
#omics_collection<-read.csv(file.choose(),sep=";")
#meta_omics<-read.csv(file.choose(),row.names=1,sep=";")
#meta_omics_twodigit<-read.csv(file.choose(),row.names=1,sep=";")
#more_cell_counts<-read.csv(file.choose(),sep=";")
#more_cell_counts_60<-subset(more_cell_counts,day >= 60)
#library(reshape)
#library(ggplot2)
#library(data.table)
#samples_omic<-c("A01","A02","A03","A04","A05","A06","A07","B08","B09","B10")

#START HERE (if workspace and libraries loaded)

#subset(omics_collection, grepl("k141_28018", gene))
#k141_28018<-subset(omics_collection, grepl("k141_28018", gene) & !grepl("k141_28018C", gene))

k141_166724<-subset(omics_collection, grepl("k141_166724", contig_id)) #& grepl("1.5.3.1", ec_number))# & !grepl("3.1.3.16", ec_number))
head(k141_166724)
nrow(k141_166724)
k141_166724single<-subset(k141_166724,grepl("soxA", gene))
k141_166724_select<-k141_166724single[,c(1,3,4,8:17)]
k141_166724_select<-droplevels(k141_166724_select)
head(k141_166724_select)
molten_k141_166724<-reshape2::melt(k141_166724_select,id=c("contig_id","ec_number","gene"))
molten_k141_166724


#meta data for plotting: treatment/control and time
merged_k141_166724_select<-merge(molten_k141_166724,meta_omics,by.x="variable",by.y="row.names",all.x=T)
head(merged_k141_166724_select)

#create variable for adjusting different ranges of data series
readfactor<-max(merged_k141_166724_select$value)
cellfactor<-max(merged_k141_166724_select$value)

#plotting only relative contig abundance
goi3<-expression(paste(,italic("k141_166724")," metagenomic relative contig abundance in water column"))
k141_166724_abu<-ggplot(merged_k141_166724_select, aes(x=day))+
	geom_bar(width=1.4,stat="identity",aes(x=day, y=value,colour="contig_read_\nabundance\n",fill=contig_id))+
	geom_vline(aes(xintercept=69.5),linetype="dashed", size=0.8)+
	geom_point(data=more_cell_counts_60, aes(y=cells_ml*readfactor/(100000000),colour="cell counts\n"),fill="grey2",shape=21,size=3.5,alpha=0.8)+
	facet_wrap(~treatment,nrow=2)+
	ggtitle(goi3)+
	scale_fill_discrete(name = "contig ID")+
	scale_colour_discrete(name = "data series")+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 20))+
	theme(plot.title = element_text(size=18,face="bold",hjust = 0.5))+
	theme(axis.title = element_text(size=12,face="bold"))+
	theme(axis.title.y = element_text(angle=0,vjust=0.5))+
	theme(axis.text=element_text(size=12))+
	theme(legend.title=element_text(size=13,face="bold"))+
	theme(legend.text=element_text(size=11))+
	xlab("days")+
	#theme(legend.position="none")+
	#guides(fill=FALSE)+
	ylab("rpm per\ncontig")
k141_166724_abu
ggsave(file="k141_166724_abundance.jpg", width=14, height=8)