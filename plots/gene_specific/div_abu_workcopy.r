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

#subset(omics_collection, grepl("soxB", gene))
#soxB<-subset(omics_collection, grepl("soxB", gene) & !grepl("soxBC", gene))

soxB<-subset(omics_collection, grepl("soxB", gene))# & grepl("1.3.7.8", ec_number))# & !grepl("3.1.3.16", ec_number))
head(soxB)
nrow(soxB)
soxB_select<-soxB[,c(1,3,4,8:17)]
soxB_select<-droplevels(soxB_select)
head(soxB_select)
levels(soxB_select$ec_number)
#check that only one ec number is selected (gene names can be ambiguous, eg soxB, soxB)
droplevels(soxB_select$ec_number)
droplevels(soxB_select$gene)
molten_soxB<-melt(soxB_select,id=c("contig_id","ec_number","gene"))
head(molten_soxB)
nrow(molten_soxB)

#meta data for plotting: treatment/control and time
merged_soxB_select<-merge(molten_soxB,meta_omics,by.x="variable",by.y="X",all.x=T)
head(merged_soxB_select)
merged_soxB_select2<-merged_soxB_select[order(merged_soxB_select[,2]),]

###remove zero count entries to correctly count abundant gene and contigs; drop levels in subset
soxB_nonzero<-subset(merged_soxB_select,value!=0.0000000)
soxB_nonzero <- droplevels(soxB_nonzero)
soxB_unique1<-as.data.table(soxB_nonzero)[, unique_contigs := as.numeric(length(unique(gsub("\\..*$", "", contig_id)))), by = variable][]
soxB_unique2<-as.data.table(soxB_unique1)[, unique_gene := length(gsub("_.*$", "", gene)), by = variable][]
soxB_agg<-aggregate(value~variable+ec_number+gsub("_.*$", "", gene)+time+day+treatment+cells_ml+glyph_mg_L+unique_contigs+unique_gene,data=soxB_unique2,sum)
soxB_ratio<-as.data.table(soxB_agg)[,ratio := unique_contigs/value, by = variable][]
#aggregate contig counts "sum"
#setting ratio to zero?
#new column (ratio*?)

#create variable for adjusting different ranges of data series
readfactor<-max(soxB_unique2$value)
cellfactor<-max(soxB_unique2$unique_contigs)

#plotting only relative contig abundance
goi3<-expression(paste(,italic("soxB")," metagenomic relative contig abundance in water column"))
soxB_abu<-ggplot(soxB_unique2, aes(x=day))+
	geom_bar(width=1.4,stat="identity",aes(x=day, y=value,colour="contig_read_\nabundance\n",fill=contig_id))+
	geom_vline(aes(xintercept=69.5),linetype="dashed", size=0.8)+
	geom_point(data=more_cell_counts_60, aes(y=cells_ml*readfactor/(100000000),colour="cell counts\n"),fill="grey2",shape=21,size=3.5,alpha=0.8)+
	facet_wrap(~treatment,nrow=2)+
	ggtitle(goi3)+
	scale_fill_discrete(name = "??")+
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
	guides(fill=FALSE)+
	ylab("rpm per\ncontig")
soxB_abu
ggsave(file="soxB_abundance.jpg", width=14, height=8)

#create variable for adjusting different ranges of data series
readfactor<-max(soxB_unique2$value)
cellfactor<-max(soxB_unique2$unique_contigs)
#plotting diversity and scaled relative contig abundance
goi<-expression(paste(,italic("soxB")," gene richness in water column"))
soxB_iso<-ggplot(soxB_unique2, aes(x=day))+
	geom_bar(width=1,stat="summary",fun.y="mean",aes(y=unique_gene,colour="sum_unique_\ngene\n"),alpha=1,fill="orange")+
	geom_line(data=more_cell_counts_60, aes(y=glyph_mg_L*(cellfactor/13),colour="glyphosate\nconcentration\n"),alpha=0.8,linetype="solid", size=1)+
	geom_line(data=more_cell_counts_60, aes(y=glyph_theor*(cellfactor/13),colour="glyphosate\ndilution\n"),alpha=0.5,linetype="F1", size=1)+
	geom_vline(aes(xintercept=70.0),linetype="dashed", size=0.8)+
	geom_point(data=more_cell_counts_60, aes(y=cells_ml*cellfactor/(50000000*1.5),colour="cell counts\n"),fill="grey2",shape=21,size=3.5,alpha=0.8)+
	facet_wrap(~treatment,nrow=2)+
	ggtitle(goi)+
	scale_fill_discrete(name = "??")+
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
	guides(fill=FALSE)+
	ylab("gene\nrichness")
soxB_iso
ggsave(file="soxB_richness.jpg", width=14, height=8)

###this plot only contains the "factorized" unique contig amount
#plot ratio and unique contigs
#goi2<-expression(paste(,italic("soxB")," diversity to contig abundance ratio in water column"))
#soxB_ratio_plot<-ggplot(soxB_ratio, aes(x=day))+
	#geom_vline(aes(xintercept=70.0),linetype="dashed", size=0.8)+
	#geom_bar(width=0.5,stat="identity",aes(x=day-0.5, y=ratio*unique_contigs,colour="ratio *\nunique\ncontigs"),alpha=1,fill="green3")+
	#facet_wrap(~treatment,nrow=2)+
	#ggtitle(goi2)+
	#scale_fill_discrete(name = "??")+
	#scale_colour_discrete(name = "data series")+
	#scale_x_continuous(breaks = scales::pretty_breaks(n = 20))+
	#theme(plot.title = element_text(size=18,face="bold",hjust = 0.5))+
	#theme(axis.title = element_text(size=12,face="bold"))+
	#theme(axis.title.y = element_text(angle=0,vjust=0.5))+
	#theme(axis.text=element_text(size=12))+
	#theme(legend.title=element_text(size=13,face="bold"))+
	#theme(legend.text=element_text(size=11))+
	#xlab("days")+
	#theme(legend.position="none")+
	#guides(fill=FALSE)+
	#ylab("unique\ncontigs")
#soxB_ratio_plot
#ggsave(file="soxB_ratio.jpg", width=14, height=8)

#plot ratio and unique contigs
#goi2<-expression(paste(,italic("soxB")," diversity to contig abundance ratio in water column"))
#soxB_ratio_plot<-ggplot(soxB_ratio, aes(x=day))+
	#geom_vline(aes(xintercept=70.0),linetype="dashed", size=0.8)+
	#geom_bar(width=0.5,stat="identity",aes(x=day+0.5, y=unique_contigs,colour="sum_unique_\ncontigs\n"),alpha=1,fill="yellow1")+
	#geom_bar(width=0.5,stat="identity",aes(x=day-0.5, y=ratio*1000,colour="ratio\n"),alpha=1,fill="green3")+
	#facet_wrap(~treatment,nrow=2)+
	#ggtitle(goi2)+
	#scale_fill_discrete(name = "??")+
	#scale_colour_discrete(name = "data series")+
	#scale_x_continuous(breaks = scales::pretty_breaks(n = 20))+
	#theme(plot.title = element_text(size=18,face="bold",hjust = 0.5))+
	#theme(axis.title = element_text(size=12,face="bold"))+
	#theme(axis.title.y = element_text(angle=0,vjust=0.5))+
	#theme(axis.text=element_text(size=12))+
	#theme(legend.title=element_text(size=13,face="bold"))+
	#theme(legend.text=element_text(size=11))+
	#xlab("days")+
	#theme(legend.position="none")+
	#guides(fill=FALSE)+
	#ylab("counts per\ncontig")
#soxB_ratio_plot