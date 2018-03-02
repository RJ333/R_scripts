####################compare functional diversity
#START HERE (if clean workspace)

##read required files and libraries
#kallisto_prokka_concoct2_metaxa_checkm<-read.csv(file.choose(),sep=";")
#meta_omics<-read.csv(file.choose(),row.names=1,sep=";")
#meta_omics_twodigit<-read.csv(file.choose(),row.names=1,sep=";")
#more_cell_counts<-read.csv(file.choose(),sep=";")
#more_cell_counts_44<-subset(more_cell_counts,day >= 44)
#more_cell_counts_60<-subset(more_cell_counts,day >= 60)
#library(reshape)
#library(ggplot2)
#library(data.table)
#samples_omic<-c("A01","A02","A03","A04","A05","A06","A07","B08","B09","B10")

#START HERE (if workspace and libraries loaded)

#nrow(subset(kallisto_prokka_concoct2_metaxa_checkm, grepl("phnA", genes)))
phnA<-subset(kallisto_prokka_concoct2_metaxa_checkm, grepl("phnA", genes))
head(phnA)
nrow(phnA)
phnA_select<-phnA[,c(1,5,6,8,10,12,14,16,18,20,22,24,26)]
phnA_select<-droplevels(phnA_select)
head(phnA_select)
molten_phnA<-melt(phnA_select,id=c("contig_id","ec_number","genes"))
head(molten_phnA)
nrow(molten_phnA)

#meta data for plotting: treatment/control and time
merged_phnA_select<-merge(molten_phnA,meta_omics,by.x="variable",by.y="row.names",all.x=T)
head(merged_phnA_select)
merged_phnA_select2<-merged_phnA_select[order(merged_phnA_select[,2]),]

###remove zero count entries to correctly count abundant genes and contigs; drop levels in subset
phnA_nonzero<-subset(merged_phnA_select,value!=0.0000000)
phnA_nonzero <- droplevels(phnA_nonzero)
phnA_unique1<-as.data.table(phnA_nonzero)[, unique_contigs := length(unique(gsub("\\..*$", "", contig_id))), by = variable][]
phnA_unique2<-as.data.table(phnA_unique1)[, unique_genes := length(gsub("_.*$", "", genes)), by = variable][]
#create variable for adjusting different ranges of data series
readfactor<-max(phnA_unique2$value)
cellfactor<-max(phnA_unique2$unique_contigs)

#plotting 
goi<-expression(paste(,italic("phnA")," metagenomic read distribution in water column"))
phnA_iso<-ggplot(phnA_unique2, aes(x=day))+
	geom_line(data=more_cell_counts_60, aes(y=glyph_mg_L*(cellfactor/13),colour="glyphosate\nconcentration\n"),alpha=0.8,linetype="solid", size=1)+
	geom_line(data=more_cell_counts_60, aes(y=glyph_theor*(cellfactor/13),colour="glyphosate\ndilution\n"),alpha=0.5,linetype="F1", size=1)+
	geom_vline(aes(xintercept=70.0),linetype="dashed", size=0.8)+
	geom_bar(width=0.5,stat="summary",fun.y="mean",aes(x=day+0.5, y=unique_contigs,colour="sum_unique_\ncontigs\n"),alpha=1,fill="yellow1")+
	geom_bar(width=0.5,stat="summary",fun.y="mean",aes(y=unique_genes,colour="sum_unique_\ngenes\n"),alpha=1,fill="orange")+
	geom_bar(width=1.4,stat="identity",aes(x=day-1, y=value/(readfactor/20),colour="contig_read_\nabundance\n",fill=contig_id))+
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
	guides(fill=FALSE)
	ylab("relative\nunits")
phnA_iso
#ggsave(file="phnA_gene_isoforms.jpg", width=14, height=8)