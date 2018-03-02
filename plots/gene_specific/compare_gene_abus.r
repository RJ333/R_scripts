####################compare functional diversity
#START HERE (if clean workspace)

##read required files and libraries
#kallisto_prokka_concoct2_metaxa_checkm<-read.csv(file.choose(),sep=";")
#meta_omics<-read.csv(file.choose(),row.names=1,sep=";")
#meta_omics_twodigit<-read.csv(file.choose(),row.names=1,sep=";")
#more_cell_counts<-read.csv(file.choose(),sep=";")
#more_cell_counts_60<-subset(more_cell_counts,day >= 60)
#library(reshape)
#library(ggplot2)
#library(data.table)
#samples_omic<-c("A01","A02","A03","A04","A05","A06","A07","B08","B09","B10")

#START HERE (if workspace and libraries loaded)

#subset(kallisto_prokka_concoct2_metaxa_checkm, grepl("phnACGHIJLMXW_pphA_soxAB_thiO", genes))
phnACGHIJLMXW_pphA_soxAB_thiO<-subset(kallisto_prokka_concoct2_metaxa_checkm, grepl("phnA|phnC|phnG|phnH|phnI|phnJ|phnL|phnM|phnW|phnX|pphA|soxA|soxB|thiO", genes) & !grepl("rpoBC", genes))
#phnACGHIJLMXW_pphA_soxAB_thiO<-subset(kallisto_prokka_concoct2_metaxa_checkm, grepl("phnG|phnH|phnI|phnJ|phnL|phnM", genes))
head(phnACGHIJLMXW_pphA_soxAB_thiO)
nrow(phnACGHIJLMXW_pphA_soxAB_thiO)
phnACGHIJLMXW_pphA_soxAB_thiO_select<-phnACGHIJLMXW_pphA_soxAB_thiO[,c(1,5,6,8,10,12,14,16,18,20,22,24,26)]
phnACGHIJLMXW_pphA_soxAB_thiO_select<-droplevels(phnACGHIJLMXW_pphA_soxAB_thiO_select)
head(phnACGHIJLMXW_pphA_soxAB_thiO_select)
molten_phnACGHIJLMXW_pphA_soxAB_thiO<-melt(phnACGHIJLMXW_pphA_soxAB_thiO_select,id=c("contig_id","ec_number","genes"))
head(molten_phnACGHIJLMXW_pphA_soxAB_thiO)
nrow(molten_phnACGHIJLMXW_pphA_soxAB_thiO)

#meta data for plotting: treatment/control and time
merged_phnACGHIJLMXW_pphA_soxAB_thiO_select<-merge(molten_phnACGHIJLMXW_pphA_soxAB_thiO,meta_omics,by.x="variable",by.y="row.names",all.x=T)
head(merged_phnACGHIJLMXW_pphA_soxAB_thiO_select)
merged_phnACGHIJLMXW_pphA_soxAB_thiO_select2<-merged_phnACGHIJLMXW_pphA_soxAB_thiO_select[order(merged_phnACGHIJLMXW_pphA_soxAB_thiO_select[,2]),]

###remove zero count entries to correctly count abundant genes and contigs; drop levels in subset
phnACGHIJLMXW_pphA_soxAB_thiO_nonzero<-subset(merged_phnACGHIJLMXW_pphA_soxAB_thiO_select,value!=0.0000000)
phnACGHIJLMXW_pphA_soxAB_thiO_nonzero <- droplevels(phnACGHIJLMXW_pphA_soxAB_thiO_nonzero)
phnACGHIJLMXW_pphA_soxAB_thiO_unique1<-as.data.table(phnACGHIJLMXW_pphA_soxAB_thiO_nonzero)[, unique_contigs := length(unique(gsub("\\..*$", "", contig_id))), by = variable][]
phnACGHIJLMXW_pphA_soxAB_thiO_unique2<-as.data.table(phnACGHIJLMXW_pphA_soxAB_thiO_unique1)[, unique_genes := length(gsub("_.*$", "", genes)), by = variable][]
#create variable for adjusting different ranges of data series
readfactor<-max(phnACGHIJLMXW_pphA_soxAB_thiO_unique2$value)
cellfactor<-max(phnACGHIJLMXW_pphA_soxAB_thiO_unique2$unique_contigs)

#plotting 
goi<-expression(paste(,italic("phnACGHIJLMXW_pphA_soxAB_thiO")," metagenomic read distribution in water column"))
phnACGHIJLMXW_pphA_soxAB_thiO_compare<-ggplot(phnACGHIJLMXW_pphA_soxAB_thiO_unique2, aes(x=day))+
	geom_vline(aes(xintercept=70.0),linetype="dashed", size=0.8)+
	geom_bar(width=2.8,stat="summary",fun.y="sum",position="dodge",aes(x=day-0.5, y=value,colour="contig_read_\nabundance\n",fill=gsub("_.*$", "",genes)))+
	geom_point(data=more_cell_counts_60, aes(y=cells_ml/30000,colour="cell counts\n"),fill="grey2",shape=21,size=3.5,alpha=0.8)+
	facet_wrap(~treatment,nrow=2)+
	ggtitle(goi)+
	scale_fill_discrete(name = "compared genes")+
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
	ylab("contigs per\nmillion")
phnACGHIJLMXW_pphA_soxAB_thiO_compare
ggsave(file="phnACGHIJLMXW_pphA_soxAB_thiO_abu_compare.jpg", width=14, height=8)

#gsub("_.*$", "", 