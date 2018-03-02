####plot mit total gene richness, species richness und gene richness und contig richness? phnM?



all_unique_contigs<-omics_collection[which(!duplicated(omics_collection[,"contig_id"])),]
nrow(all_unique_contigs)
all_unique_contigs2<-all_unique_contigs[,c(1,3,4,8:17)]

molten_all_unique_contigs<-melt(all_unique_contigs2,id=c("contig_id","ec_number","gene"))
#meta_omics<-read.csv(file.choose(),row.names=1,sep=";")
merged_all_unique_contigs<-merge(molten_all_unique_contigs,meta_omics,by.x="variable",by.y="row.names",all.x=T)

###remove zero count entries to correctly count abundant gene and contigs; drop levels in subset
all_unique_nonzero<-subset(merged_all_unique_contigs,value!=0.0000000)
all_unique_nonzero <- droplevels(all_unique_nonzero)
all_unique1<-as.data.table(all_unique_nonzero)[, unique_contigs := length(unique(gsub("\\..*$", "", contig_id))), by = variable][]
all_unique2<-as.data.table(all_unique1)[, unique_gene := length(gsub("_.*$", "", gene)), by = variable][]
all_unique_agg<-aggregate(value~variable+ec_number+gsub("_.*$", "", gene)+time+day+treatment+cells_ml+glyph_mg_L+unique_contigs+unique_gene,data=all_unique2,sum)
all_unique_ratio<-as.data.table(all_unique_agg)[,ratio := unique_contigs/value, by = variable][]

##plot all
#all_plot<-subset(merged_all_unique_contigs, grepl("all", gene)& !grepl("allC", gene))
all_unique_contigs_plot<-merged_all_unique_contigs[order(merged_all_unique_contigs[,1]),]
nrow(all_unique_contigs_plot)

div_dna_water<-subset(diversity_16S,habitat=="water" & nucleic_acid == "dna" & days >= 60)
#Gallaecimonas_water_dna_plot<-subset(final_table_tax, habitat == "water" & genus == "Gallaecimonas" &nucleic_acid=="dna"& days >= 60)

all_unique_contigs_barplot<-ggplot(all_unique_contigs_plot, aes(x = day,y=total_contig_richness))+						
	#geom_bar(width = 1.7, stat = "identity")+
	geom_bar(data=all_unique2, width=1,stat="summary",fun.y="mean",aes(x=day, y=unique_contigs,colour="contig\nrichness\n"),alpha=1,fill="yellow1")+
	geom_line(aes(colour="total_contig_richness\n"),size=1.2)+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 20))+
	#ylim(0,1000000)+
	#geom_text(data=subset(all_unique_contigs_plot, value > 50),aes(label=contig_id),size=rel(2.7), position = position_stack(vjust = 0))+
	#geom_label_repel(data=subset(all_unique_contigs_plot, value > 50),position="stack",min.segment.length = unit(0.1, "lines"),aes(label=gene,size=1.4))+
	ggtitle("contig, gene and species richness")+
	theme(plot.title = element_text(size=18,face="bold",hjust = 0.5))+
	theme(axis.title = element_text(size=12,face="bold"))+
	theme(axis.title.y = element_text(angle=90,vjust=0.5))+
	theme(axis.text=element_text(size=12))+
	theme(legend.title=element_text(size=13,face="bold"))+
	theme(legend.text=element_text(size=11))+
	geom_vline(aes(xintercept=70),linetype="dashed", size=1.2)+
	xlab("day")+
	ylab("relative units")+
	geom_point(data=more_cell_counts_60,aes(y=cells_ml/900,colour="cell_counts\n"),size=3,alpha=1)+
	geom_line(aes(y=total_gene_richness*25,colour="total_gene_richness\n"))+
	geom_line(aes(y=total_contig_richness_with_annotation,colour="contig_richness\nwith_annotation\n"))+
	geom_line(aes(y=total_contig_richness_without_annotation,colour="contig_richness\nwithout_annotation\n"))+
	#geom_point(data=div_dna_water,aes(x=days,y=richness*1000,colour="species_richness\n"))+
	stat_summary(data=div_dna_water,aes(x=days,y=richness*1000,colour="species_richness\n"),fun.y="mean", geom="line", size=1.2)+
	#geom_point(data=Gallaecimonas_water_dna_plot,aes(x=days, y=value*30000,colour="Gallaecimonas sp.\nDNA"))+
	#stat_summary(data=Gallaecimonas_water_dna_plot,aes(x=days, y=value*30000,colour="Gallaecimonas sp.\nDNA"),fun.y="mean", geom="line", size=1.5)+
	#theme(legend.position="none")+
	facet_wrap( ~treatment,nrow=2)
all_unique_contigs_barplot
 ggsave(file="contig_and_gene_and_species_richness.jpg", width=14, height=8)
 