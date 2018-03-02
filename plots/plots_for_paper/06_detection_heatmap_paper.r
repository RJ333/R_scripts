comp<-read.csv(file.choose(),sep=";") #im long format
cols1<-c("yellow","red","blue")
ggplot(comp,aes(x=days))+ #but label should be days
	geom_tile(data=subset(comp,method=="16S"&value!=1),aes(y=method,fill="16S"),colour="white")+
	geom_tile(data=subset(comp,method=="cells"&value!=1),aes(y=method,fill="cells"),colour="white")+
	geom_tile(data=subset(comp,method=="metagenome"&value!=1),aes(y=method,fill="metagenome"),colour="white")+
	geom_tile(data=subset(comp,value==1),aes(y=method,fill="no detectable\nreaction"),colour="white")+
	geom_vline(aes(xintercept=6.5),linetype="dotted", size=2)+
	coord_fixed(ratio = 0.5)+
	scale_fill_manual(values=c("16S"="orange","cells"="red","metagenome"="blue","no detectable\nreaction"="black"),
							name="blabla",
							breaks=c("16S","cells","metagenome","no detectable\nreaction"),
							labels=c("Amplicon data","Cell counts","Metagenomic data","no detectable\nreaction"))+
	#scale_x_continuous(breaks = scales::pretty_breaks(n = 20),lim=c(4, 19))
	ggsave(file="compare_approaches_paper.jpg", width=14, height=8)

##########################
comp2<-read.csv(file.choose(),sep=";") #16S workspace, final_plot_data2

comp3<-read.csv(file.choose(),sep=";",row.names=1)
comp3<-comp3[,-8]
comp3_melt<-melt(comp3,id="new_days")
comp3_melt$variable
comp3_melt$variable2<-factor(comp3_melt$variable,labels=c("CCA plot","increased species richness (DNA)","changed diversity indices","increased cell numbers","proposed memory effect","single OTU responses"))
comp4<-subset(comp3_melt,variable!="CCA_plot")

####################################clean version
ggplot(comp4,aes(x=new_days))+
	stat_bin_2d(geom="tile",bins=78,data=subset(comp4,variable=="significant_cell_counts"&value!=1),aes(y=variable,fill="significant_cell_counts"),colour="black")+
	stat_bin_2d(geom="tile",bins=78,data=subset(comp4,variable!="glyphosate_conc"&value==1),aes(y=variable,fill="no detectable reaction"),colour="black")+
	stat_bin_2d(geom="tile",bins=78,data=subset(comp4,variable=="species_richness_dna"&value!=1),aes(y=variable,fill="species_richness_dna"),colour="black")+
	stat_bin_2d(geom="tile",bins=78,data=subset(comp4,variable=="pielou.shannon"&value!=1),aes(y=variable,fill="pielou_shannon"),colour="black")+
	stat_bin_2d(geom="tile",bins=78,data=subset(comp4,variable=="gene_instance_remainder"&value!=1),aes(y=variable,fill="gene_instance_remainder"),colour="black")+
	stat_bin_2d(geom="tile",bins=78,data=subset(comp4,variable=="single_otu_plots"&value==7),aes(y=variable,fill="single_otu_plots_most"),colour="black")+
	stat_bin_2d(geom="tile",bins=78,data=subset(comp4,variable=="single_otu_plots"&value==8),aes(y=variable,fill="single_otu_plots_some"),colour="black")+
	geom_vline(aes(xintercept=0),colour="black",linetype="dotted", size=1.5)+
	theme(axis.text.y=element_blank(),
			axis.ticks.y=element_blank(),
			axis.title.y=element_blank(),
			panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
			panel.grid.minor=element_blank(),plot.background=element_blank())+
	theme(axis.title = element_text(size=20,face="bold"))+
	theme(axis.text=element_text(size=17,face="bold"))+
	theme(legend.position="none")+
	xlab("Days")+
	coord_fixed(ratio = 5)+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
	scale_fill_manual(values=c("significant_cell_counts"="grey25","no detectable\nreaction"="white","species_richness_dna"="grey25","pielou_shannon"="grey25","gene_instance_remainder"="grey25","single_otu_plots_most"="grey25","single_otu_plots_some"="grey50"),
							name="Approaches")
							ggsave(file="comparison_detection_length.png", width=20, height=9)
								
#################################dev version
ggplot(comp4,aes(x=new_days))+ #but label should be days
	stat_bin_2d(geom="tile",bins=78,data=subset(comp4,variable=="significant_cell_counts"&value!=1),aes(y=variable,fill="significant_cell_counts"),colour="black")+
	#stat_bin_2d(geom="tile",bins=78,data=subset(comp4,variable=="CCA_plot"&value==2),aes(y=variable,fill="CCA_plot_strong"),colour="black")+
	stat_bin_2d(geom="tile",bins=78,data=subset(comp4,variable!="glyphosate_conc"&value==1),aes(y=variable,fill="no detectable reaction"),colour="black")+
	stat_bin_2d(geom="tile",bins=78,data=subset(comp4,variable=="species_richness_dna"&value!=1),aes(y=variable,fill="species_richness_dna"),colour="black")+
	stat_bin_2d(geom="tile",bins=78,data=subset(comp4,variable=="pielou.shannon"&value!=1),aes(y=variable,fill="pielou_shannon"),colour="black")+
	stat_bin_2d(geom="tile",bins=78,data=subset(comp4,variable=="gene_instance_remainder"&value!=1),aes(y=variable,fill="gene_instance_remainder"),colour="black")+
	stat_bin_2d(geom="tile",bins=78,data=subset(comp4,variable=="single_otu_plots"&value==7),aes(y=variable,fill="single_otu_plots_most"),colour="black")+
	stat_bin_2d(geom="tile",bins=78,data=subset(comp4,variable=="single_otu_plots"&value==8),aes(y=variable,fill="single_otu_plots_some"),colour="black")+
	#stat_bin_2d(geom="tile",bins=78,data=subset(comp4,variable=="glyphosate_conc"&value==9),aes(y=variable,fill="glyphosate_conc_high"),colour="black")+
	#stat_bin_2d(geom="tile",bins=78,data=subset(comp4,variable=="glyphosate_conc"&value==10),aes(y=variable,fill="glyphosate_conc_low"),colour="black")+
	geom_vline(aes(xintercept=0),colour="black",linetype="dotted", size=2)+
	#ggtitle("Duration of detected reactions to Glyphosate")+
	theme(plot.title = element_text(size=18,face="bold",hjust = 0.5))+
	theme(axis.text.y=element_blank(),
			axis.ticks.y=element_blank(),
			axis.title.y=element_blank(),
			panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
			panel.grid.minor=element_blank(),plot.background=element_blank())+
	xlab("Days")+
	coord_fixed(ratio = 1.5)+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 20))+
	#scale_y_continuous(breaks=c("significant_cell_counts","CCA_plot","no detectable\nreaction","species_richness_dna","pielou_shannon","gene_instance_remainder","single_otu_plots_most","single_otu_plots_some","glyphosate_conc_high","glyphosate_conc_low"),labels=letters[1:10])+
	theme(legend.position="none")+
	#geom_label(data=subset(comp4,new_days==-4),aes(y=variable),size=2.5)+
	scale_fill_manual(values=c("significant_cell_counts"="grey25","no detectable\nreaction"="white","species_richness_dna"="grey25","pielou_shannon"="grey25","gene_instance_remainder"="grey25","single_otu_plots_most"="grey25","single_otu_plots_some"="grey50"),
							name="Approaches",
							breaks=c("significant_cell_counts","no detectable\nreaction","species_richness_dna","pielou_shannon","gene_instance_remainder","single_otu_plots_most","single_otu_plots_some"))+
							annotate( "label", label = "Species richness\n(16S rRNA gene)", x = -7, y = 5, col = "black",size=3,hjust=0 ) +
							annotate( "label", label = "Single OTU", x = -7, y = 4, col = "black",size=3,hjust=0 )+
							annotate( "label", label = "Total cell counts", x = -7, y = 3, col = "black",size=3,hjust=0 )+
							annotate( "label", label = "Diversity indices", x = -7, y = 2, col = "black",size=3,hjust=0 )+
							annotate( "label", label = "Gene richness", x = -7, y = 1, col = "black",size=3,hjust=0 )
							#annotate( "label", label = "Changes in CCA", x = -7, y = 1, col = "black",size=3,hjust=0)
							#labels=c("Blooming","CCA","CCA_weak","no detectable\nreaction","DNA species richness","Diversity indices","memory effect","OTU reactions","OTU extended reactions","Glyphosate > 1 mg/L","Glyphosate < 1 mg/L"))
								ggsave(file="compare_approaches_paper_ohne_CCA.png", width=24, height=12)
								
#annotation!								
plot <- ggplot( diamonds ) + 
    geom_histogram( aes( carat ), bins = 30 ) +
    annotate( "text", label = "label here", x = 1, y = 750, col = "red" ) +
    annotate( "text", label = "and another", x = 2, y = 550, col = "blue" )

plot <- plot + 
    xlim( 0, 3 ) + 
    ggtitle( "Main title" ) +
    xlab( "label x" ) + 
    ylab( "label y" )

plot