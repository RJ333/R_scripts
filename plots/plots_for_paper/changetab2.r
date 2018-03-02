change_tab2
cols = c("#FEECDE", "#FDD0A1", "#FDAD6B", "#FD8C3B", "#E5540C")
change_tab2_cdna<-subset(change_tab2, nucleic_acid=="cdna")	
change_tab2_cdna$otu <- factor(change_tab2_cdna$otu, levels=unique(change_tab2_cdna[order(change_tab2_cdna$change,decreasing=F),]$otu))
ggplot(change_tab2_cdna, aes(y=otu))+
	geom_tile(aes(x=nucleic_acid,fill = change))+
	geom_point(aes(x="",size=abu),shape=16,colour="black")+
	scale_fill_manual(values=cols,labels=c("7 - 10","10 - 20","20 - 50 ","50 - 100","100 - 150"),name="range of\nchange")+
	scale_size_discrete(labels=c("0.005-0.01","0.01-0.1","0.1-0.5","0.5-3","3-15"),name="maximum relative\nabundance",range = c(2, 8))+
	theme_classic()+
	xlab("")+
	ylab("OTU")+
	ggtitle("Changes in OTU abundance")+
	coord_fixed(ratio = 0.4)+
	scale_x_discrete(breaks=c("cdna",""),labels=c("cdna\nchange","max_rel\n_abu"))+
	theme(plot.title = element_text(size=18,face="bold"))+
	theme(axis.title = element_text(size=12,face="bold"))+
	theme(axis.title.y = element_text(angle=0,vjust=0.5))+
	theme(axis.text=element_text(size=12))+
	theme(legend.title=element_text(size=13,face="bold"))+
	theme(legend.text=element_text(size=11))+
	theme(axis.text.y = element_text(face = "italic"))+
	#theme(axis.text.x = element_blank())+
	theme(axis.ticks = element_blank())