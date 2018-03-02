pd <- position_dodge(1)
more_cell_counts<-read.csv(file.choose(),sep=";")
#######einfacher plot für zellzahlen und glyphosatkonzentration
more_cell_counts_62<-subset(more_cell_counts,day>61)
more_cell_counts_66<-subset(more_cell_counts,day>65)
#for glyphosate:
more_cell_counts_62_glyph<-subset(more_cell_counts,day>61&treatment=="glyph")
ggplot(more_cell_counts_62_glyph,aes(x=day))+
	geom_point(aes(y=glyph_theor,colour="glyphosate\ndilution\n"),alpha=1, size=1.5)+
	geom_point(aes(y=glyph_mg_L,colour="glyphosate\nconcentration\n"),alpha=1, size=1.5)+
	geom_errorbar(aes(ymin=glyph_mg_L-glyph_se, ymax=glyph_mg_L+glyph_se), width=.2)+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 20))+
	ggtitle("glyphosate")+
	theme(plot.title = element_text(size=18,face="bold",hjust = 0.5))+
	theme(axis.title = element_text(size=12,face="bold"))+
	#theme(axis.text.y = element_text(angle=90))+
	theme(axis.text=element_text(size=12))+
	theme(legend.title=element_text(size=13,face="bold"))+
	theme(legend.text=element_text(size=11))+
	theme(axis.title.y = element_text(angle=90,vjust=0.5))+
	xlab("day")+
	ylab("glyphosate mg/L")+
	facet_wrap( ~treatment,nrow=2)

 ggsave(file="cells_counts_glyphR.jpg", width=14, height=8)
 
 ggplot(more_cell_counts_62,aes(x=day))+
	geom_errorbar(aes(ymin=cells_ml-cells_se, ymax=cells_ml+cells_se), width=.2)+
	geom_point(data=subset(more_cell_counts_62,group=="normal"),aes(y=cells_ml,colour="cells/mL\n"),size=3,alpha=0.8,color="black")+
	geom_point(data=subset(more_cell_counts_62,group=="bloom"),aes(y=cells_ml,colour="bloom"),size=3.5,alpha=0.8,color="blue")+
	geom_point(data=subset(more_cell_counts_62,group=="second"),aes(y=cells_ml,colour="bloom2"),size=3.5,alpha=0.8,color="green")+
	geom_vline(data=subset(more_cell_counts_62, treatment =="glyph"), aes(xintercept=70),linetype="dashed", size=1.2,color="blue")+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 20))+
	ggtitle("Cell counts in water columns")+
	theme(plot.title = element_text(size=18,face="bold",hjust = 0.5))+
	theme(axis.title = element_text(size=12,face="bold"))+
	#theme(axis.text.y = element_text(angle=90))+
	theme(axis.text=element_text(size=12))+
	theme(legend.title=element_text(size=13,face="bold"))+
	theme(legend.text=element_text(size=11))+
	theme(axis.title.y = element_text(angle=90,vjust=0.5))+
	xlab("day")+
	ylab("cells/mL")+
	facet_wrap( ~treatment,nrow=2)
ggsave(file="cellcounts_manuscript.jpg", width=10, height=8)
	

adsorp<-read.csv(file.choose(),sep=";")
ggplot(adsorp,aes(x=hours))+
	geom_line(aes(y=area_ratio,colour=setup),size=1.5,alpha=0.8,position=pd)+ #farbe stimmt mit legende noch nich überein
	geom_errorbar(aes(ymin=area_ratio-SD, ymax=area_ratio+SD),position=pd, width=1)+
	geom_point(aes(y=area_ratio),size=2,alpha=0.8,position=pd)+
	ggtitle("glyphosate adsorption")+
	theme(plot.title = element_text(size=18,face="bold",hjust = 0.5))+
	theme(axis.title = element_text(size=12,face="bold"))+
	theme(axis.text=element_text(size=12))+
	theme(legend.title=element_text(size=13,face="bold"))+
	theme(legend.text=element_text(size=11))+
	theme(axis.title.y = element_text(angle=90,vjust=0.5))+
	xlab("hour")+
	ylab("AreaGly/AreaIS")
ggsave(file="glyph_adsorption.jpg", width=10, height=8)