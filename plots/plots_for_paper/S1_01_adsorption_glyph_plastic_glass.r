############# clean final version ################################################################

pd <- position_dodge(2)
ggplot(adsorp,aes(x=hours,y=rel_area,group=setup))+
	geom_line(aes(lty=setup),position=pd,size=1.5,alpha=1)+ 
	geom_errorbar(aes(ymin=rel_area-rel_SD, ymax=rel_area+rel_SD,lty=setup),position=pd, width=5, size=1.5)+
	scale_linetype_manual(values=c("liquid"= 6,"sediment"=3,"incubated"=1),
							name="Bottle setup with \nstandard deviation",
							breaks=c("liquid","sediment","incubated"),
							labels=c("ABW","ABW + quartz sediment","ABW + quartz sediment +\nbiota"))+
	theme_bw()+
	theme(panel.grid.major.x=element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor.x=element_line(colour = NA, size = 0.5))+
	scale_x_continuous(breaks=c(0,4,24,48,72))+
	theme(axis.title = element_text(size=20,face="bold"))+
	theme(axis.text=element_text(size=17,face="bold"))+
	theme(legend.position="none")+
	theme(axis.title.y = element_text(angle=90,vjust=0.5))+
	theme(strip.text.x=element_text(size=20,face="bold"))+
	xlab("Hours")+
	ylab("relative ratio  [%]")+
	facet_wrap(~material2,ncol=2,labeller=label_value)
ggsave(file="glyph_adsorption_paper.png", width=20, height=8.75)




#pd hat funktioniert nachdem x y und group in den obersten call verschoben wurden
pd <- position_dodge(2)
#adsorp<-read.csv(file.choose(),sep=";")
ggplot(adsorp,aes(x=hours,y=rel_area,group=setup))+
	geom_line(aes(lty=setup),position=pd,size=1,alpha=1)+ 
	geom_errorbar(aes(ymin=rel_area-rel_SD, ymax=rel_area+rel_SD,lty=setup),position=pd, width=5, size=1)+
	#geom_point(aes(y=rel_area,shape=setup),size=1.5)+
	scale_linetype_manual(values=c("liquid"= 6,"sediment"=3,"incubated"=1),
							name="Bottle setup with \nstandard deviation",
							breaks=c("liquid","sediment","incubated"),
							labels=c("LS","LS + HS","LS + HS +\nbiofilm"))+
	theme_bw()+
	theme(panel.grid.major.x=element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor.x=element_line(colour = NA, size = 0.5))+
	scale_x_continuous(breaks=c(0,4,24,48,72))+
	#ggtitle("Glyphosate adsorption in glass and\nplastic bottles with different setups")+
	theme(plot.title = element_text(size=18,face="bold",hjust = 0.5))+
	#theme(axis.title = element_text(size=12,face="bold"))+
	theme(axis.text=element_text(size=12))+
	#theme(legend.title=element_text(size=13,face="bold"))+
	#theme(legend.text=element_text(size=11))+
	theme(legend.position="bottom")+
	theme(legend.background = element_rect(fill="grey90",linetype="solid"))+
	theme(legend.text=element_text(size=11))+
	theme(axis.title.y = element_text(angle=90,vjust=0.5))+
	xlab("Hours")+
	ylab("relative ratio  [%]")+
	facet_wrap(~material2,ncol=2,labeller=label_value)
ggsave(file="glyph_adsorption_paper.png", width=14, height=8)

adsorp$material2<-factor(adsorp$material,labels=c("Glass bottle","Polypropylene bottle"))