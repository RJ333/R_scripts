adsorp<-read.csv(file.choose(),sep=";")
ggplot(adsorp,aes(x=hours))+
	geom_line(aes(y=rel_area,group=setup),size=0.5,alpha=0.6)+ 
	geom_errorbar(aes(ymin=rel_area-rel_SD, ymax=rel_area+rel_SD,colour=setup),position=pd, width=3)+
	#geom_point(aes(y=rel_area,shape=setup),size=1.5)+
	#scale_shape_manual(values=c("liquid"=15,"sediment"=16,"incubated"=17),
							#name="setup",
							#breaks=c("liquid","sediment","incubated"),
							#labels=c("LS\n","LS + HS II\n","LS + HS II +\nbiofilm\n"))+
	scale_colour_manual(values=c("liquid"="darkgreen","sediment"="red","incubated"="purple"),
							name="setup with \nstandard error",
							breaks=c("liquid","sediment","incubated"),
							labels=c("LS\n","LS + HS II\n","LS + HS II +\nbiofilm\n"))+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 9))+
	ggtitle("glyphosate adsorption")+
	theme(plot.title = element_text(size=18,face="bold",hjust = 0.5))+
	theme(axis.title = element_text(size=12,face="bold"))+
	theme(axis.text=element_text(size=12))+
	theme(legend.title=element_text(size=13,face="bold"))+
	theme(legend.text=element_text(size=11))+
	theme(axis.title.y = element_text(angle=90,vjust=0.5))+
	xlab("hours")+
	ylab("relative ratio")
ggsave(file="glyph_adsorption_paper.jpg", width=9, height=8)