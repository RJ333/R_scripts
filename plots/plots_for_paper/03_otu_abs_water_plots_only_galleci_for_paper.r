#write.csv(abs_abu_molten_tax,file="abs_abu_molten_tax.csv")
###########otu plots with glyphosate addition
############# clean final version ################################################################
require(scales)
require(gridExtra)
require(grid)
require(ggplot2)
########## first the Gallaecimonas plot
Gallaecimonas_abs_water_plot<-subset(abs_abu_molten_tax, habitat == "water" & grepl('Gallaecimonas',genus) & days > 40)
species_title_Gallaecimonas<-expression(paste(,italic("Gallaecimonas")," sp."))
Gallaecimonas_abs_water_plot$treatment2<-factor(Gallaecimonas_abs_water_plot$treatment,labels=c("Control","Treatment"))

a<-ggplot(data=Gallaecimonas_abs_water_plot, aes(x=days-69,y=value,group=nucleic_acid,lty=nucleic_acid))+ 
	#coord_cartesian(ylim = c(0, 1))+
	annotate("segment", x = -8.0, xend = 18, y = 4300000, yend = 4300000, colour = "grey50", size=2, alpha=1,lineend = "butt")+
	annotate("segment", x = 18.0, xend = 18, y = 4440000, yend = 4160000, colour = "grey50", size=2, alpha=1,lineend = "butt")+
	annotate("segment", x = -8.0, xend = -8, y = 4440000, yend = 4160000, colour = "grey50", size=2, alpha=1,lineend = "butt")+
	geom_vline(aes(xintercept=0),linetype="dashed", size=1.2)+
	geom_point(data=subset(Gallaecimonas_abs_water_plot,treatment=="control"),aes(colour=treatment),alpha=1)+
	stat_summary(data=subset(Gallaecimonas_abs_water_plot,treatment=="control"),aes(colour=treatment),fun.y="mean", geom="line", size=2,alpha=1)+
	stat_summary(data=subset(Gallaecimonas_abs_water_plot,treatment=="glyph"),aes(colour=treatment),fun.y="mean", geom="line", size=2)+
	geom_point(data=subset(Gallaecimonas_abs_water_plot,treatment=="glyph"),aes(colour=treatment))+
	scale_linetype_manual(values=c("dna"=1,"cdna"=6),
						name="Nucleic acid  ",
						breaks=c("cdna","dna"),
						labels=c("16S rRNA","16S rRNA gene"))+
	scale_colour_manual(values=c("glyph"="black","control"="grey50"),
						name="Microcosm  ",
						breaks=c("glyph","control"),
						labels=c("Treatment","Control"))+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
	theme_bw()+
	scale_y_continuous(label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})+
	theme(axis.text=element_text(size=18))+
	theme(panel.grid.major=element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor=element_line(colour = NA, size = 0.5))+
	theme(legend.position="none")+
	theme(axis.title=element_blank())
	
	
#########################################################################then combine with gridExtra	
c<-expression(bold(atop("Absolute abundance in cell equivalents",bgroup("[",relative~abundance~x~cells~mL^{-1},"]"))))

grid.arrange(a, nrow = 1, bottom=textGrob(expression(bold("Days")),gp=gpar(fontsize=22)),left=textGrob(c,rot=90,gp=gpar(fontsize=22)))

ggsave(file="fig_03_OTU_plots_paper_Gallaeci.png", width=12, height=8.75,
	grid.arrange(a, nrow = 1, 
		bottom=textGrob(expression(bold("Days")),gp=gpar(fontsize=22)),
		left=textGrob(c,rot=90,gp=gpar(fontsize=22))))








######################################################################################DEV code

##########Gallaecimonas
Gallaecimonas_abs_water_plot<-subset(abs_abu_molten_tax, habitat == "water" & grepl('Gallaecimonas',genus) & days > 40)
species_title_Gallaecimonas<-expression(paste(,italic("Gallaecimonas")," sp."))
Gallaecimonas_abs_water_plot$treatment2<-factor(Gallaecimonas_abs_water_plot$treatment,labels=c("Control","Treatment"))
a<-ggplot(data=Gallaecimonas_abs_water_plot, aes(x=days-69,y=value,group=nucleic_acid,lty=nucleic_acid))+ 
	#coord_cartesian(ylim = c(0, 1))+
		annotate("segment", x = -8.0, xend = 18, y = 3500000, yend = 3500000, colour = "grey80", size=2, alpha=1,lineend = "butt")+
	annotate("segment", x = 18.0, xend = 18, y = 3650000, yend = 3350000, colour = "grey80", size=2, alpha=1,lineend = "butt")+
	annotate("segment", x = -8.0, xend = -8, y = 3650000, yend = 3350000, colour = "grey80", size=2, alpha=1,lineend = "butt")+
	geom_vline(aes(xintercept=0),linetype="dashed", size=1.2)+
	geom_point(data=subset(Gallaecimonas_abs_water_plot,treatment=="glyph"),aes())+
	geom_point(data=subset(Gallaecimonas_abs_water_plot,treatment=="control"),aes())+
	stat_summary(data=subset(Gallaecimonas_abs_water_plot,treatment=="glyph"),aes(colour=treatment),fun.y="mean", geom="line", size=1.5)+
	stat_summary(data=subset(Gallaecimonas_abs_water_plot,treatment=="control"),aes(colour=treatment),fun.y="mean", geom="line", size=1.5)+
	scale_linetype_manual(values=c("dna"=1,"cdna"=6),
						name="Nucleic acid  ",
						breaks=c("cdna","dna"),
						labels=c("16S rRNA","16S rRNA gene"))+
	scale_colour_manual(values=c("glyph"="black","control"="grey60"),
						name="Microcosm  ",
						breaks=c("glyph","control"),
						labels=c("Treatment","Control"))+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 20))+
	#ggtitle(species_title_Gallaecimonas)+
	theme_bw()+
	scale_y_continuous(label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})+
	#theme(plot.title = element_text(size=18,face="bold",hjust=0.5))+
	#theme(axis.title = element_text(size=12,face="bold"))+
	theme(axis.title.y = element_text(angle=90,vjust=0.5))+
	theme(axis.text=element_text(size=12))+
	#theme(legend.title=element_text(size=13,face="bold"))+
	theme(legend.text=element_text(size=11))+
	theme(panel.grid.major=element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor=element_line(colour = NA, size = 0.5))+
	theme(legend.position="bottom")+
	theme(legend.background = element_rect(fill="grey90",linetype="solid"))+
	labs(x="Days",
			 y=expression(atop("Absolute abundance in cell equivalents",bgroup("[",relative~abundance~x~cells~mL^{-1},"]"))))

ggsave(file="absolute_genus_Gallaecimonas_water_for_paper_oneplot_bar.png", width=10, height=12)

#Escherichia
Escherichia_abs_water_plot<-subset(abs_abu_molten_tax, habitat == "water" & grepl('Escherichia',genus) & days > 40)
species_title_Escherichia<-expression(paste(,italic("Escherichia")," sp."))
Escherichia_abs_water_plot$treatment2<-factor(Escherichia_abs_water_plot$treatment,labels=c("Control","Treatment"))
b<-ggplot(data=Escherichia_abs_water_plot, aes(x=days-69,y=value,group=nucleic_acid,lty=nucleic_acid))+ 
	#coord_cartesian(ylim = c(0, 1))+
		annotate("segment", x = -8.0, xend = 18, y = 80000, yend = 80000, colour = "grey80", size=2, alpha=1,lineend = "butt")+
	annotate("segment", x = 18.0, xend = 18, y = 82000, yend = 78000, colour = "grey80", size=2, alpha=1,lineend = "butt")+
	annotate("segment", x = -8.0, xend = -8, y = 82000, yend = 78000, colour = "grey80", size=2, alpha=1,lineend = "butt")+
	geom_vline(aes(xintercept=0),linetype="dashed", size=1.2)+
	geom_point(data=subset(Escherichia_abs_water_plot,treatment=="glyph"),aes())+
	geom_point(data=subset(Escherichia_abs_water_plot,treatment=="control"),aes())+
	stat_summary(data=subset(Escherichia_abs_water_plot,treatment=="glyph"),aes(colour=treatment),fun.y="mean", geom="line", size=1.5)+
	stat_summary(data=subset(Escherichia_abs_water_plot,treatment=="control"),aes(colour=treatment),fun.y="mean", geom="line", size=1.5)+
	scale_linetype_manual(values=c("dna"=1,"cdna"=6),
						name="Nucleic acid  ",
						breaks=c("cdna","dna"),
						labels=c("16S rRNA","16S rRNA gene"))+
	scale_colour_manual(values=c("glyph"="black","control"="grey60"),
						name="Microcosm  ",
						breaks=c("glyph","control"),
						labels=c("Treatment","Control"))+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 20))+
	#ggtitle(species_title_Escherichia)+
	theme_bw()+
	scale_y_continuous(label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})+
	#theme(plot.title = element_text(size=18,face="bold",hjust=0.5))+
	#theme(axis.title = element_text(size=12,face="bold"))+
	theme(axis.title.y = element_text(angle=90,vjust=0.5))+
	theme(axis.text=element_text(size=12))+
	#theme(legend.title=element_text(size=13,face="bold"))+
	theme(legend.text=element_text(size=11))+
	theme(panel.grid.major=element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor=element_line(colour = NA, size = 0.5))+
	theme(legend.position="bottom")+
	theme(legend.background = element_rect(fill="grey90",linetype="solid"))+
	labs(x="Days",
			 y=expression(atop("Absolute abundance in cell equivalents",bgroup("[",relative~abundance~x~cells~mL^{-1},"]"))))
ggsave(file="absolute_genus_Escherichia_water_for_paper_oneplot_bar.png", width=10, height=12)

ggsave(file="OTU_plots_paper_Ecoli_Gallaeci.png", width=20, height=8.75,arrangeGrob(a,b,nrow=1))
###Benennung der beiden OTUs in plot sowie Skalenanpassung je OTU fehlt noch

#Gall_Esch
Gall_Esch_abs_water_plot<-subset(abs_abu_molten_tax, habitat == "water" & grepl('Gallaecimonas|Escherichia',genus) & days > 40)
species_title_Gall_Esch<-expression(paste(,italic("Gall_Esch")," sp."))
Gall_Esch_abs_water_plot$treatment2<-factor(Gall_Esch_abs_water_plot$treatment,labels=c("Control","Treatment"))
ggplot(data=Gall_Esch_abs_water_plot, aes(x=days-69,y=value,colour=nucleic_acid,group=nucleic_acid,lty=nucleic_acid))+ 
	#coord_cartesian(ylim = c(0, 1))+
	geom_vline(data=subset(Gall_Esch_abs_water_plot,treatment=="glyph"),aes(xintercept=0),linetype="dashed", size=1.2)+
	geom_point(aes(),colour="black")+
	stat_summary(aes(colour=nucleic_acid),colour="black",fun.y="mean", geom="line", size=1.5)+
	facet_wrap(~treatment2*genus,nrow=2,ncol=2)+
	scale_linetype_manual(values=c("dna"=1,"cdna"=4),
						name="Nucleic acid  ",
						breaks=c("cdna","dna"),
						labels=c("16S rRNA","16S rDNA"))+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 20))+
	#ggtitle(species_title_Gall_Esch)+
	theme_bw()+
	scale_y_continuous(label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})+
	#theme(plot.title = element_text(size=18,face="bold",hjust=0.5))+
	#theme(axis.title = element_text(size=12,face="bold"))+
	theme(axis.title.y = element_text(angle=90,vjust=0.5))+
	theme(axis.text=element_text(size=12))+
	#theme(legend.title=element_text(size=13,face="bold"))+
	theme(legend.text=element_text(size=11))+
	theme(panel.grid.major=element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor=element_line(colour = NA, size = 0.5))+
	theme(legend.position="bottom")+
	theme(legend.background = element_rect(fill="grey90",linetype="solid"))+
	labs(x="Days",
			 y=expression(atop("Absolute abundance in cell equivalents",bgroup("[",relative~abundance~x~cells~mL^{-1},"]"))))
ggsave(file="absolute_genus_Gall_Esch_water_for_paper.png", width=14, height=8)
