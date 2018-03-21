pd <- position_dodge(1)
more_cell_counts <- read.csv(file.choose(), sep = ";")
#######einfacher plot für zellzahlen und glyphosatkonzentration
more_cell_counts_62 <- subset(more_cell_counts, day > 61)
more_cell_counts_66 <- subset(more_cell_counts, day > 65)
more_cell_counts_0 <- subset(more_cell_counts, new_day >= 0)

	


############################	
library(ggplot2)
library(scales)

ggplot(more_cell_counts_0,aes(x=new_day,colour=treatment,linetype=treatment,group=treatment))+
	geom_errorbar(aes(ymin=cells_ml-cells_se, ymax=cells_ml+cells_se), colour="black",width=2,size=1.2)+
	geom_line(aes(y=cells_ml,group=treatment,linetype=treatment),colour="black",size=2,alpha=0.8)+
	#geom_point(aes(y=glyph_theor*2000000,shape="glyph"),colour="black",alpha=1, size=4)+
	#geom_point(aes(y=glyph_mg_L*2000000,shape="glyph_deg"),colour="black",alpha=1, size=4)+
	#geom_errorbar(aes(ymin=(glyph_mg_L-glyph_se)*2000000, ymax=(glyph_mg_L+glyph_se)*2000000),colour="black", width=2,size=1.2)+
	#geom_vline(data=subset(more_cell_counts_0, treatment =="glyph"), aes(xintercept=0),linetype="dashed", size=1.5)+
		theme_bw()+
		#scale_y_continuous(label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))}, sec.axis=sec_axis(~./2000000,name=expression(bold(paste("Glyphosate concentration  ",bgroup("[",mg~L^{-1},"]"))))))+
		scale_linetype_manual(values=c("control"=4,"glyph"=1),
			name="Microcosm",
			breaks=c("control","glyph"),
			labels=c("Control","Treat-\nment"))+
		scale_shape_manual(values=c("glyph"=2,"glyph_deg"=17),
			name="Glyphosate decrease by",
			breaks=c("glyph","glyph_deg"),
			labels=c("calculated\ndilution","measured\nconcentration"))+
		scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
		theme(panel.grid.major=element_line(colour = NA, size = 0.2))+
		theme(panel.grid.minor=element_line(colour = NA, size = 0.5))+
		theme(axis.title = element_text(size=20,face="bold"))+
	theme(axis.text=element_text(size=17,face="bold"))+
	theme(legend.position="none")+
	theme(axis.title.y = element_text(angle=90,vjust=0.5))+
		labs(x="Days",
			 y=expression(bold(paste("Total cell count  ",bgroup("[",cells~mL^{-1},"]")))))
			 
ggsave(file="cellcounts_only_after.png", width=14, height=8)