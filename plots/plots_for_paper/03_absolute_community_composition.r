#creating grey scale palette
grey_ordered<-as.vector(scale_fill_grey(start=0.1,end=0.9))

#to get colours from a plot you can use
g <- ggplot_build(test_groesser_0.5)
grey_ordered<-unique(g$data[[1]]["fill"])

#shuffled colours
#To shuffle the order of the colors, you could use
sample_subset<-subset(abs_abu_molten_mean_tax,days>40&value>20000)
sample_subset<-droplevels(sample_subset)
require(scales)
n <- length(levels(sample_subset$order)) # number of colors
cols <- hue_pal(h = c(0, 360) + 15, 
                c = 100, l = 65, 
                h.start = 0, direction = 1)(n)[order(sample(1:n, n))] # color palette in random order
#or
n <- length(levels(sample_subset$order)) # number of colors
cols2 = rainbow(30, s=.7, v=1)[sample(1:26,26)]
col<-as.vector(c(cols,grey_ordered)
cols2<-col[sample(1:49,n)]
#ask to create extra legend on the side with cDNA and DNA
#ask for placing number of replicates below
	
sample_subset$treatment2<-factor(sample_subset$treatment,labels=c("Control","Treatment"))
more_cell_counts_44$treatment2<-factor(more_cell_counts_44$treatment,labels=c("Control","Treatment"))
cols2 = rainbow(30, s=.7, v=1)[sample(1:26,26)]
test_groesser_0.5<-ggplot(sample_subset, aes(x = new_days))+
	scale_fill_manual(values=cols2)+
	geom_bar(data=subset(sample_subset,nucleic_acid=="dna"&treatment=="glyph"),aes(x=new_days-0.7, y = value), width = 1.4, stat = "sum")+
	geom_bar(data=subset(sample_subset,nucleic_acid=="dna"&treatment=="glyph"),aes(x=new_days-0.7, y = value, fill=order), width = 1.3, stat = "identity")+
	geom_bar(data=subset(sample_subset,nucleic_acid=="cdna"&treatment=="glyph"),aes(x=new_days+0.7, y = value), width = 1.4, stat = "sum")+
	geom_bar(data=subset(sample_subset,nucleic_acid=="cdna"&treatment=="glyph"),aes(x=new_days+0.7, y = value, fill=order), width = 1.3, stat = "identity")+
	#geom_point(data=subset(more_cell_counts_44,treatment2=="Treatment"),aes(x=new_day,y=glyph_mg_L*5000000,group=treatment),size=4,alpha=1,shape=21,colour="black",fill="white")+
	geom_vline(data=subset(sample_subset,treatment=="glyph"),aes(xintercept=1.5),linetype="dashed", size=1.2)+
	#geom_text(data=sample_subset,aes(x=new_days,y=value,label=nucleic_acid, group=nucleic_acid), position=position_dodge(width=1), size=4)+
	guides(colour=FALSE,
		   fill=guide_legend(ncol=1,
						keyheight=1.5,
						label.theme=element_text(size=10,
												face="italic",
												angle=0),
						(title = NULL)))+
	theme_bw()+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 20))+
	theme(plot.title = element_text(size=18,face="bold",hjust=0.5))+
	theme(legend.background = element_rect(fill="grey90",linetype="solid"))+
	theme(panel.grid.major=element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor=element_line(colour = NA, size = 0.5))+
	scale_y_continuous(label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})+
	labs(x="Days",
			 y=expression(atop("Absolute abundance in cell equivalents",bgroup("[",relative~abundance~x~cells~mL^{-1},"]"))))
	#facet_wrap(~treatment2,nrow=1,ncol=2)
test_groesser_0.5
ggsave(file="absolute_community_order_level_for_manuscript.png", width=14, height=8)				


#missing ARKICE90 order
levels(sample_subset$order)[match("",sample_subset$order)]<- "ARKICE-90"