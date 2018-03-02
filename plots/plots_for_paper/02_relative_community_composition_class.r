require(ggplot2)
require(scales)
############# clean final version ################################################################
#data preparation, setting threshold for samples and clusters which are included in plot
sample_subset<-subset(rel_abu_molten_mean_tax,days>40&value>0.15&habitat=="water"&treatment=="glyph") #relative abundance above 0.15 %
sample_subset<-droplevels(sample_subset)
sample_subset$new_days<-sample_subset$days-69 #set day of glyphosate addition to 0

length(levels(sample_subset$class)) # number of colors needed for all classes
length(levels(sample_subset$order)) # number of colors needed for all orders

levels(sample_subset$order)[match("",sample_subset$order)]<- "ARKICE-90" #add missing value
#sort orders based on class, first alpha, then beta, then gammaproteos
sample_subset$order<-factor(sample_subset$order, levels=c("Caulobacterales","Rhizobiales","Rhodobacterales","Rhodospirillales","Sphingomonadales",
"Burkholderiales","Hot Creek 32","Methylophilales","Alteromonadales","Gammaproteobacteria Incertae Sedis","Pseudomonadales","Sphingobacteriales","ARKICE-90","Flavobacteriales"))
#define fill colours for all orders
fill_values<-c("Flavobacteriales"="green",
	"Rhizobiales"="red",
	"Pseudomonadales"="grey30",
	"Rhodobacterales"="black",
	"Caulobacterales"="lightblue",
	"Gammaproteobacteria Incertae Sedis"="white",
	"Burkholderiales"="pink",
	"Sphingobacteriales"="grey75",
	"Hot Creek 32"="grey",
	"Alteromonadales"="orange",
	"Sphingomonadales"="purple",
	"Rhodospirillales"="yellow",
	"ARKICE-90"="blue2",
	"Methylophilales"="green3")

#plotting all selected clusters in bar plot ordered by class and displaying orders over time for DNA and RNA
test_groesser_0.5_class<-ggplot(sample_subset, aes(x = new_days,group=order))+
	scale_fill_manual(breaks=levels(sample_subset$order),values=fill_values)+
	geom_bar(data=subset(sample_subset,nucleic_acid=="dna"&treatment=="glyph"),aes(x=new_days-0.5, y = value), fill="black", width = 0.9, stat = "sum")+
	geom_bar(data=subset(sample_subset,nucleic_acid=="dna"&treatment=="glyph"),aes(x=new_days-0.5, y=value, fill=order), width = .6, stat = "identity")+
	geom_bar(data=subset(sample_subset,nucleic_acid=="cdna"&treatment=="glyph"),aes(x=new_days+0.5, y = value), fill="black", width = 0.9, stat = "sum")+
	geom_bar(data=subset(sample_subset,nucleic_acid=="cdna"&treatment=="glyph"),aes(x=new_days+0.5, y = value, fill=order), width = 0.6, stat = "identity")+
	geom_vline(data=subset(sample_subset,treatment=="glyph"),aes(xintercept=1.5),linetype="dashed", size=1.2)+
	guides(colour=FALSE, size=FALSE, width=FALSE,
		   fill=guide_legend(ncol=1,
						keyheight=1.5,
						label.theme=element_text(size=15,
												face="italic",
												angle=0),
						(title = NULL)))+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
	scale_y_continuous(expand = c(0,0))+
	theme_bw()+
	theme(axis.text=element_text(size=17))+
	theme(axis.title=element_text(size=20,face="bold"))+
	theme(legend.background = element_rect(fill="grey90",linetype="solid"))+
	theme(panel.grid.major=element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor=element_line(colour = NA, size = 0.5))+
	labs(x="Days", y="Relative abundance [%]")+
  annotate("text", x = -27, y = 90, label = "a)" , color="black", size=6 , angle=0, fontface="bold")+
  annotate("text", x = -22.5, y = 90, label = "b)" , color="black", size=6 , angle=0, fontface="bold")
test_groesser_0.5_class
ggsave(file="fig_02_relative_community_order_level_sorted_for_manuscript.png", width=16, height=8)

#############################################################dev version


sample_subset<-subset(rel_abu_molten_mean_tax,days>40&value>0.15&habitat=="water"&treatment=="glyph")
sample_subset<-droplevels(sample_subset)
sample_subset$new_days<-sample_subset$days-69
length(levels(sample_subset$class)) # number of colors
length(levels(sample_subset$order)) # number of colors

require(scales)
n <- length(levels(sample_subset$class)) # number of colors
n <- length(levels(sample_subset$order)) # number of colors

	
sample_subset$treatment2<-factor(sample_subset$treatment,labels=c("Control","Treatment"))
more_cell_counts_44$treatment2<-factor(more_cell_counts_44$treatment,labels=c("Control","Treatment"))


test_groesser_0.5_class<-ggplot(sample_subset, aes(x = new_days,group=class))+
	#scale_fill_manual(values=c("red","green","orange","blue","black","pink","yellow","grey","green2","white","beige","grey50"))+
	#geom_bar(data=subset(sample_subset,nucleic_acid=="dna"&treatment=="glyph"),aes(x=new_days-1.25, y = value, fill=class), width = 0.5, stat = "sum")+
	geom_bar(data=subset(sample_subset,nucleic_acid=="dna"&treatment=="glyph"),aes(x=new_days-0.5, y = value, fill=order), width = 1., stat = "identity")+
	#geom_bar(data=subset(sample_subset,nucleic_acid=="cdna"&treatment=="glyph"),aes(x=new_days+1.25, y = value, fill=class), width = 0.5, stat = "sum")+
	geom_bar(data=subset(sample_subset,nucleic_acid=="cdna"&treatment=="glyph"),aes(x=new_days+0.5, y = value, fill=order), width = 1., stat = "identity")+
	#geom_point(data=subset(more_cell_counts_44,treatment2=="Treatment"),aes(x=new_day,y=glyph_mg_L,group=treatment),size=4,alpha=1,shape=21,colour="black",fill="white")+
	geom_vline(data=subset(sample_subset,treatment=="glyph"),aes(xintercept=1.5),linetype="dashed", size=1.2)+
	#geom_text(data=sample_subset,aes(x=new_days,y=value,label=nucleic_acid, group=nucleic_acid), position=position_dodge(width=1), size=4)+
	guides(colour=FALSE, size=FALSE, width=FALSE,
		   fill=guide_legend(ncol=1,
						keyheight=1.5,
						label.theme=element_text(size=10,
												face="italic",
												angle=0),
						(title = NULL)))+
	theme_bw()+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 7))+
	theme(plot.title = element_text(size=18,face="bold",hjust=0.5))+
	theme(legend.background = element_rect(fill="grey90",linetype="solid"))+
	theme(panel.grid.major=element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor=element_line(colour = NA, size = 0.5))+
	#scale_y_continuous(label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))})+
	labs(x="Days",
			 y=expression(atop("relative abundance")))+
  annotate("text", x = -27, y = 90, label = "a)" , color="black", size=6 , angle=0, fontface="bold")+
  annotate("text", x = -23, y = 90, label = "b)" , color="black", size=6 , angle=0, fontface="bold")
  #geom_text(data=subset(sample_subset,nucleic_acid=="dna"&treatment=="glyph"&new_days==0),aes(y=value,label=order),position="stack")+
	#facet_wrap(~treatment2,nrow=1,ncol=2)
test_groesser_0.5_class
ggsave(file="relative_community_class_level_for_manuscript.png", width=14, height=8)				


#missing ARKICE90 order
levels(sample_subset$order)[match("",sample_subset$order)]<- "ARKICE-90"

#relevel factor levels, takes the "ref" and puts it first
sample_subset$class<-relevel(sample_subset$class,ref="Gammaproteobacteria")
sample_subset$class<-relevel(sample_subset$class,ref="Betaproteobacteria")
sample_subset$class<-relevel(sample_subset$class,ref="Alphaproteobacteria")
levels(sample_subset$class)

#sort by class

sample_subset<-sample_subset[order(sample_subset$class),]

#for plotting: "group=class" important