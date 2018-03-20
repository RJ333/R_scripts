# creates a plot of averaged phn operon gene richness in relation to start 
# for treatment and control and the same for the species richness, but this without confidence intervall

library(ggplot2)
library(reshape2)

# omics_neu as workspace, read in phn_op_richness.csv, statistics by online available function summarySE
test <- read.csv(file.choose(), row.names = 1, sep = ",") 
test_summary <- summarySE(test, measurevar = "gene_richness_relative", groupvars = c("days" ,"treatment")) #creating statistical data based on script summarySE
fuck <- c("group_func", "group_func", "group_func", "group_func", "group_func", "group_func", "group_func", "group_func", "group_func", "group_func")
test_summary["group"] <- fuck 

# read in cell count data
more_cell_counts["new_days"] <- more_cell_counts$day - 69
more_cell_counts_0 <- subset(more_cell_counts, new_days >= 0)
write.csv(more_cell_counts_0, file = "more_cell_counts_0.csv")  # add column cells_rel
more_cell_counts_0 <- read.csv(file.choose(), row.names = 1, sep = ";") 

# read in taxonomy richness
tax_rich <- read.delim(file.choose(), header = TRUE)
tax_rich
fuck2 <- c(rep("group_tax", 32))
tax_rich["group"] <- fuck2 
tax_rich_omics <- subset(tax_rich, days >= 0)

test_mean_ci<-ggplot(test_summary, aes(x = days, y = gene_richness_relative, group = treatment, colour = treatment))+
geom_line(linetype = "dashed")+
geom_ribbon(aes(ymax = gene_richness_relative + ci, ymin = gene_richness_relative - ci), alpha = 0.1)+
geom_line(data = tax_rich_omics, aes(y = richness_rel))
test_mean_ci


#plot with averaged groups and confidence interval
my_y_title <- expression(paste(italic("phn"), " operon and species richness in relation to day 0 [%]"))
my_legend <- expression(paste(italic("phn"), " operon richness"))

phn_mean_ci <- ggplot(test_summary, aes(x = days, y = gene_richness_relative, group = treatment, colour = treatment))+
geom_line(aes(linetype = group), size = 2)+
geom_ribbon(aes(ymax = gene_richness_relative + ci, ymin = gene_richness_relative - ci), alpha = 0.1, colour = NA)+
geom_line(data = tax_rich_omics, aes(y = richness_rel, linetype = group), size = 1.5)+
#geom_point(data = more_cell_counts_0, aes(x = new_days, y = cells_rel, group = treatment), size = 3, alpha = 0.8)+
scale_colour_manual(values = c("glyph" = "black", "control" = "grey50"),
						name = "  ",
						breaks = c("glyph", "control"),
						labels = c("Treatment", "Control"))+
scale_linetype_manual(values = c("group_func" = 1, "group_tax" = 3),
						name = "",
						breaks = c("group_func", "group_tax"),
						labels = c(my_legend, "16S rRNA genes richness"))+								
theme_bw()+
	theme(panel.grid.major=element_line(colour = NA, size = 0.2))+
	theme(panel.grid.minor=element_line(colour = NA, size = 0.5))+
	scale_x_continuous(breaks = scales::pretty_breaks(n = 10))+
	theme(axis.title = element_text(size=20))+
	theme(axis.text.x = element_text(angle=0,vjust=0.5))+
	theme(axis.title.y = element_text(angle=90,vjust=0.5))+
	theme(axis.text=element_text(size=17))+
	theme(legend.position="right")+
	theme(strip.background = element_blank(),
		strip.text.x = element_blank())+
	xlab("Days after glyphosate addition")+
	ylab(my_y_title)+
	guides(lty = guide_legend(keywidth = 1.5, keyheight = 1))
phn_mean_ci
ggsave(file = "phn_and_tax.png", width = 14, height = 12)

my_y_title <- expression(paste(italic("phn"), " operon and species richness in relation to day 0 [%]"))
.... + labs(y=my_y_title)