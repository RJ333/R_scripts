# creates a plot of averaged phn operon gene richness in relation to start 
# for treatment and control and the same for the species richness, but this without confidence intervall

library(ggplot2)
library(reshape2)

# omics_neu as workspace, read in phn_op_richness.csv, statistics by online available function summarySE
not_phn <- read.csv(file.choose(), row.names = 1, sep = ",") 

# ydiF, malF, ktrB from "phnM_like" too small
# norB from "nitrogen" too small
# pphA, phnR, phnA, phnT from "more_phosphonate"
# chtA not found at all, chiA always 2
# all glycosyl hydrolases too snall
# tfdB fdm celB cenC in some samples not present

# remove to small genes, replace NA with 0
not_phn[is.na(not_phn)] <- 0
not_phn <- subset(not_phn, !grepl("fdm|ydiF|norB|malF|ktrB|pphA|phnR|phnA|phnT|chtA|chiA|xynB|bcsZ|celB|cenC", gene))

# creating statistical data based on script summarySE, it does not understand "NA"! dataframe[is.na(dataframe)] <- 0
not_phn_summary <- summarySE(not_phn, measurevar = "gene_richness_relative", groupvars = c("days" ,"treatment")) 
fuck <- c("group_func", "group_func", "group_func", "group_func", "group_func", "group_func", "group_func", "group_func", "group_func", "group_func")
not_phn_summary["group"] <- fuck 

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


#plot with averaged groups and confidence interval
my_y_title <- expression(paste("not ",italic("phn"), " operon and species richness in relation to day 0 [%]"))
my_legend <- expression(paste("not ",italic("phn"), " operon richness"))

not_phn_mean_ci <- ggplot(not_phn_summary, aes(x = days, y = gene_richness_relative, group = treatment, colour = treatment))+
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
not_phn_mean_ci
ggsave(file = "not_phn_and_tax.png", width = 14, height = 12)

my_y_title <- expression(paste(italic("phn"), " operon and species richness in relation to day 0 [%]"))
.... + labs(y=my_y_title)