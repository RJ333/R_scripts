##VennDiagram
library(data.table)
library(VennDiagram)
library(ggplot2)
library(reshape2)

outersect <- function(x, y, ...) {
  big.vec <- c(x, y, ...)
  duplicates <- big.vec[duplicated(big.vec)]
  setdiff(big.vec, unique(duplicates))
}
threeway.Venn <- function(A,B,C,cat.names = c("A", "B", "C")){
  grid.newpage()
  area1 <- length(A)
  area2 <- length(B)
  area3 <- length(C)
  n12 <- length(intersect(A,B))
  n23 <- length(intersect(B,C))
  n13 <- length(intersect(A,C))
  n123 <- length(intersect(intersect(A, B), intersect(B,C )))
venn.plot <- draw.triple.venn(
  area1 = area1,
  area2 = area2,
  area3 = area3,
  n12 = n12,
  n23 = n23,
  n13 = n13,
  n123 = n123,
  category = cat.names,
  cat.pos = c(0,180,0),
  fill = c("blue", "red", "green"),
  alpha = .3,
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.col = c("blue", "red", "green")
)
grid.draw(venn.plot)
}
fourway.Venn <- function(A,B,C,D,cat.names = c("water treatment", "water control", "biofilm treatment","biofilm control")){
  grid.newpage()
  area1 <- length(A)
  area2 <- length(B)
  area3 <- length(C)
  area4 <- length(D)
  n12<-length(Reduce(intersect, list(A,B)))
  n13<-length(Reduce(intersect, list(A,C)))
  n14<-length(Reduce(intersect, list(A,D)))
  n23<-length(Reduce(intersect, list(B,C)))
  n24<-length(Reduce(intersect, list(B,D)))
  n34<-length(Reduce(intersect, list(C,D)))
  n123<-length(Reduce(intersect, list(A,B,C)))
  n124<-length(Reduce(intersect, list(A,B,D)))
  n134<-length(Reduce(intersect, list(A,C,D)))
  n234<-length(Reduce(intersect, list(B,C,D)))
  n1234<-length(Reduce(intersect, list(A,B,C,D)))
  
venn.plot <- draw.quad.venn(
  area1 = area1,
  area2 = area2,
  area3 = area3,
  area4 = area4,
  n12 = n12,
  n13 = n13,
  n14 = n14,
  n23 = n23,
  n24 = n24,
  n34 = n34,
  n123 = n123,
  n124 = n124,
  n134 = n134,
  n234 = n234,
  n1234 = n1234,
  category = cat.names,
  cat.pos = c(0,180,0,270),
  fill = c("blue", "red", "green","yellow"),
  alpha = .3,
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.col = c("blue", "red", "green","yellow")
)
grid.draw(venn.plot)
}


#multiple intersects:
Reduce(intersect, list(a,b,c,d))
##VennDiagram needs to know how many entries each group and the overlapping sections contain
final_table_tax<-read.csv(file.choose(),row.names=1,sep=";")
final_table_wholetax<-final_table_tax[,c(1:8)]
final_table_wholetax_mean<-aggregate(value~variable+time+days+treatment+habitat+nucleic_acid, data = final_table_wholetax, mean)


#######################################################################
#####subsetting 2 groups:
water_otus<-subset(final_table_wholetax_mean,habitat=="water"&value!=0)
water_unique_otus<-water_otus[which(!duplicated(water_otus[,"variable"])),]
nrow(water_unique_otus)
#536

biofilm_otus<-subset(final_table_wholetax_mean,habitat=="biofilm"&value!=0)
biofilm_unique_otus<-biofilm_otus[which(!duplicated(biofilm_otus[,"variable"])),]
nrow(biofilm_unique_otus)
#494

length(intersect(water_unique_otus$variable,biofilm_unique_otus$variable))
#377
####################################################################

water_dna_otus<-subset(final_table_wholetax_mean,habitat=="water"&nucleic_acid=="dna"&value!=0)
water_dna_unique_otus<-water_dna_otus[which(!duplicated(water_dna_otus[,"variable"])),]
nrow(water_dna_unique_otus)
#457

water_cdna_otus<-subset(final_table_wholetax_mean,habitat=="water"&nucleic_acid=="cdna"&value!=0)
water_cdna_unique_otus<-water_cdna_otus[which(!duplicated(water_cdna_otus[,"variable"])),]
nrow(water_cdna_unique_otus)
#361

length(intersect(water_dna_unique_otus$variable,water_cdna_unique_otus$variable))
#377

#####subsetting 4 groups:
water_control_otus<-subset(final_table_wholetax_mean,treatment=="control"&habitat=="water"&value!=0)
water_control_unique_otus<-water_control_otus[which(!duplicated(water_control_otus[,"variable"])),]
nrow(water_control_unique_otus)
#442

biofilm_control_otus<-subset(final_table_wholetax_mean,treatment=="control"&habitat=="biofilm"&value!=0)
biofilm_control_unique_otus<-biofilm_control_otus[which(!duplicated(biofilm_control_otus[,"variable"])),]
nrow(biofilm_control_unique_otus)
#415

water_glyph_otus<-subset(final_table_wholetax_mean,treatment=="glyph"&habitat=="water"&value!=0)
water_glyph_unique_otus<-water_glyph_otus[which(!duplicated(water_glyph_otus[,"variable"])),]
nrow(water_glyph_unique_otus)
#436

biofilm_glyph_otus<-subset(final_table_wholetax_mean,treatment=="glyph"&habitat=="biofilm"&value!=0)
biofilm_glyph_unique_otus<-biofilm_glyph_otus[which(!duplicated(biofilm_glyph_otus[,"variable"])),]
nrow(biofilm_glyph_unique_otus)
#406
######################################################################

length(intersect(water_control_unique_otus$variable,water_glyph_unique_otus$variable))
#342


####pair venn plot
grid.newpage()
venn.plot <- draw.pairwise.venn(area1        = 457,
                                area2        = 361,
                                cross.area   = 282,
                                scaled       = F,
                                category     = c("water_DNA", "water_cDNA"),
                                fill         = c("blue", "red"),
                                alpha        = 0.3,
                                lty          = "blank",
                                cex          = 2,
                                cat.cex      = 2,
                                cat.pos      = c(285, 105),
                                cat.dist     = 0.09,
                                cat.just     = list(c(-1, -1), c(1, 1)),
                                ext.pos      = 30,
                                ext.dist     = -0.05,
                                ext.length   = 0.85,
                                ext.line.lwd = 2,
                                ext.line.lty = "dashed")
grid.draw(venn.plot)
dev.copy(png,"Venn_water_nucleic_acids.png")
dev.off()
###triple venn plot
threeway.Venn(water_control_unique_otus$variable,water_glyph_unique_otus$variable,biofilm_glyph_unique_otus$variable)
dev.copy(png,"plotname.png")
dev.off()
###quadruple venn plot
fourway.Venn(water_glyph_unique_otus$variable,water_control_unique_otus$variable,biofilm_glyph_unique_otus$variable,biofilm_control_unique_otus$variable)
dev.copy(png,"plotname.png")
dev.off()

####wie intersect otus etc rausziehen?
is.element(x, y)?