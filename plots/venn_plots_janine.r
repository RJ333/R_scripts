####venn plots janine
##VennDiagram
library(VennDiagramm)
#functions for plotting and calculating
outersect <- function(x, y, ...) {
  big.vec <- c(x, y, ...)
  duplicates <- big.vec[duplicated(big.vec)]
  setdiff(big.vec, unique(duplicates))
}
threeway.Venn <- function(A,B,C,cat.names = c("gut", "filter25", "ganz")){
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

##befehle für schnittmengen
setdiff(gut_otus$variable,filter25_otus$variable)
length(setdiff(gut_otus$variable,filter25_otus$variable))
length(setdiff(filter25_otus$variable,gut_otus$variable))

#multiple intersects:
Reduce(intersect, list(a,b,c,d))

##VennDiagram needs to know how many entries each group and the overlapping sections contain
janines_final<-read.csv(file.choose(),row.names=1,sep=";") #daten einlesen

janines_final$depth<-as.numeric(as.character(janines_final$depth)) # je nach suchkriterium faktor oder numeric nutzen, mit str() checken


#######################################################################
#####subsetting 2 groups:
water_otus<-subset(janines_final,habitat=="water"&value!=0)
water_unique_otus<-water_otus[which(!duplicated(water_otus[,"variable"])),]
nrow(water_unique_otus)
#536

biofilm_otus<-subset(janines_final,habitat=="biofilm"&value!=0)
biofilm_unique_otus<-biofilm_otus[which(!duplicated(biofilm_otus[,"variable"])),]
nrow(biofilm_unique_otus)
#494

length(intersect(water_unique_otus$variable,biofilm_unique_otus$variable))
#377
####################################################################


#####subsetting 3 groups:
gut_otus<-subset(janines_final,material=="gut"&value!=0)
gut_unique_otus<-gut_otus[which(!duplicated(gut_otus[,"variable"])),]
nrow(gut_unique_otus)
#442

filter25_otus<-subset(janines_final,material == "filter" & depth >= 25 & depth <= 45 &extraction=="REPLI" & value!=0)
filter25_unique_otus<-filter25_otus[which(!duplicated(filter25_otus[,"variable"])),]
nrow(filter25_unique_otus)
#415

ganz_otus<-subset(janines_final,material=="ganzA" | material== "ganzT"&value!=0)
ganz_unique_otus<-ganz_otus[which(!duplicated(ganz_otus[,"variable"])),]
nrow(ganz_unique_otus)
#436

######################################################################
#plotting
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
A<-gut_otus$variable
B<-filter25_otus$variable
C<-ganz_otus$variable
threeway.Venn(A=C,B=B,C=A)
dev.copy(png,"janine1.png")
dev.off()

###quadruple venn plot
fourway.Venn(water_glyph_unique_otus$variable,water_control_unique_otus$variable,biofilm_glyph_unique_otus$variable,biofilm_control_unique_otus$variable)
dev.copy(png,"plotname.png")
dev.off()